"""One deterministic phase-I solve of the unchanged v9 source constraints.

No pose sweep: a source-start trust-region sequence minimizes the exact
same profile/plane violation, using mm-scaled variables and checked central
Jacobians. Per-vertex inequalities equal v9's minimum-of-twelve constraint.
"""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.tube_sections import redirect_local_sections

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'
SOURCE=BASE/'one-coupled-foraminal-candidate-v9'
OUT=BASE/'one-coupled-numerical-phase1-v11'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def clean(q):
    if isinstance(q,np.ndarray):return q.tolist()
    if isinstance(q,np.generic):return q.item()
    raise TypeError(type(q).__name__)

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False)
    assert digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')==EXPECTED
    z=np.load(SOURCE/'actual-coupled-constraints-before-solve.npz');K=z['metric'];B={k:z['basis_'+k]for k in ['l','r','b']};count=len(K);L=np.linalg.cholesky(K)
    # Normalize every coefficient by its actual maximum centre displacement.
    # Coordinates and residuals are mm, rotation is an equivalent tangential
    # displacement over the actual root support lever, not raw radians.
    response=np.vstack([q@L for q in B.values()]);column=np.abs(response).max(0);assert (column>0).all();P=L/column[None]
    data={}
    for k in ['l','r','b']:
        d=np.load(BASE/(k+'-literal-source-stations.npz')if k!='b'else BASE/'necessary-basilar-source-closure-v5/basilar-literal-source-stations.npz');data[k]=dict(centers=d['centers'],profiles=d['physical_vertices'][d['rings']])
    pivot=data['b']['centers'][0];lever_mm=max(float((np.linalg.norm(d['centers']-pivot,axis=1)*B[k][:,0]).max())for k,d in data.items())*1000;assert lever_mm>0
    raw=json.loads((SOURCE/'actual-coupled-witnesses-before-solve.json').read_bytes());planes=[]
    for q in raw['rows']:
        w=q.get('sourceWitness');n=np.array(w['outwardFaceUnitNormal']if w else q['normal']);p=np.array(w['point_m']if w else q['point_m']);planes.append((q['owner'],q['station'],n,p))
    cache={}
    def evaluate(x):
        if cache.get('x')is not None and np.array_equal(cache['x'],x):return cache['value']
        coef=(P@x[:-3].reshape(count,3))*.001;R=Rotation.from_rotvec(x[-3:]/lever_mm).as_matrix();targets={};rotations={};profiles={}
        for key,bb in B.items():
            d=data[key];old=d['centers'];common=bb[:,0];target=old+bb@coef+common[:,None]*((old-pivot)@(R-np.eye(3)).T)
            if key!='b':target[:171]=old[:171]
            else:target[19:]=old[19:]
            rr=redirect_local_sections(d['profiles'][None],target,dtype=float)[1];relative=Rotation.from_matrix(R@rr.transpose(0,2,1)).as_rotvec();rr=Rotation.from_rotvec(common[:,None]*relative).as_matrix()@rr
            if key!='b':rr[:170]=np.eye(3);rr[203:]=R
            else:rr[:2]=R;rr[20:]=np.eye(3)
            targets[key]=target;rotations[key]=rr;profiles[key]=np.einsum('sij,svj->svi',rr,d['profiles']-old[:,None])+target[:,None]
        vertex=np.array([((profiles[key][j]-p)@n)-.0001 for key,j,n,p in planes])*1000
        value=(vertex.ravel(),targets,rotations,profiles,coef,R);cache.update(x=x.copy(),value=value);return value
    def fun(x):return evaluate(x)[0]
    def jac(x,h):
        out=np.empty((len(planes)*12,len(x)))
        for j in range(len(x)):
            step=np.zeros_like(x);step[j]=h;out[:,j]=(fun(x+step)-fun(x-step))/(2*h)
        return out
    x=np.zeros(3*count+3);c=fun(x);initial=float(max(0,-c.min()));J=jac(x,1e-4);Jhalf=jac(x,5e-5);error=float(abs(J-Jhalf).max());scale=float(abs(Jhalf).max());checked=error<=max(1e-5,scale*2e-5)
    atomic_json(OUT/'derivative-and-scaling-check.json',dict(initialMaximumViolation_mm=initial,mmEquivalentRotationLever_mm=lever_mm,coefficientColumnScales=column,centralStep=1e-4,halfStep=5e-5,maximumJacobianDifference=error,maximumJacobianMagnitude=scale,checked=checked,
        v9MinimumOfTwelveExactlyEquivalentToAllTwelveVertexInequalities=True,sourceConstraintSnapshotSha256=digest_file(SOURCE/'actual-coupled-constraints-before-solve.npz'),sourceWitnessSha256=digest_file(SOURCE/'actual-coupled-witnesses-before-solve.json')))
    if not checked:raise ValueError('The numerical derivative consistency gate failed')
    history=[];trust=.5;status='iteration limit';m=len(x)
    for iteration in range(12):
        if time.monotonic()-start>43:status='internal numerical deadline';break
        c=fun(x);before=float(max(0,-c.min()))
        if before<=1e-6:status='feasible';break
        J=Jhalf if iteration==0 else jac(x,5e-5)
        # Delta, nonnegative violation slack t, absolute-delta auxiliaries.
        A=np.c_[-J,-np.ones(len(c)),np.zeros((len(c),m))];rhs=c
        absA=np.r_[np.c_[np.eye(m),np.zeros(m),-np.eye(m)],np.c_[-np.eye(m),np.zeros(m),-np.eye(m)]]
        A=np.r_[A,absA];rhs=np.r_[rhs,np.zeros(2*m)];cost=np.r_[np.zeros(m),1.,np.full(m,1e-5)]
        lp=linprog(cost,A_ub=A,b_ub=rhs,bounds=[(-trust,trust)]*m+[(0,None)]+[(0,None)]*m,method='highs',options={'time_limit':max(1.,min(4.,43-(time.monotonic()-start))),'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
        if not lp.success:status='linearized phase-I failed: '+lp.message;history.append(dict(iteration=iteration,lpStatus=int(lp.status),message=lp.message));break
        step=lp.x[:m];predicted=max(0.,float(-(c+J@step).min()));accepted=False;trials=[]
        for alpha in [1.,.5,.25,.125,.0625]:
            value=fun(x+alpha*step);after=float(max(0.,-value.min()));trials.append(dict(alpha=alpha,maximumViolation_mm=after))
            if after<before-max(1e-8,before*1e-6):x+=alpha*step;accepted=True;break
        history.append(dict(iteration=iteration,maximumViolationBefore_mm=before,predictedMaximumViolation_mm=predicted,trust_mm_equivalent=trust,stepNorm_mm_equivalent=float(np.linalg.norm(step)),accepted=accepted,trials=trials,elapsedSeconds=time.monotonic()-start))
        atomic_json(OUT/('progress-'+str(iteration)+'.json'),dict(history=history,elapsedSeconds=time.monotonic()-start))
        atomic_npz(OUT/('checkpoint-'+str(iteration)+'.npz'),scaled_variables=x,trust_mm_equivalent=np.array(trust))
        if iteration==0:
            prior=json.loads((BASE/'one-coupled-numerical-phase1-v10/progress.json').read_bytes())['history'][0]
            fields=['maximumViolationBefore_mm','predictedMaximumViolation_mm','trust_mm_equivalent','stepNorm_mm_equivalent']
            assert all(prior[k]==history[0][k]for k in fields)and prior['trials']==history[0]['trials']
            atomic_json(OUT/'unchanged-prefix-replay.json',dict(firstDeterministicLPTraceReplayedExactly=True,sourceCheckpointFailureSha256=digest_file(BASE/'one-coupled-numerical-phase1-v10/failed-stage.json'),anatomicalInputsAndInitialPoseUnchanged=True))
        if not accepted:
            trust*=.25
            if trust<1e-4:status='no useful checked trust-region direction';break
        elif alpha==1. and after<before*.75:trust=min(1.,trust*1.5)
        else:trust=max(1e-4,trust*.6)
    final,targets,rotations,profiles,coef,R=evaluate(x);violation=float(max(0,-final.min()));passed=violation<=1e-6
    if passed:status='feasible'
    # This is the one resulting proper SE(3) pose, evaluated on the actual
    # source profiles and exact original planes, never a surrogate sphere.
    atomic_npz(OUT/'one-proper-result-fields.npz',scaled_variables=x,coefficients_m=coef,root_rotation=R,root_pivot_m=pivot,**{'target_'+k:q for k,q in targets.items()},**{'rotation_'+k:q for k,q in rotations.items()},**{'actual_profiles_'+k:q for k,q in profiles.items()},actual_vertex_plane_residuals_mm=final.reshape(-1,12))
    result=dict(schema='unchanged-v9-single-source-start-numerical-phase1/v10',sourceCatalogueSha256=EXPECTED,sourceConstraintsSha256=digest_file(SOURCE/'actual-coupled-constraints-before-solve.npz'),sourceWitnessesSha256=digest_file(SOURCE/'actual-coupled-witnesses-before-solve.json'),method='checked central Jacobians; mm-scaled source-start trust-region phase-I LP; actual proper SE3 evaluation',history=history,
        status=status,sourceMaximumViolation_mm=initial,finalMaximumViolation_mm=violation,feasibleAtDeclaredOneNanometreNumericalResidual=passed,rootRotationVector_rad=Rotation.from_matrix(R).as_rotvec(),rootProperDeterminant=float(np.linalg.det(R)),rootOrthogonalityError=float(abs(R.T@R-np.eye(3)).max()),rootTranslation_m=coef[0],variablesChanged=bool(np.any(x)),
        anatomicalPlanesAndClearanceAndSourceCalibreUnchanged=True,noRandomInitialPoses=True,noGeometryCandidateOrAtlasChange=True,wholeContextAndSelfGatesStillRequired=True,acceptedForInstallation=False,elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)));print(json.dumps(dict(evidence=str(OUT/'evidence.json'),status=status,initialViolation_mm=initial,finalViolation_mm=violation,iterations=len(history),elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':
    try:main()
    except Exception as e:
        if OUT.exists():atomic_json(OUT/'failed-stage.json',dict(error=type(e).__name__,message=str(e),notInstallable=True))
        raise
