"""One source-bound positive-J basal correction below a frozen tongue body.

The shared M/APO/SL/TV carriers remain literal. One C1 vector field acts on the
entire12-part dependency closure below their protected plane; no fixed bone is
moved. This construction never installs and still needs actual geometry audits.
"""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
CORE={'TA98TONGUE:male:mucosa','TA98TONGUE:male:aponeurosis','TA98TONGUE:male:intrinsic:superior-longitudinal','TA98TONGUE:male:intrinsic:transverse-vertical-shared'}
def read(p):return json.loads(Path(p).read_text())
def require(x,message):
    if not x:raise ValueError(message)
def field(v,low,high,delta):
    u=np.clip((v[:,1]-low)/(high-low),0.,1.);w=1-u*u*(3-2*u);dy=-6*u*(1-u)/(high-low);out=v+w[:,None]*delta
    j=np.broadcast_to(np.eye(3),(len(v),3,3)).copy();j[:,:,1]+=dy[:,None]*delta[None,:];return out,j,w
def main(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use a new protected basal stage');fit_path=W/'female-axis-tongue-fit-v1/fit.json';prior=read(fit_path)
    for path,sha in prior['source_bindings'].items():require(digest_file(path)==sha,'Frozen source/target changed: '+path)
    origin_path=W/'female-common-tongue-fit-v1/fit.json';origin=read(origin_path);h=next(r for r in origin['landmarks']if r['name']=='hyoid-body-anterior-midline');matrix=np.asarray(prior['matrix']);sp=np.array(h['source']['point'])@matrix[:3,:3].T+matrix[:3,3];tp=np.array(h['target']['point']);delta=tp-sp
    loaded={r['partId']:Mesh.load(r['mesh'])for r in prior['parts']};low=float(sp[1]);high=min(float(loaded[pid].vertices[:,1].min())for pid in CORE);bound=1-1.5*delta[1]/(high-low)
    require(delta[1]>0 and high>low and bound>0,'Source-bound basal field is not globally orientation preserving')
    contract=dict(schema='female-tongue-protected-basal-contract/v1',source_fit=str(fit_path),source_fit_sha256=digest_file(fit_path),source_hyoid_witness=h['source'],target_hyoid_witness=h['target'],transformed_source_hyoid_point=sp.tolist(),target_hyoid_point=tp.tolist(),delta_m=delta.tolist(),
        field='p + delta * (1 - smoothstep((Y-low)/(high-low))); clamped outside [low,high]',low_Y_m=low,protected_Y_m=high,analytic_determinant_lower_bound=float(bound),globally_invertible=True,
        protected_core_carriers=sorted(CORE),fixed_skeleton=True,shared_positions_and_normals_use_one_identical_field=True,source_bound_target_not_parameter_search=True,
        remaining=['Actual discrete self/orientation checks','Female fixed bone/soft context checks','Explicit local extrinsic attachment correction where measured gaps/crossings remain'],installed=False,anatomical_acceptance=False)
    output.mkdir(parents=True);atomic_json(output/'contract.json',contract);rows=[];preservation=[]
    for i,row in enumerate(prior['parts']):
        source=loaded[row['partId']];v,j,w=field(source.vertices,low,high,delta);changed=w>0;attrs=dict(source.attributes)
        if 'normals'in attrs and np.any(changed):
            raw=attrs['normals'].copy();n=np.linalg.solve(j[changed].transpose(0,2,1),raw[changed].astype(float)[...,None])[...,0];n/=np.linalg.norm(n,axis=1)[:,None];raw[changed]=np.rint(n*32767).astype('<i2');attrs['normals']=raw
        m=Mesh(v,source.faces.copy(),source.frame,attrs,dict(source.source,protected_basal_contract=str(output/'contract.json')));path=output/('candidate-'+str(i).zfill(2)+'.npz');m.save(path)
        if row['partId']in CORE:
            literal=np.array_equal(m.vertices,source.vertices)and np.array_equal(m.faces,source.faces)and all(np.array_equal(m.attributes[k],source.attributes[k])for k in source.attributes);require(literal,'Protected core geometry or raw chart words changed');preservation.append(dict(partId=row['partId'],literal=True,geometry_sha256=m.geometry_sha))
        rows.append(dict(row,mesh=str(path),sha256=digest_file(path),geometry_sha256=m.geometry_sha,pre_basal_mesh=row['mesh'],pre_basal_sha256=row['sha256'],basal_controls_changed=int(changed.sum()),maximum_displacement_m=float(np.linalg.norm(v-source.vertices,axis=1).max()),minimum_sampled_jacobian=float(np.linalg.det(j).min())))
    report=dict(prior,schema='installed-tongue-female-protected-basal-fit/v1',parts=rows,body_placement_review_pending=False,body_placement='Selected by actual root Blender review48993c; four core carrier geometries preserved literally',
        basal_contract=str(output/'contract.json'),basal_contract_sha256=digest_file(output/'contract.json'),local_corrections=[contract],protected_core_preservation=preservation,
        source_bindings={**prior['source_bindings'],str(fit_path):digest_file(fit_path),str(origin_path):digest_file(origin_path),str(Path(__file__).resolve()):digest_file(__file__)},one_proper_uniform_transform=False,
        continuous_cage_acceptance_is_not_discrete_anatomical_acceptance=True,qualifying_candidate=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'fit.json',report);print(json.dumps(dict(output=str(output),contract=contract,part_changes=[dict(partId=r['partId'],changed=r['basal_controls_changed'],max_mm=r['maximum_displacement_m']*1000,Jmin=r['minimum_sampled_jacobian'])for r in rows],elapsed_seconds=report['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
