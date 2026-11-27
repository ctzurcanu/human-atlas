"""Second saved stage: preserve section frames, test whole volume, pack once."""
from pathlib import Path
import sys,time,json,importlib.util,types,copy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.tube_sections import bend_diagnostics
from geometry_tools.algorithms.tube_prism_constraints import StationPrismConstraints

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1';OUT=B/'artery-source-profile-transport-v42'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345';FRAME='model:local-ta98:assembled-metres'

def load(name):
    sp=importlib.util.spec_from_file_location(name.replace('-','_'),ROOT/'server'/name);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def rotation(a,b,weight):
    cross=np.cross(a,b);s=np.linalg.norm(cross);c=float(a@b)
    if s<1e-12:return np.eye(3)
    axis=cross/s;angle=np.arctan2(s,c)*weight;x,y,z=axis;K=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    return np.eye(3)+np.sin(angle)*K+(1-np.cos(angle))*(K@K)

def preserve_words(positions,surfaces,words,origin):
    """Literal corner fields on every unchanged actual parent/branch triangle."""
    a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();lookup={}
    for pid in ['TA98HEPATIC:male:artery','TA98CYSTICBLOOD:male:artery']:
        m=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT);normal=m.attributes['normals'][m.faces].astype(np.int64)*m.attributes['normalSigns'][m.faces,None]
        for triangle,n in zip(m.triangles,normal):
            keys=[tuple(map(tuple,np.roll(triangle,k,axis=0)))for k in range(3)];k=min(range(3),key=keys.__getitem__);lookup[keys[k]]=np.roll(n,k,axis=0)
    encoded=(positions-origin).astype('<f4').astype(float)+origin;replayed=0
    for sid,f in surfaces.items():
        for j,t in enumerate(encoded[f]):
            keys=[tuple(map(tuple,np.roll(t,k,axis=0)))for k in range(3)];k=min(range(3),key=keys.__getitem__)
            if keys[k]in lookup:words[sid][j]=np.roll(lookup[keys[k]],-k,axis=0).astype('<i2');replayed+=1
    atomic_json(OUT/'literal-actual-corner-replay.json',dict(sourceCatalogueSha256=EXPECTED,literalActualTrianglesWithOriginalCornerWords=replayed,coordinateMatchingIsExact=True))

def partition_and_pack():
    text=(ROOT/'server/partition-cystic-artery-local-stage3.py').read_text().replace('artery-exterior-corridor-v37','artery-source-profile-transport-v42')
    mod=types.ModuleType('partition_v42');mod.__file__=str(ROOT/'server/partition-cystic-artery-local-stage3.py');exec(compile(text,mod.__file__,'exec'),mod.__dict__);mod.main()
    text=(ROOT/'server/pack-cystic-artery-local-stage3.py').read_text().replace('artery-exterior-corridor-v37','artery-source-profile-transport-v42')
    text=text.replace('packed=pack_finite_cells(positions,surfaces,words,cells,', 'preserve_words(positions,surfaces,words,origin)\n    packed=pack_finite_cells(positions,surfaces,words,cells,')
    mod=types.ModuleType('pack_v42');mod.__file__=str(ROOT/'server/pack-cystic-artery-local-stage3.py');mod.preserve_words=preserve_words;exec(compile(text,mod.__file__,'exec'),mod.__dict__);sys.argv=['pack','artery'];mod.main()
    folder=OUT/'bank/artery';generated=json.loads((folder/'records.json').read_text());a=json.loads(CAT.read_text());current={p['id']:p for p in a['parts']};geometry={'positions','normals','normalSigns','indices','indexCount','vertexCount','geometryOrigin','triangles','bounds','tissueDrawGroups','pairedTissue','sectionSolid'}
    records={}
    for pid,p in generated.items():
        q=copy.deepcopy(current[pid])
        for k in ['geometryScale','geometryInstance','sourceOffset','uvs','colors']:q.pop(k,None)
        for k in geometry:
            if k in p:q[k]=p[k]
        records[pid]=q
    atomic_json(folder/'records.json',records,overwrite=True)
    atomic_json(OUT/'metadata-preservation.json',dict(sourceCatalogueSha256=EXPECTED,allTenSourceRecordIds=list(records),currentSourceMetadataRetained=True,geometryFieldsOnlyReplaced=sorted(geometry),termCountUnchanged=True))

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED;assert not OUT.exists();OUT.mkdir();v14=load('construct-admissible-cystic-departure-v14.py');v10=v14.v10
    a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();source=np.load(B/'admissible-departure-v14/sections.npz');old_outer=source['outer_sections'];old_fluid=source['fluid_sections'];old_centers=source['centers'];N=len(old_centers);variable=np.arange(8,22);fixed=np.setdiff1d(np.arange(N),variable)
    candidate=np.load(B/'artery-valid-envelope-corridor-v41/sections.npz');centers=candidate['centers'];delta=centers-old_centers;assert np.all(delta[fixed]==0)
    oldt=np.gradient(old_centers,axis=0);oldt/=np.linalg.norm(oldt,axis=1)[:,None];newt=np.gradient(centers,axis=0);newt/=np.linalg.norm(newt,axis=1)[:,None]
    contract=dict(schema='cystic-artery-source-profile-transport/v42',sourceCatalogueSha256=EXPECTED,sourceSectionsSha256=digest_file(B/'admissible-departure-v14/sections.npz'),centerCorridorSourceSha256=digest_file(B/'artery-valid-envelope-corridor-v41/sections.npz'),movableRings=variable.tolist(),fixedRings=fixed.tolist(),first8RootSectionsAndRing7AnnulusLiteral=True,terminalGBEntryLiteral=True,sourceShapesAndSizesRetainedByProperTransform=True,actualGlobalDuctOracle=str(B/'biliary-actual-envelope-v39/proof.json'),fixedVeinBankSha256='fb7943d2f5c20f4e7848a386649ef1b646ad725da9eadc82d54a1173ad3392ef',allOtherOrgansAndNativeParentFarFixed=True,stage=2,maximumSavedStages=3,installed=False)
    atomic_json(OUT/'assembly-contract.json',contract)
    contexts=[(pid,reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT))for pid in ['TA98CYSTICBLOOD:male:vein','TA98HEPATIC:male:portal','TA98HEPATIC:male:hepatic_vein','TA98HEPATIC:male:biliary','FJ3080','TA98PART:FJ3079:A05.8.01.061']]+[('GB-innermost',Mesh.load(B/'first-similarity/fixed-mucosa-free-surface.npz'))]
    oracle=Mesh.load(B/'biliary-actual-envelope-v39/actual-biliary-outer-computational-envelope.npz');oracleproof=json.loads((B/'biliary-actual-envelope-v39/proof.json').read_text());assert oracleproof['passed']and next(m for p,m in contexts if p=='TA98HEPATIC:male:biliary').geometry_sha==oracleproof['actualWallGeometrySha256']
    root=old_centers[0];ds=make_solid((oracle.vertices-root)*1000,oracle.faces)[0];host=Mesh.load(B/'fixed-target-20.npz');blood=Mesh.load(B/'native-parent-envelopes-v9/artery-blood.npz');lineage=json.loads((B/'native-parent-envelopes-v9/artery-inner-lineage.json').read_text());_,ho,hb,residual=v10.native_outer(host,blood,lineage,root)
    weights=[]
    for alpha in [0.,.25,.5,.75]:
        w=np.zeros(N);w[variable]=alpha;weights.append((f'constant-source-frame-weight-{alpha}',w))
    taper=np.zeros(N);u=np.clip((np.arange(N)-11)/4,0,1);taper[variable]=(u*u*(3-2*u))[variable];weights.insert(0,('source-root-frame-smoothly-transported-after-ring11',taper))
    trials=[];accepted=None
    for label,weight in weights:
        if time.monotonic()-began>27:break
        outer=old_outer.copy();fluid=old_fluid.copy()
        for j in variable:
            R=rotation(oldt[j],newt[j],weight[j]);outer[j]=(old_outer[j]-old_centers[j])@R.T+centers[j];fluid[j]=(old_fluid[j]-old_centers[j])@R.T+centers[j]
        assert np.array_equal(outer[fixed],old_outer[fixed])and np.array_equal(fluid[fixed],old_fluid[fixed])
        ov,of,_=v10.swept_profiles(outer);iv,iff,_=v10.swept_profiles(fluid);side=of[:-36];fn=np.cross(ov[side[:,1]]-ov[side[:,0]],ov[side[:,2]]-ov[side[:,0]]);radial=ov[side].mean(1)-np.repeat(centers,36,axis=0)[:len(side)]
        if np.median(np.sum(fn*radial,axis=1))<0:side=side[:,::-1]
        nominal=centers[:,None]+.82*(outer-centers[:,None]);prism=StationPrismConstraints(np.array([outer.reshape(-1,3),nominal.reshape(-1,3)]),side,np.repeat(np.arange(N),18),np.ones(N,bool),allow_source_faults=True);bad=np.flatnonzero(~prism.sourceHealthy[0]);invalid=np.unique(np.repeat(np.arange(N),18)[side[bad]].min(1)).tolist();trial=dict(frameTransport=label,invalidContinuousPrismIntervals=invalid)
        trials.append(trial)
        if invalid:continue
        bo=v10.previous.solid(ov,of,root);bb=v10.previous.solid(iv,iff,root);occupation=float((bo^ds).volume());trial['fullContinuousOuterEnvelopeDuctOccupation_mm3']=occupation
        if occupation!=0:continue
        whole=(ho+bo)-(hb+bb);branchsolid=(bo-hb)-bb;network=Mesh(*v10.previous.arrays(whole,root),FRAME);branch=Mesh(*v10.previous.arrays(branchsolid,root),FRAME);top=topology(network);bloodtop=topology(Mesh(*v10.previous.arrays(hb+bb,root),FRAME));escape=[]
        for i,s in enumerate((bb-bo).decompose()):
            p,f=v10.previous.arrays(s,root);near=int(np.linalg.norm(centers-p.mean(0),axis=1).argmin());volume=float(s.volume())/1e9
            if near<N-2 and volume>1e-23:escape.append(dict(component=i,nearestStation=near,volume_m3=volume))
        checks=[]
        for pid,m in contexts:
            data,c=pair_contacts(branch,m,0.,max_pairs=500000,max_tests=5000000);assert c['complete_within_pair_budget'];checks.append(dict(id=pid,geometrySha256=m.geometry_sha,**c))
        trial.update(contextCrossings={c['id']:c['crossing_face_pairs']for c in checks},lateralBloodEscape=escape)
        okay=not escape and all(c['crossing_face_pairs']==c['coplanar_overlap_face_pairs']==0 for c in checks)and top['closed_oriented_manifold']and len(top['components'])==1 and top['components'][0]['closed_orientable_genus']==1 and len((hb^bb).decompose())==1 and len(bloodtop['components'])==1 and bloodtop['components'][0]['closed_orientable_genus']==0
        if okay:accepted=(outer,fluid,network,branch,top,bloodtop,checks,trial);break
    atomic_json(OUT/'bounded-frame-transport-trials.json',dict(trials=trials,elapsedSeconds=time.monotonic()-began,installed=False))
    if accepted is None:
        atomic_json(OUT/'gates.json',dict(allRequiredGeometryGatesPassed=False,stoppedBeforePartition=True,trials=trials,installed=False,elapsedSeconds=time.monotonic()-began));print(json.dumps(dict(passed=False,trials=trials,elapsedSeconds=time.monotonic()-began)));return
    outer,fluid,network,branch,top,bloodtop,checks,trial=accepted;network.save(OUT/'artery-joined-wall.npz');branch.save(OUT/'artery-branch-wall.npz');atomic_npz(OUT/'sections.npz',outer_sections=outer,fluid_sections=fluid,centers=centers,source_outer_sections=old_outer,source_fluid_sections=old_fluid,fixed_stations=fixed,station_displacements_m=delta)
    bend=bend_diagnostics(centers,outer[None]);g=dict(allRequiredGeometryGatesPassed=True,sourceCatalogueSha256=EXPECTED,frameTransport=trial['frameTransport'],completeOuterDuctOccupation_mm3=0.,continuousPrismInvalidIntervals=[],joinedTopology=top,composedBloodTopology=bloodtop,contextChecks=checks,sourceSectionShapeAndCalibrePreserved=True,first8RootSectionsLiteral=True,terminalAnnulusAndGBEntryLiteral=True,maximumStationTurnDegrees=float(np.rad2deg(bend['station_turn_rad']).max()),maximumCurvatureSupport=float(bend['curvature_support'].max()),assemblyContract=str(OUT/'assembly-contract.json'),installed=False,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'gates.json',g)
    if time.monotonic()-began>31:atomic_json(OUT/'pack-held.json',dict(reason='Reserved stage budget insufficient for pack',elapsedSeconds=time.monotonic()-began));return
    partition_and_pack();assert digest_file(CAT)==EXPECTED;print(json.dumps(dict(passed=True,bank=str(OUT/'bank/artery/cells-and-regions.bin'),frameTransport=trial['frameTransport'],elapsedSeconds=time.monotonic()-began)))
if __name__=='__main__':main()
