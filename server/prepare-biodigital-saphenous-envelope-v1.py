"""Measured medial-thigh vein envelope candidate; no geometry installation."""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.oriented_rays import OrientedRayIndex
from geometry_tools.algorithms.compact_positive_fields import compose_compact_radial_maps

started=time.monotonic();B=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
out=B/'saphenous-envelope-v1';out.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
base=B/'lower-limb-joints-v4';candidate=base/'candidate-v1'
contract=json.loads((base/'contract.json').read_text())
rows=json.loads((candidate/'instance-ownership.json').read_text())['instances']
source=np.load(candidate/'venous.npz');vertices=source['vertices'];faces=source['faces']
probe_path=base/'localized-skin-regression-v1/moved-vertex-classification.npz'
probe=np.load(probe_path);n=len(probe['after_points']);outside=probe['classification'][n:]==0
controls=probe['after_points'][outside];assert len(controls)==35
fields=[];section_rows=[];flatcore_rows=[]
for side in ['right','left']:
    sid=f'human_20_male_cardiovascular_system-{side}_great_saphenous_vein_ID'
    row=next(r for r in rows if r['sourceInstanceId']==sid)
    f=faces[row['candidateFaceStart']:row['candidateFaceStart']+row['candidateFaceCount']]
    tri=vertices[f];stations=[]
    for y in [.555,.58,.60,.625,.65,.675,.70,.715]:
        cut=[]
        for a,b in [(0,1),(1,2),(2,0)]:
            hit=(tri[:,a,1]<y)&(tri[:,b,1]>=y)|(tri[:,b,1]<y)&(tri[:,a,1]>=y)
            t=tri[hit];u=(y-t[:,a,1])/(t[:,b,1]-t[:,a,1]);cut.append(t[:,a]+u[:,None]*(t[:,b]-t[:,a]))
        points=np.concatenate(cut);assert len(points)>=6,(sid,y)
        center=(points.min(0)+points.max(0))/2
        stations.append(dict(y_m=y,x_m=float(center[0]),z_m=float(center[2])))
        section_rows.append(dict(sourceId=sid,Y_m=y,controls=len(points),bounds_m=[points.min(0).tolist(),points.max(0).tolist()]))
    sign=np.sign(np.mean([s['x_m']for s in stations]));assert sign!=0
    field=dict(kind='separable-corridor-translation',centerline_stations=stations,center_y_m=.633,
        core_axes_m=[.008,.024,.008],support_axes_m=[.038,.085,.038],
        y_core_interval_m=[.605,.660],y_support_interval_m=[.555,.715],
        translation_m=[float(sign*.011),0.,0.])
    # All literal main-vein controls in the full-translation Y core must stay
    # within its transverse core, rather than being squeezed by the fade.
    from geometry_tools.algorithms.compact_positive_fields import _corridor_charts
    px,_,pz,_=_corridor_charts(tuple((s['y_m'],s['x_m'],s['z_m'])for s in stations))
    v=vertices[row['candidateVertexStart']:row['candidateVertexStart']+row['candidateVertexCount']]
    core=v[(v[:,1]>=.605)&(v[:,1]<=.660)]
    offsets=np.abs(np.c_[core[:,0]-px(core[:,1]),core[:,2]-pz(core[:,1])])
    assert np.all(offsets<np.array([.008,.008])),(sid,offsets.max(0))
    flatcore_rows.append(dict(sourceId=sid,vertices=len(core),maximum_XZ_offset_m=offsets.max(0).tolist()))
    fields.append(field)
skin_path=B/'global-envelope-v1/fixed-skin/skin.npz'
skin=Mesh.load(skin_path);boundary_path=B/'global-envelope-v1/fixed-skin/envelope-topology-v1/boundary-uncertainty-v1/uncovered-boundary-segments.npz';boundary=np.load(boundary_path)
index=OrientedRayIndex(skin,nthreads=2,uncertainty_segments=boundary['primary_segments_m'],uncertainty_component_ids=boundary['primary_boundary_component_index'])
proposed,jac,proof=compose_compact_radial_maps(controls,fields)
classified,inside_proof=index.classify(proposed,probe['directions'])
nearest,nearest_proof=SurfaceIndex(skin,nthreads=2).closest_points(proposed,certify=True)
ready=bool(np.all(classified['classification']==1)and nearest['distance_m'].min()>=.0025)
np.savez_compressed(out/'controls.npz',before=controls,proposed=proposed,jacobian=jac,
    source_instance_indices=probe['source_instance_indices'][outside],source_vertex_ids=probe['original_source_vertex_ids'][outside],
    nearest_skin_points=nearest['points'],nearest_skin_faces=nearest['face_ids'],certified_nearest_distance_m=nearest['distance_m'],classification=classified['classification'])
(out/'field-recipe.json').write_text(json.dumps(dict(fields=fields,sections=section_rows,flatcore=flatcore_rows),indent=2)+'\n')
report=dict(schema='measured-medial-thigh-saphenous-envelope/v1',baseline=str(candidate),
    inputs={str(p):sha(p)for p in [base/'contract.json',candidate/'venous.npz',candidate/'instance-ownership.json',probe_path,skin_path,boundary_path]},
    sourceIds=[r['sourceId']for r in flatcore_rows],controls=len(controls),readyForOneTrial=ready,
    proposedInside=int((classified['classification']==1).sum()),proposedOutside=int((classified['classification']==0).sum()),
    proposedUncertain=int((classified['classification']<0).sum()),minimumCertifiedSkinDistance_m=float(nearest['distance_m'].min()),
    maximumMovement_m=float(np.linalg.norm(proposed-controls,axis=1).max()),sourceMainCorePreserved=True,
    outsideSourceFaceAndAttachedBranchChecksPending=True,sourceGeometryModified=False,installed=False,
    rayProof=inside_proof,nearestProof=nearest_proof,elapsedSeconds=time.monotonic()-started)
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
if ready:
    contract['post_fields']+=fields
    _,_,all_proof=compose_compact_radial_maps(np.zeros((1,3)),contract['post_fields'])
    contract['analytic_proof']['analytic_determinant_lower_bound']=.19428726760129453*all_proof['analytic_determinant_lower_bound']
    contract['stagePurpose']='Resolve measured medial-thigh great-saphenous skin protrusions; preserve deep knee correction and all source members'
    contract['regionalEnvelopeEvidence']={'report':'saphenous-envelope-v1/report.json','report_sha256':sha(out/'report.json'),'skin_sha256':sha(skin_path)}
    (out/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
print(json.dumps({k:report[k]for k in ['readyForOneTrial','controls','proposedInside','proposedOutside','proposedUncertain','minimumCertifiedSkinDistance_m','maximumMovement_m','elapsedSeconds']}))
