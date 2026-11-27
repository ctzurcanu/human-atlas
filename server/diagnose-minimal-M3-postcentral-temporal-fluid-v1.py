"""Root-passed three-course topology only; secondary leaves omitted diagnostically."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.audit import audit,topology
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance
import manifold3d as md
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';C=P/'three-course-finite-educational-v1';I=P/'temporal-source-incidence-v2';FIVE=I/'five-course-fluid-diagnostic-v1';OUT=I/'minimal-three-course-fluid-diagnostic-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);g=json.loads((C/'generation-complete.json').read_text());origin=np.array(g['origin_m']);sources=[]
 for role in ['M3-fluid','main-fluid']:
  item=next(r for r in g['native_inputs']if r['role']==role);meta=next(r for r in g['classifiers']if r['role']==role);sources.append((role,C/item['checkpoint_input'],meta))
 temporal_input=next(FIVE.glob('*-make-temporal0-fluid-input.npz'));temporal=Mesh.load(FIVE/'temporal0-diagnostic-fluid.npz');tz=np.load(I/'complete-source-charts.npz');distal=tz['temporal_vertices'][tz['component0_rings'][168]].mean(0);degree=np.bincount(temporal.faces.ravel(),minlength=len(temporal.vertices));tips=np.flatnonzero((np.arange(len(degree))>=len(np.unique(tz['temporal_faces'][tz['component0_source_face_ids']])))&(degree==12));tip=int(tips[np.argmin(np.linalg.norm(temporal.vertices[tips]-distal,axis=1))]);assert np.linalg.norm(temporal.vertices[tip]-distal)<2e-6;diskids=np.flatnonzero(np.any(temporal.faces==tip,axis=1)).tolist();assert len(diskids)==12;temporal_meta=dict(source_faces=4032,mouths=[dict(kind='external-fluid-classifier',mouth=0,newFaceIds=diskids)]);sources.append(('temporal-long-fluid',temporal_input,temporal_meta));solids=[];records=[];op=0
 def checkpoint(name,arrays=None):
  nonlocal op
  op+=1;stem=f'{op:02d}-{name}'
  if arrays is not None:atomic_npz(OUT/(stem+'-input.npz'),**arrays)
  atomic_json(OUT/(stem+'-pending.json'),dict(operation=name,elapsed_seconds=time.monotonic()-start,input_archive=stem+'-input.npz'if arrays is not None else None));return stem
 for role,path,meta in sources:
  q=np.load(path);stem=checkpoint('make-'+role,dict(positions_mm=q['positions_mm'],faces=q['faces']));solid,oid=make_solid(q['positions_mm'],q['faces']);solids.append(solid);records.append(dict(role=role,original_id=oid,source_checkpoint=str(path),source_side_faces=meta['source_faces'],mouths=meta['mouths']));atomic_json(OUT/(stem+'-done.json'),dict(status=str(solid.status()),original_id=oid))
 stem=checkpoint('union-root-passed-three-fluid-courses');joined=md.Manifold.batch_boolean(solids,md.OpType.Add);assert joined.status()==md.Error.NoError;checkpoint('extract-topology-only-union');out=joined.to_mesh64();ids,fids=provenance(out);vertices=np.asarray(out.vert_properties)[:,:3]/1000.+origin;faces=np.asarray(out.tri_verts,dtype=np.int64);mesh=Mesh(vertices,faces,FRAME);mesh.save(OUT/'diagnostic-fluid-union.npz');atomic_npz(OUT/'diagnostic-provenance.npz',original_ids=ids,source_face_ids=fids,frame=np.array(FRAME));quality=audit(mesh,require_closed=True,expected_components=1);atomic_json(OUT/'union-intrinsic-audit.json',quality);t=quality['topology'];tree=quality['passed']and len(t['components'])==1 and t['components'][0]['closed_orientable_genus']==0
 # Radii are read from the existing unmodified charts, not changed to make
 # the failed short-course fan fit.
 z=np.load(I/'complete-source-charts.npz');radii=[]
 for name,key,station in [('temporal-long-parent','component0_rings',144),('temporal-short-child','component1_rings',48)]:
  ring=z['temporal_vertices'][z[key][station]];center=ring.mean(0);R=float(np.linalg.norm(ring-center,axis=1).mean());thickness=min(.00012*np.sqrt(R/.001),.35*R);radii.append(dict(role=name,station=station,outer_radius_m=R,derived_educational_thickness_m=thickness,inner_radius_m=R-thickness,formula='min(0.12mm*sqrt(R/1mm),0.35R)',native_radius_variation_m=float(np.ptp(np.linalg.norm(ring-center,axis=1)))))
 mouths=[]
 for row in records:
  for mouth in row['mouths']:
   if mouth['kind']!='external-fluid-classifier':continue
   count=int(np.sum((ids==row['original_id'])&np.isin(fids,mouth['newFaceIds'])));mouths.append(dict(role=row['role'],mouth=mouth['mouth'],surviving_classifier_disk_faces=count))
 tree=tree and len(mouths)==4 and all(r['surviving_classifier_disk_faces']>0 for r in mouths)
 result=dict(schema='minimal-root-passed-M3-postcentral-temporal-fluid-topology/v1',diagnostic_only=True,installable_candidate=False,no_wall_or_fit_constructed=True,omitted_only_for_this_diagnostic=['secondary postcentral source course','short temporal source course'],complete_five_course_failure_unchanged=True,root_inputs=records,expected_terminal_classifier_mouths=4,actual_terminal_classifier_mouths=mouths,expected_tree_genus=0,actual_topology=t,self_audit_passed=quality['passed'],extra_handle_count=t['components'][0]['closed_orientable_genus']if len(t['components'])==1 else None,tree_topology_test_passed=bool(tree),short_temporal_failed_root_radii=radii,anatomical_acceptance=False,scope='Topology answers whether the actual long-temporal/main-postcentral overlap creates a handle. Genus0 alone cannot authorize outside-parent tissue merging, cortex acceptance or installation.',elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'topology-proof.json',result);print(json.dumps(dict(tree_topology_test_passed=tree,actual_components=t['components'],self_audit_passed=quality['passed'],actual_terminal_classifier_mouths=mouths,short_temporal_failed_root_radii=radii,elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
