"""Check emitted lung cells, airway connectivity and named territories."""
from pathlib import Path
import argparse,json,sys,time,itertools
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import topology
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.shared_cell_disjoint import shared_cell_disjoint
from geometry_tools.open3d_backend import library
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--catalogue',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--allow-apical-tissue-change',action='store_true');args=ap.parse_args();start=time.monotonic()
path=ROOT/args.catalogue;before_path=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{args.sex}/atlas.json';a=json.loads(path.read_text());before=json.loads(before_path.read_text());parts={p['id']:p for p in a['parts']};bp={p['id']:p for p in before['parts']};reader=AtlasReader();frame=f'model:local-{args.sex}-ta98:assembled-metres';out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
def read(pid,new=True):return reader.read_part(a if new else before,parts[pid]if new else bp[pid],frame=frame,catalogue_path=path if new else before_path)
ids=a['pulmonarySegmentIntegration']['partIds'];branches=a['pulmonaryBranchIntegration']['partIds'];cells={i:read(i)for i in ids}
# Recenter on the common storage origin and use a binary scale. This retains
# every represented coordinate; no rounding, snapping or volume waiver.
origin=np.array([0.,1.3521,0.])if args.sex=='male'else np.zeros(3);binary_scale=8.
solids={i:make_solid((m.vertices-origin)*binary_scale,m.faces)[0]for i,m in cells.items()};pairs=[]
for i,j in itertools.combinations(ids,2):
 overlap=solids[i]^solids[j];row=dict(first=i,second=j,intersection_triangles=int(overlap.num_tri()),intersection_m3=float(overlap.volume())/binary_scale**3)
 if row['intersection_triangles']:
  proof=shared_cell_disjoint(cells[i],cells[j],max_pairs=1000000,max_tests=8000000)
  # Cells meeting only at edges have no paired interface facets. Certify
  # their complete touching-only contacts and all face centroids instead.
  if proof.get('reason')=='no-explicit-opposed-interface':
   from geometry_tools.contact import pair_contacts
   import igl
   _,contacts=pair_contacts(cells[i],cells[j],0.,max_pairs=1000000,max_tests=8000000);winding=[]
   for first,second in [(cells[i],cells[j]),(cells[j],cells[i])]:
    w=igl.winding_number(second.vertices-origin,second.faces,(first.triangles-origin).mean(1));winding.append(dict(faces_checked=len(w),maximum_absolute_winding=float(abs(w).max()),finite=bool(np.isfinite(w).all())))
   proof=dict(method='complete-touching-only-contacts-and-bidirectional-generalized-winding',contacts=contacts,winding=winding,passed=contacts['complete_within_pair_budget']and contacts['crossing_face_pairs']==contacts['coplanar_overlap_face_pairs']==0 and all(r['finite']and r['maximum_absolute_winding']<1e-7 for r in winding))
  row['independent_sheet_or_overlap_certificate']=proof;row['passed']=proof['passed']
 else:row['passed']=True
 pairs.append(row);atomic_json(out/'occupied-pairs.json',dict(catalogue_sha256=digest_file(path),pairs=pairs,complete=len(pairs)==len(ids)*(len(ids)-1)//2),overwrite=True)
 if not row['passed']:print('Cell pair failed',i,j,proof.get('reason'),flush=True)
print('Occupied cell pairs',len(pairs),'failures',sum(not r['passed']for r in pairs),flush=True)
def combined(meshes):
 offset=np.r_[0,np.cumsum([len(m.vertices)for m in meshes])];return Mesh(np.concatenate([m.vertices for m in meshes]),np.concatenate([m.faces+offset[i]for i,m in enumerate(meshes)]),frame)
old_tree=topology(combined([read(i,False)for i in branches]));new_tree=topology(combined([read(i)for i in branches]));assert len(old_tree['components'])==len(new_tree['components'])==1;assert old_tree['open_edges']==new_tree['open_edges'];assert new_tree['zero_area_faces']==0
upper=1.475 if args.sex=='male'else 1.44;fixed_controls=0
for pid in ids+branches:
 old,new=read(pid,False),read(pid);assert np.array_equal(old.faces,new.faces);mask=old.vertices[:,1]>=upper
 if not(args.allow_apical_tissue_change and pid in ids):assert np.array_equal(old.vertices[mask],new.vertices[mask]);fixed_controls+=int(mask.sum())
e=json.loads(Path(before['pulmonaryBranchIntegration']['namedBranchDestinationAudit']).read_text());o3d=library()
def distance(mesh,q):
 origin=mesh.vertices.mean(0);scale=np.ptp(mesh.vertices,axis=0).max();s=o3d.t.geometry.RaycastingScene(nthreads=2);s.add_triangles(o3d.core.Tensor(((mesh.vertices-origin)/scale).astype(np.float32)),o3d.core.Tensor(mesh.faces.astype(np.uint32)));return s.compute_signed_distance(o3d.core.Tensor(((q-origin)/scale).astype(np.float32)),nthreads=2,nsamples=3).numpy().astype(float)*scale
destinations=[]
for row in e['destinations']:
 name=row['branch'].removeprefix('HRA:VH_F_');bid=f'LOCAL:{args.sex}:pulmonary:bronchus:'+name;pid=f'TA98PULMONARY:{args.sex}:'+row['expected_lobe']+':'+row['expected_region'];assert bid in parts and pid in cells,(bid,pid)
 old,new=read(bid,False),read(bid);oq=np.r_[old.vertices,old.triangles.mean(1)];nq=np.r_[new.vertices,new.triangles.mean(1)];source_distance=distance(read(pid,False),oq);mask=source_distance<=1e-7;target_distance=distance(cells[pid],nq[mask]);outside=int((target_distance>1e-6).sum());destinations.append(dict(branch=bid,target=pid,source_in_region_samples=int(mask.sum()),candidate_outside_expected_region_samples=outside,maximum_candidate_outside_m=float(target_distance.max(initial=0)),tolerance_m=1e-6,passed=outside==0))
report=dict(catalogue_sha256=digest_file(path),source_catalogue_sha256=digest_file(before_path),cells=len(ids),occupied_volume_pairs=pairs,airway_components_before=len(old_tree['components']),airway_components_after=len(new_tree['components']),airway_open_edges_preserved=new_tree['open_edges'],upper_source_controls_literal=fixed_controls,named_destinations=destinations,passed=all(r['passed']for r in destinations+pairs),anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start);atomic_json(out/'report.json',report);print(json.dumps({k:v for k,v in report.items()if k not in ['occupied_volume_pairs','named_destinations']},indent=2));print('Destination failures',json.dumps([r for r in destinations if not r['passed']]),flush=True)
assert report['passed'],'Named airway destinations require correction'
