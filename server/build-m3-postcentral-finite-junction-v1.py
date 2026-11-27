"""One bounded, source-bound three-port construction; no catalogue writes.

The inner profiles are an educational construction assumption, NOT measured
media/adventitial thickness or a recovered native finite arterial wall.
"""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology,winding_contains
from geometry_tools.algorithms.branch_junction_envelopes import close_branch
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops,signed_volume
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance
import manifold3d as md

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3'
P=B/'m3-postcentral-port-feasibility-v1'
OUT=P/'finite-junction-educational-v1'
FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);z=np.load(P/'source-port-charts.npz');plan=json.loads((P/'feasibility.json').read_text());source=B/'source-atlas.json'
 assert digest_file(source)==plan['source_catalogue_sha256'];assert digest_file(P/'source-port-charts.npz')==plan['source_chart_archive_sha256']
 atomic_json(OUT/'construction-contract.json',dict(source_catalogue_sha256=digest_file(source),source_charts_sha256=digest_file(P/'source-port-charts.npz'),frame=FRAME,installed=False,scope='Only M3 stations 33..41 and postcentral component 0 stations 176..180. Other postcentral component and its 58 pre-existing crossings remain unchanged/pending.',inner_profiles='Derived circular radial scaling of literal 12-control outer station charts',thickness_formula='t(R)=min(0.12 mm * sqrt(R / 1 mm),0.35 R)',thickness_status='Tunable educational construction assumption; not clinical measurement or a named histological coat.',trend_source='https://pubmed.ncbi.nlm.nih.gov/10503931/',trend_limit='Supports only increasing relative media thickness as human pial arteries narrow; does not supply this formula.',fixed_outer_faces_outside_local_domain=True,fixed_source_attributes_outside_local_domain=True,fixed_brain_and_skeleton=True))
 # Save literal decoded source arrays/attributes, never inferred index pairing.
 a=json.loads(source.read_bytes());parts={p['id']:p for p in a['parts']};reader=AtlasReader();source_rows=[]
 for role,row in zip(['parent','branch'],plan['source_parts']):
  mesh=reader.read_part(a,parts[row['id']],frame=FRAME);assert mesh.geometry_sha==row['geometry_sha256'] and mesh.attributes_sha==row['attributes_sha256'];mesh.save(OUT/(role+'-source-untouched.npz'));source_rows.append(dict(role=role,id=row['id'],geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha))
 specs=[];radii=[]
 for role,rkey,lo,hi in [('parent','parent_station_rings',33,41),('branch','branch_root_station_rings',176,180)]:
  oldv=z[role+'_vertices'];oldf=z[role+'_faces'];faceids=z[role+'_local_patch_face_ids'];used=np.unique(oldf[faceids]);lookup=np.full(len(oldv),-1,int);lookup[used]=np.arange(len(used));v=oldv[used].copy();f=lookup[oldf[faceids]];inner=v.copy();rings=z[rkey][lo:hi+1];profiles=[]
  for station,ring in enumerate(rings,start=lo):
   q=oldv[ring];c=q.mean(0);length=np.linalg.norm(q-c,axis=1);R=float(length.mean());t=min(.00012*np.sqrt(R/.001),.35*R);assert 0<t<R
   inner[lookup[ring]]=c+(q-c)*((R-t)/R);profiles.append(dict(station=station,radius_m=R,thickness_m=t,inner_radius_m=R-t,outer_radius_max_deviation_m=float(abs(length-R).max())))
  loops=boundary_loops(f);assert len(loops)==2
  anchors={}
  if role=='branch':
   root=set(lookup[z[rkey][180]].tolist());which=[i for i,loop in enumerate(loops)if set(loop.tolist())==root];assert len(which)==1;anchors[which[0]]=oldv[z[rkey][180]].mean(0)
  outer,fluid=close_branch(v,inner,f,anchors,origin=v.mean(0))
  for name,cell in [('outer',outer),('fluid',fluid)]:
   if signed_volume(cell['positions'],cell['faces'])<0:cell['faces']=cell['faces'][:,[0,2,1]]
   topo=topology(Mesh(cell['positions'],cell['faces'],FRAME));assert topo['closed_oriented_manifold'] and len(topo['components'])==1
   cell.update(role=role+'-'+name,source_face_ids=faceids,source_used_vertex_ids=used)
  specs.extend([outer,fluid]);radii.append(dict(role=role,profiles=profiles));atomic_npz(OUT/(role+'-derived-paired-chart.npz'),outer_positions=v,inner_positions=inner,faces=f,source_face_ids=faceids,source_used_vertex_ids=used,frame=np.array(FRAME))
 # Every point on the artificial internal root fan must lie inside parent
 # fluid. A convex parent cross-section and actual root fan certify this
 # finite classifier continuation, not a gap inferred from unsigned distance.
 parent_fluid=specs[1];branch_outer=specs[2];branch_fluid=specs[3];anchor=z['branch_vertices'][z['branch_root_station_rings'][180]].mean(0)
 proof=winding_contains(Mesh(parent_fluid['positions'],parent_fluid['faces'],FRAME),np.r_[anchor[None],branch_fluid['positions'][np.unique(branch_fluid['faces'][branch_fluid['junction_faces']])]],boundary_tolerance_m=1e-10)
 assert all(proof['inside']) and not any(proof['boundary_ambiguous'])
 atomic_json(OUT/'derived-profile-proof.json',dict(profiles=radii,root_anchor_m=anchor.tolist(),root_anchor_parent_station_36_residual_m=float(np.linalg.norm(anchor-z['parent_vertices'][z['parent_station_rings'][36]].mean(0))),root_fluid_fan_vertices_containment=proof,source_rows=source_rows))
 # Translate to millimetres near origin for native robustness. Each input,
 # result and pending-operation manifest is durable BEFORE the native call.
 origin=anchor;solids={};ids={};cells={};op=0
 def checkpoint(name,arrays=None,extra=None):
  nonlocal op
  op+=1;stem=f'{op:02d}-{name}'
  if arrays is not None:atomic_npz(OUT/(stem+'-input.npz'),**arrays)
  atomic_json(OUT/(stem+'-pending.json'),dict(operation=name,elapsed_seconds=time.monotonic()-start,input_archive=stem+'-input.npz'if arrays is not None else None,**(extra or {})))
  return stem
 for cell in specs:
  role=cell['role'];cells[role]=cell;vp=(cell['positions']-origin)*1000.;stem=checkpoint('make-'+role,dict(positions_mm=vp,faces=cell['faces'],disk_faces=cell['disk_faces'],junction_faces=cell['junction_faces']))
  solids[role],ids[role]=make_solid(vp,cell['faces']);atomic_json(OUT/(stem+'-done.json'),dict(status=str(solids[role].status()),original_id=ids[role]))
 results={}
 for kind in ['outer','fluid']:
  stem=checkpoint('union-'+kind,extra=dict(input_roles=['parent-'+kind,'branch-'+kind],input_checkpoints_complete=True));solid=md.Manifold.batch_boolean([solids['parent-'+kind],solids['branch-'+kind]],md.OpType.Add);assert solid.status()==md.Error.NoError;results[kind]=solid
  checkpoint('extract-'+kind,extra=dict(preceding_union=stem));out=solid.to_mesh64();origins,faceids=provenance(out);world=np.asarray(out.vert_properties)[:,:3]/1000.+origin;faces=np.asarray(out.tri_verts,dtype=np.int64);atomic_npz(OUT/(kind+'-union.npz'),vertices=world,faces=faces,original_id=origins,source_face_id=faceids,frame=np.array(FRAME));atomic_json(OUT/(stem+'-done.json'),dict(status=str(solid.status()),volume_m3=float(solid.volume())/1e9))
 stem=checkpoint('subtract-fluid-from-outer',extra=dict(input_archives=['outer-union.npz','fluid-union.npz']));wall=results['outer']-results['fluid'];assert wall.status()==md.Error.NoError
 checkpoint('extract-wall',extra=dict(preceding_difference=stem));out=wall.to_mesh64();origins,faceids=provenance(out);world=np.asarray(out.vert_properties)[:,:3]/1000.+origin;faces=np.asarray(out.tri_verts,dtype=np.int64);atomic_npz(OUT/'wall-raw.npz',vertices=world,faces=faces,original_id=origins,source_face_id=faceids,frame=np.array(FRAME));atomic_json(OUT/(stem+'-done.json'),dict(status=str(wall.status()),volume_m3=float(wall.volume())/1e9))
 classifier_rows=[]
 for role,cell in cells.items():classifier_rows.append(dict(role=role,original_id=ids[role],source_faces=len(cell['source_face_ids']),disk_faces=cell['disk_faces'].tolist(),junction_faces=cell['junction_faces'].tolist(),mouths=cell['mouths']))
 atomic_json(OUT/'generation-complete.json',dict(scope='Uninstalled educational finite three-port local candidate, pending independent intrinsic checks and actual render/sections.',installed=False,classifiers=classifier_rows,source_rows=source_rows,origin_m=origin.tolist(),wall_faces=len(faces),wall_volume_m3=float(wall.volume())/1e9,outer_volume_m3=float(results['outer'].volume())/1e9,fluid_volume_m3=float(results['fluid'].volume())/1e9,elapsed_seconds=time.monotonic()-start));print(json.dumps(dict(output=str(OUT),wall_faces=len(faces),elapsed_seconds=time.monotonic()-start)))

if __name__=='__main__':main()
