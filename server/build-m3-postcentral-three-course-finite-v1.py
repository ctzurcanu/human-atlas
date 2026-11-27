"""Complete three-native-course finite candidate; educational inner profiles.

Four native terminal annuli replace the local prototype's artificial rims.
Every native-operation input is checkpointed; failures are terminal, no retry.
"""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology,winding_contains
from geometry_tools.algorithms.branch_junction_envelopes import close_branch
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops,signed_volume
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance
import manifold3d as md
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';OUT=P/'three-course-finite-educational-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);z=np.load(F/'three-course-source-charts.npz');feas=json.loads((F/'feasibility.json').read_text());assert feas['finite_construction_feasible'];origin=z['M3_vertices'][z['M3_rings'][36]].mean(0);sources={'M3':Mesh.load(P/'finite-junction-educational-v1/parent-source-untouched.npz'),'main':Mesh.load(P/'finite-junction-educational-v1/branch-source-untouched.npz')};source_rows={name:dict(geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha)for name,m in sources.items()}
 atomic_json(OUT/'construction-contract.json',dict(source_charts_sha256=digest_file(F/'three-course-source-charts.npz'),source_feasibility_sha256=digest_file(F/'feasibility.json'),frame=FRAME,installed=False,three_courses=True,native_terminal_mouths=4,thickness_formula='t(R)=min(0.12 mm * sqrt(R / 1 mm),0.35 R)',thickness_status='Tunable educational construction assumption, not clinical measurement, source-recovered finite wall or exact named histological coat.',trend_source='https://pubmed.ncbi.nlm.nih.gov/10503931/',trend_limit='Relative media thickness rises as human pial arteries narrow; this article does not supply the construction formula.',fixed_outer_domains=dict(M3_except=[33,41],main_except=[[115,123],[176,180]],twig_except=[19,24]),no_artificial_internal_patch_rims=True,fixed_brain_skeleton_and_other_vessels=True,source_rows=source_rows))
 specs=[];profile_rows=[]
 for role,vkey,fkey,rkey,faceidskey in [('M3','M3_vertices','M3_faces','M3_rings',None),('main','branch_vertices','branch_faces','main_rings','main_source_face_ids'),('twig','branch_vertices','branch_faces','twig_rings','twig_source_face_ids')]:
  oldv=z[vkey];oldf=z[fkey];faceids=np.arange(len(oldf))if faceidskey is None else z[faceidskey];used=np.unique(oldf[faceids]);lookup=np.full(len(oldv),-1,int);lookup[used]=np.arange(len(used));outer=oldv[used].copy();inner=outer.copy();faces=lookup[oldf[faceids]];profiles=[]
  for station,ring in enumerate(z[rkey]):
   q=oldv[ring];center=q.mean(0);R=float(np.linalg.norm(q-center,axis=1).mean());t=min(.00012*np.sqrt(R/.001),.35*R);assert 0<t<R;inner[lookup[ring]]=center+(q-center)*((R-t)/R);profiles.append(dict(station=station,outer_radius_m=R,thickness_m=t,inner_radius_m=R-t))
  anchors={};loops=boundary_loops(faces);assert len(loops)==2
  if role!='M3':
   ring=z[rkey][-1];ids=set(lookup[ring].tolist());which=[i for i,loop in enumerate(loops)if set(loop.tolist())==ids];assert len(which)==1;anchors[which[0]]=oldv[ring].mean(0)
  outercell,fluidcell=close_branch(outer,inner,faces,anchors,origin=origin)
  for kind,cell in [('outer',outercell),('fluid',fluidcell)]:
   if signed_volume(cell['positions'],cell['faces'])<0:cell['faces']=cell['faces'][:,[0,2,1]]
   t=topology(Mesh(cell['positions'],cell['faces'],FRAME));assert t['closed_oriented_manifold']and len(t['components'])==1
   cell.update(role=role+'-'+kind,course=role,source_face_ids=faceids,source_used_vertex_ids=used);specs.append(cell)
  atomic_npz(OUT/(role+'-derived-paired-chart.npz'),outer_positions=outer,inner_positions=inner,faces=faces,source_face_ids=faceids,source_used_vertex_ids=used,rings=lookup[z[rkey]],frame=np.array(FRAME));profile_rows.append(dict(course=role,profiles=profiles))
 # The actual internal closing fans must be strictly in their parent lumen;
 # they are construction classifiers, never anatomical end membranes.
 containment=[]
 for parent,child in [('M3','main'),('main','twig')]:
  pc=next(s for s in specs if s['role']==parent+'-fluid');cc=next(s for s in specs if s['role']==child+'-fluid');points=cc['positions'][np.unique(cc['faces'][cc['junction_faces']])];proof=winding_contains(Mesh(pc['positions'],pc['faces'],FRAME),points,boundary_tolerance_m=1e-10);assert all(proof['inside'])and not any(proof['boundary_ambiguous']);containment.append(dict(parent=parent,child=child,proof=proof))
 atomic_json(OUT/'derived-profile-proof.json',dict(profiles=profile_rows,actual_root_fans_containment=containment))
 op=0;solids={};ids={};native_inputs=[]
 def checkpoint(name,arrays=None,extra=None):
  nonlocal op
  op+=1;stem=f'{op:02d}-{name}'
  if arrays is not None:atomic_npz(OUT/(stem+'-input.npz'),**arrays)
  atomic_json(OUT/(stem+'-pending.json'),dict(operation=name,elapsed_seconds=time.monotonic()-start,input_archive=stem+'-input.npz'if arrays is not None else None,**(extra or {})));return stem
 for cell in specs:
  role=cell['role'];vp=(cell['positions']-origin)*1000.;stem=checkpoint('make-'+role,dict(positions_mm=vp,faces=cell['faces'],source_face_ids=cell['source_face_ids'],disk_faces=cell['disk_faces'],junction_faces=cell['junction_faces']));solids[role],ids[role]=make_solid(vp,cell['faces']);native_inputs.append(dict(role=role,checkpoint_input=stem+'-input.npz'));atomic_json(OUT/(stem+'-done.json'),dict(status=str(solids[role].status()),original_id=ids[role]))
 results={}
 for kind in ['outer','fluid']:
  roles=[r+'-'+kind for r in ['M3','main','twig']];stem=checkpoint('union-'+kind,extra=dict(input_roles=roles));solid=md.Manifold.batch_boolean([solids[r]for r in roles],md.OpType.Add);assert solid.status()==md.Error.NoError;results[kind]=solid;checkpoint('extract-'+kind,extra=dict(preceding_union=stem));mesh=solid.to_mesh64();origins,faceids=provenance(mesh);atomic_npz(OUT/(kind+'-union.npz'),vertices=np.asarray(mesh.vert_properties)[:,:3]/1000.+origin,faces=np.asarray(mesh.tri_verts,dtype=np.int64),original_id=origins,source_face_id=faceids,frame=np.array(FRAME));atomic_json(OUT/(stem+'-done.json'),dict(status=str(solid.status()),volume_m3=float(solid.volume())/1e9))
 stem=checkpoint('subtract-fluid-from-outer',extra=dict(input_archives=['outer-union.npz','fluid-union.npz']));wall=results['outer']-results['fluid'];assert wall.status()==md.Error.NoError;checkpoint('extract-wall',extra=dict(preceding_difference=stem));mesh=wall.to_mesh64();origins,faceids=provenance(mesh);atomic_npz(OUT/'wall-raw.npz',vertices=np.asarray(mesh.vert_properties)[:,:3]/1000.+origin,faces=np.asarray(mesh.tri_verts,dtype=np.int64),original_id=origins,source_face_id=faceids,frame=np.array(FRAME));atomic_json(OUT/(stem+'-done.json'),dict(status=str(wall.status()),volume_m3=float(wall.volume())/1e9))
 rows=[dict(role=s['role'],course=s['course'],original_id=ids[s['role']],source_faces=len(s['source_face_ids']),source_face_ids=s['source_face_ids'].tolist(),disk_faces=s['disk_faces'].tolist(),junction_faces=s['junction_faces'].tolist(),mouths=s['mouths'])for s in specs]
 atomic_json(OUT/'generation-complete.json',dict(installed=False,scope='Complete native three-course finite educational candidate, pending source replay, intrinsic checks and actual renders/sections.',source_rows=source_rows,origin_m=origin.tolist(),native_inputs=native_inputs,classifiers=rows,wall_faces=len(mesh.tri_verts),wall_volume_m3=float(wall.volume())/1e9,outer_volume_m3=float(results['outer'].volume())/1e9,fluid_volume_m3=float(results['fluid'].volume())/1e9,elapsed_seconds=time.monotonic()-start));print(json.dumps(dict(wall_faces=len(mesh.tri_verts),elapsed_seconds=time.monotonic()-start,output=str(OUT))))

if __name__=='__main__':main()
