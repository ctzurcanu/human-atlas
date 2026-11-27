"""Shared regional fields over pre-existing finite anatomical wall domains."""
from pathlib import Path
import argparse,sys,time,json
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--calculation',default='tissue-domains-v76');ap.add_argument('--stage',default='regional-domains-v77');ap.add_argument('--regularize-regional-ties',action='store_true');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);calc=b/args.calculation;assert json.loads((calc/'construction.json').read_text())['passed'];z=np.load(calc/'tetrahedra.npz');p=z['positions'];m=Mesh.load(b/'native-pv-mouths-v69/whole-wall.npz');o=z['origin'];charts=np.load(b/'lining-charts-v71/charts.npz');inner=charts['endocardial_faces'];owner=charts['face_owner'];normal=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None];dist=[];closest=[];normals=[]
for chamber in range(4):
 mask=np.flatnonzero(inner&(owner==chamber));scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((m.vertices-o).astype('float32')),o3d.core.Tensor(m.faces[mask].astype('uint32')));hit=scene.compute_closest_points(o3d.core.Tensor((p-o).astype('float32')),nthreads=1);q=hit['points'].numpy().astype(float)+o;dist.append(np.linalg.norm(p-q,axis=1));closest.append(q);normals.append(normal[mask[hit['primitive_ids'].numpy()]])
de=np.column_stack(dist);baseline=-de.min(1)
def septum(i,j,limit):
 a=closest[i];bb=closest[j];delta=bb-a;length=np.linalg.norm(delta,axis=1);unit=delta/np.maximum(length,1e-30)[:,None];projection=np.einsum('ij,ij->i',p-a,unit);opposition=np.einsum('ij,ij->i',normals[i],normals[j]);g=np.minimum.reduce([.5*(limit-de[:,i]-de[:,j]),projection,length-projection,.003*(-.3-opposition)]);g[length<.001]=np.minimum(g[length<.001],-.01);return g
scores=np.c_[-de,baseline+septum(0,2,.012),baseline+septum(1,3,.025),baseline+septum(0,3,.010)]
if args.regularize_regional_ties:scores-=.000001*np.sqrt(np.array([2.,3.,5.,7.,11.,13.,17.]))
else:scores=np.rint(scores/.000001)*.000001
core_ids=np.unique(z['tetrahedra'][z['domains']==0]);candidate=core_ids[np.argmax(scores[core_ids],axis=1)==5];assert len(candidate);path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());vm=AtlasReader().read_part(a,next(x for x in a['parts']if x['id']=='ZA:Non-coronary leaflet'),frame=m.frame,catalogue_path=path);apex=m.vertices[np.argmin(m.vertices[:,1]-.3*m.vertices[:,0])];direction=apex-vm.vertices.mean(0);guide=vm.vertices.mean(0)+direction/np.linalg.norm(direction)*.008;center=p[candidate[np.argmin(np.linalg.norm(p[candidate]-guide,axis=1))]];fields=[];zero=[];comparisons={}
for i in range(7):
 for j in range(i+1,7):comparisons[i,j]=len(fields);fields.append(scores[:,i]-scores[:,j]);zero.append(1)
membrane=len(fields);fields.append(.003-np.linalg.norm(p-center,axis=1));zero.append(-1)
regions=['right-atrium','right-ventricle','left-atrium','left-ventricle','interatrial-septum','interventricular-muscular-septum','atrioventricular-septum','interventricular-membranous-septum'];constraints=[];names=[];cell_domains=[]
for domain,kind in enumerate(['myocardium','endocardium','epicardium']):
 for region,key in enumerate(regions):
  i=5 if region==7 else region;refs=[(comparisons[min(i,j),max(i,j)],1 if i<j else-1)for j in range(7)if j!=i]
  if region in [5,7]:refs.append((membrane,1 if region==7 else-1))
  constraints.append(refs);names.append(('fibrous'if domain==0 and region==7 else kind)+':'+key);cell_domains.append(domain)
atomic_npz(out/'fields.npz',positions=p,tetrahedra=z['tetrahedra'],domains=z['domains'],fields=np.column_stack(fields),zero_side_priority=np.array(zero),cell_domains=np.array(cell_domains),origin=o);atomic_json(out/'cells.json',dict(names=names,constraints=constraints,regions=regions));report=dict(prepared=True,source_calculation_sha256=digest_file(calc/'tetrahedra.npz'),source_wall_geometry_sha256=m.geometry_sha,regions=regions,names=names,tetrahedra=len(z['tetrahedra']),planes=len(fields),preexisting_finite_tissue_domains=True,coatings_not_extrapolated_from_nodal_distance=True,membranous_center_m=center.tolist(),membranous_guide_m=guide.tolist(),membranous_guide_distance_m=float(np.linalg.norm(center-guide)),regional_boundaries_provisional=True,regional_tie_regularization=args.regularize_regional_ties,maximum_regional_score_bias_difference_m=float(.000001*(np.sqrt(17)-np.sqrt(2)))if args.regularize_regional_ties else 0.,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
