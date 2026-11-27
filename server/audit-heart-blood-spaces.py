"""Inspect actual installed myocardial boundaries and computation-only root gates.

The gates never become tissue, viewer geometry or selectable anatomy. Their
purpose is to distinguish the external space from the two blood circuits.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import manifold3d as md
import igl
from scipy.spatial import ConvexHull
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.algorithms.annular_pipe_envelopes import face_components,signed_volume
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.oriented_shell_regions import occupied_shell_regions

ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--gate-stations');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());m=a['finiteCardiacIntegration'];pp={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';o=np.array([0.,1.375,.03125])
parts=[reader.read_part(a,pp[pid],frame=frame,catalogue_path=path)for pid in m['wallPartIds']]
positions=np.concatenate([p.vertices for p in parts]);offsets=np.cumsum([0]+[len(p.vertices)for p in parts[:-1]]);faces=np.concatenate([p.faces+offset for p,offset in zip(parts,offsets)])
v,inv=np.unique(positions,axis=0,return_inverse=True);faces=inv[faces];keys=np.sort(faces,axis=1);_,first,count=np.unique(keys,axis=0,return_index=True,return_counts=True);assert count.max()==2
boundary=faces[first[count==1]];wall=Mesh(v,boundary,frame);wall.save(out/'installed-boundary.npz');tissue=md.Manifold(md.Mesh64(v-o,boundary.astype('uint64')));assert str(tissue.status())=='Error.NoError'
nc,labels=face_components(boundary);components=[]
for k in range(nc):
 f=boundary[labels==k];ids=np.unique(f);edges=np.unique(np.sort(np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]],axis=1),axis=0);euler=len(ids)-len(edges)+len(f);components.append(dict(faces=len(f),genus=(2-euler)//2,volume_mL=signed_volume(v,f)*1e6))
gates=[];gate_records=[];atrium=np.array([.021,1.40,-.004]);whole=np.array([.02,1.38,.03]);declared=json.loads(Path(args.gate_stations).read_text())if args.gate_stations else{}
roots=['ZA:Ascending aorta','ZA:Pulmonary trunk','ZA:Superior vena cava','TA98HEPATIC:male:region:A12.3.09.001','ZA:Left inferior pulmonary vein','ZA:Left superior pulmonary vein','ZA:Right inferior pulmonary vein','ZA:Right superior pulmonary vein']
for pid in roots:
 # The hepatic replacement is a closed branched region. Its retained native
 # port establishes the same documented historical ring, not a new vessel.
 ring_source='ZA:Inferior vena cava (thoracic part)'if pid.startswith('TA98HEPATIC:')else pid
 root=reader.read_part(a,pp[ring_source],frame=frame,catalogue_path=path);rv,ri=np.unique(root.vertices,axis=0,return_inverse=True);rings=boundary_loops(ri[root.faces]);reference=atrium if'pulmonary vein'in pid else whole
 loop=min(rings,key=lambda q:np.linalg.norm(rv[q].mean(0)-reference));points=rv[loop];c=points.mean(0);_,_,vh=np.linalg.svd(points-c,full_matrices=False);n=vh[-1];uv=(points-c)@vh[:2].T;half=.0008;barrier_source='native mouth plane'
 if pid in declared:
  row=declared[pid];points=np.array(row['points_m']);c=np.array(row['point_m']);n=np.array(row['normal']);n/=np.linalg.norm(n);axis=np.eye(3)[np.argmin(abs(n))];u=np.cross(n,axis);u/=np.linalg.norm(u);vv=np.cross(n,u);vh=np.array([u,vv,n]);uv=(points-c)@vh[:2].T;uv=uv[ConvexHull(uv).vertices];half=.00005;barrier_source='convex envelope of measured named-root lumen contour; 2% diagnostic radial overlap'
 area=np.sum(uv[:,0]*np.roll(uv[:,1],-1)-uv[:,1]*np.roll(uv[:,0],-1))
 if area<0:uv=uv[::-1]
 cap=np.asarray(md.triangulate([uv],epsilon=0.,allow_convex=False),int);assert len(cap)==len(uv)-2,pid
 # A diagnostic gate extends 1.5% into the source boundary and is 1.6 mm
 # thick. This is a numerical separator, never an anatomical root radius.
 uv*=1.02 if pid in declared else 1.015;poly=c+uv@vh[:2];N=len(poly);verts=np.r_[poly-half*n,poly+half*n];ff=[]
 for j in range(N):
  k=(j+1)%N;ff.extend([[j,k,N+k],[j,N+k,N+j]])
 ff=np.r_[np.array(ff),cap[:,[0,2,1]],cap+N]
 if signed_volume(verts,ff)<0:ff=ff[:,[0,2,1]]
 gate=md.Manifold(md.Mesh64(verts-o,ff.astype('uint64')));assert str(gate.status())=='Error.NoError'
 overlap=float((gate^tissue).volume()*1e6);gates.append(gate);gate_records.append(dict(part_id=pid,ring_source=ring_source,center_m=c.tolist(),normal=n.tolist(),source_controls=points.tolist(),barrier_source=barrier_source,half_thickness_m=half,maximum_source_nonplanarity_m=float(abs((points-c)@n).max()),gate_tissue_overlap_mL=overlap,diagnostic_only=True,rendered=False));atomic_json(out/f'gate-{len(gates)}.json',gate_records[-1])
low=v.min(0)-.025;high=v.max(0)+.025;box=md.Manifold.cube(high-low).translate(low-o);free=box-tissue-md.Manifold.batch_boolean(gates,md.OpType.Add);assert str(free.status())=='Error.NoError';spaces=free.decompose();rows=[];shells=[];seeds={'right-atrium':[-.022,1.39,.022],'right-ventricle':[.025,1.36,.06],'left-atrium':[.018,1.405,-.012],'left-ventricle':[.05,1.36,.025]};q=np.array(list(seeds.values()))
for i,space in enumerate(spaces):
 z=space.to_mesh64();sv=z.vert_properties[:,:3]+o;sf=z.tri_verts;row=dict(boundary_shell=i,oriented_volume_mL=float(space.volume()*1e6),triangles=len(sf));rows.append(row);mesh=Mesh(sv,sf,frame);shells.append(mesh);mesh.save(out/f'fluid-boundary-{i}.npz')
regions=occupied_shell_regions(shells,q)
for i,row in enumerate(regions):row.update(component=i,seeds=[name for k,name in enumerate(seeds)if k in row['query_indices']],volume_mL=row['volume_m3']*1e6)
assign={name:[r['component']for r in regions if name in r['seeds']]for name in seeds};valid=all(len(ids)==1 for ids in assign.values());separated=valid and assign['right-atrium']==assign['right-ventricle']and assign['left-atrium']==assign['left-ventricle']and assign['left-atrium']!=assign['right-atrium'];result=dict(computed=True,passed=bool(separated and len(regions)==3 and all(r['gate_tissue_overlap_mL']>0 for r in gate_records)),catalogue_sha256=digest_file(path),wall_geometry_sha256=wall.geometry_sha,wall_components=components,space_components=regions,boundary_shells=rows,seed_components=assign,circuits_separated_by_this_gate_configuration=bool(separated),gates=gate_records,gate_dimensions_provisional=True,actual_valve_obstruction_checked=False,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False)
atomic_json(out/'blood-spaces.json',result);print(json.dumps({k:v for k,v in result.items()if k!='gates'}),flush=True)
raise SystemExit(0 if result['passed']else 2)
