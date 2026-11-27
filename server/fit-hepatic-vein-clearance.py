"""Bounded native checks for a small, fixed-cava composition of two shears."""
from pathlib import Path
import copy,json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.open3d_backend import library
from geometry_tools.hepatic_course_field import HepaticCourseField,HepaticVeinShear
base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'hepatic-vein-clearance-v196';out.mkdir(exist_ok=False);started=time.monotonic();a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres'
def solid(v,f):
    u,inv=np.unique(v,axis=0,return_inverse=True);s=md.Manifold(md.Mesh64(np.ascontiguousarray(u),np.ascontiguousarray(inv[f],dtype=np.uint64)))
    if s.status()!=md.Error.NoError:raise ValueError('Finite-wall Boolean input failed')
    return s
others={}
for role,version in [('artery','v170'),('portal','v181'),('biliary','v181')]:
    r=json.loads((base/('hepatic-contained-banks-'+version)/(role+'-evidence.json')).read_text());l=r['layout'];raw=Path(r['bank']).read_bytes();assert digest_file(r['bank'])==r['bankSha256'];v=np.ndarray((l['vertexCount'],3),'<f4',buffer=raw,offset=l['positions']).astype(float)+l['geometryOrigin'];f=np.ndarray((l['indexCount']//3,3),'<u4',buffer=raw,offset=l['indices']);others[role]=solid(v,f)
source=reader.read_part(a,parts['TA98HEPATIC:male:hepatic_vein'],frame=frame);cava=reader.read_part(a,parts['TA98HEPATIC:male:region:A12.3.09.001'],frame=frame);protection=HepaticCourseField(a,base/'hepatic-contained-course-v163/field.npz',base/'separated-courses-v105/evidence.json')
old=base/'hepatic-contained-banks-v191';r=json.loads((old/'hepatic_vein-evidence.json').read_text());l=r['layout'];raw=Path(r['bank']).read_bytes();f=np.ndarray((l['indexCount']//3,3),'<u4',buffer=raw,offset=l['indices']);sv=np.load(old/'hepatic_vein-parents.npz')['source_vertices'];edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0);p=np.unique(np.r_[sv[np.unique(f)],sv[f].mean(1),sv[edges].mean(1)],axis=0)
lib=library();caval_distance=protection.cava.compute_distance(lib.core.Tensor((p-protection.query_origin).astype(np.float32)),nthreads=1).numpy();keep=(protection.weight(p)>.15)&(caval_distance>.012);p=p[keep]
body=Mesh.load(base/'vessel-containment-v157/liver-envelope.npz');origin=body.vertices.mean(0);scene=lib.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(lib.core.Tensor((body.vertices-origin).astype(np.float32)),lib.core.Tensor(body.faces.astype(np.uint32)))
record=json.loads((base/'hepatic-vein-lift-v190/field.json').read_text());record.update(startX=float(cava.vertices[:,0].max())+.012,endX=.045,cavalOutletGuardMetres=.012,
    outletReference='https://pubmed.ncbi.nlm.nih.gov/17957095/',outletPolicy='Retain unchanged native caval outlet pedicles within 12 mm of the actual caval wall; they are not intrahepatic branches.')
history=[];chosen=None
choices=sorted([(y*.001,z*.001)for y in [3.5,4.,4.5,5.]for z in [0.,-.5,-1.,-1.5,.5,1.]],key=lambda q:q[0]**2+q[1]**2)
for dy,dz in choices:
    record['shiftYZMetres']=[dy,dz];field=HepaticVeinShear(record);assert np.array_equal(field(cava.vertices),cava.vertices)
    points=field(p);d=scene.compute_signed_distance(lib.core.Tensor((points-origin).astype(np.float32)),nthreads=1,nsamples=7).numpy();s=solid(field(source.vertices),source.faces);volumes={role:float((s^other).volume())for role,other in others.items()};row=dict(shiftYZMetres=[dy,dz],maximumOutsideMetres=float(d.max(initial=0)),outsideSamples=int((d>3e-7).sum()),overlapVolumesMetresCubed=volumes);history.append(row);print(json.dumps(row),flush=True)
    if row['outsideSamples']==0 and all(v==0 for v in volumes.values()):chosen=copy.deepcopy(record);break
atomic_json(out/'search.json',dict(history=history,geometryAccepted=False,seconds=time.monotonic()-started))
if chosen is None:raise ValueError('No fixed-cava shear candidate passed the bounded clearance search')
chosen.update(search=str(out/'search.json'),cavalCoordinatesLiteral=True,determinantIdenticallyOne=True);atomic_json(out/'field.json',chosen)
