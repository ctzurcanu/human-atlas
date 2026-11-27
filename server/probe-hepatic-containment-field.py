"""A shared, root-preserving course field; measure before accepting geometry."""
from pathlib import Path
import json,sys,time,argparse
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.open3d_backend import library

ap=argparse.ArgumentParser();ap.add_argument('--version',default='v158');ap.add_argument('--strength',type=float,default=1.);ap.add_argument('--caval-protect',type=float,default=.003);args=ap.parse_args()
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/('hepatic-contained-course-'+args.version);out.mkdir(exist_ok=False)
a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());r=AtlasReader();parts={p['id']:p for p in a['parts']};o3d=library()
body=Mesh.load(base/'vessel-containment-v157/liver-envelope.npz');scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor(body.vertices.astype(np.float32)),o3d.core.Tensor(body.faces.astype(np.uint32)))
cava=r.read_part(a,parts['TA98HEPATIC:male:region:A12.3.09.001'],frame=body.frame);cv=o3d.t.geometry.RaycastingScene(nthreads=1);cv.add_triangles(o3d.core.Tensor(cava.vertices.astype(np.float32)),o3d.core.Tensor(cava.faces.astype(np.uint32)))
z=np.load(base/'superior-fit-v8/field.npz');floor=float(z['floor_y']);roof=CloughTocher2DInterpolator(z['xz'],z['source_roof_y'],fill_value=floor+.005);correction=CloughTocher2DInterpolator(z['xz'],z['correction_y'],fill_value=0.)
anchors=json.loads((base/'separated-courses-v105/evidence.json').read_text())['anchors']
def smooth(t):t=np.clip(t,0,1);return t*t*(3-2*t)
def mapped(v):
    p=np.array(v,float,copy=True);xz=p[:,[0,2]];span=np.maximum(roof(xz)-floor,.005);d=np.clip(correction(xz),0,.55*span);depth=np.clip((p[:,1]-floor)/span,0,1)
    delta=np.zeros_like(p);delta[:,1]=-args.strength*d*smooth(depth)
    # The anterior medial protrusion followed a local course offset. Apply
    # one smooth 2.2 mm posterior correction to every channel in that region.
    delta[:,2]=-.0022*np.exp(-.5*((p[:,0]+.025)/.016)**2-.5*((p[:,1]-1.264)/.018)**2)*smooth((p[:,2]-.075)/.025)
    distance=cv.compute_distance(o3d.core.Tensor(p.astype(np.float32)),nthreads=1).numpy().astype(float)
    weight=smooth((distance-args.caval_protect)/.012)
    for anchor in anchors.values():
        center=np.asarray(anchor['center_m']);radius=anchor['radius_m'];weight*=smooth((np.linalg.norm(p-center,axis=1)-radius)/.012)
    return p+weight[:,None]*delta
old=json.loads((base/'vessel-containment-v157/evidence.json').read_text());rows=[];mapped_whole={}
for role in ['artery','portal','hepatic_vein','biliary']:
    p=parts['TA98HEPATIC:male:'+role];m=r.read_part(a,p,frame=body.frame);v=mapped(m.vertices);mapped_whole[role]=(m,v);Mesh(v,m.faces,m.frame,m.attributes).save(out/(role+'.npz'))
    print(role,'mapped',float(np.linalg.norm(v-m.vertices,axis=1).max()),flush=True)
for row in old['measuredIntrahepaticBranches']:
    p=parts[row['id']];m=r.read_part(a,p,frame=body.frame);v=mapped(m.vertices);used=np.unique(m.faces);signed=scene.compute_signed_distance(o3d.core.Tensor(v[used].astype(np.float32)),nthreads=1,nsamples=7).numpy();cd=cv.compute_distance(o3d.core.Tensor(m.vertices[used].astype(np.float32)),nthreads=1).numpy()
    item=dict(term=row['term'],name=row['name'],outside_vertices=int(sum(signed>3e-7)),maximum_outside_m=float(max(0,signed.max())),
        outside_away_from_caval_join=int(sum((signed>3e-7)&(cd>.003))),maximum_outside_away_from_caval_join_m=float(max(0,np.max(signed[cd>.003],initial=0))),
        original_maximum_outside_m=row['maximum_outside_m'])
    rows.append(item);print(json.dumps(item),flush=True)
atomic_json(out/'probe.json',dict(source='Current installed hepatic networks; whole liver and all adjoining structures fixed',
    field='Source liver superior-envelope field plus local medial posterior correction; common to all four channels; protected hilar anchors and caval trunk',
    branches=rows,parameters=vars(args),geometryAccepted=False,installed=False,elapsedSeconds=time.monotonic()-start))
