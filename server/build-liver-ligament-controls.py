"""Small, source-bound fibrous cords following the two actual fissure charts."""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader
from geometry_tools.open3d_backend import library
from geometry_tools.algorithms.tube_sections import centripetal_curve,rotation_minimizing_frames
from geometry_tools.algorithms.tube_continuity import loft_faces
from geometry_tools.tissue import pack

start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'ligament-courses-v153';out.mkdir(exist_ok=False)
a=json.loads((base/'liver-surface-controls-v143/with-landmark.json').read_text());host=next(p for p in a['parts']if p['id']=='TA98PART:ZA:Liver:capsule');m=AtlasReader().read_part(a,host,frame='model:local-ta98:assembled-metres')
outer=m.faces[:host['tissueDrawGroups'][0]['count']//3];o3d=library();scene=o3d.t.geometry.RaycastingScene(nthreads=1)
scene.add_triangles(o3d.core.Tensor(m.vertices.astype(np.float32)),o3d.core.Tensor(outer.astype(np.uint32)))
points=[];surfaces={};cells={};evidence={}
for name,number in [('venosum',10),('teres',14)]:
    g=next(g for g in a['liverFeatureIntegration']['termGroups']if g['term']==f'A05.8.01.{number:03d}');centers=m.triangles[g['surface']['faces']].mean(1);center=centers.mean(0)
    axis=np.linalg.svd(centers-center,full_matrices=False)[2][0]
    if axis[1]<0:axis=-axis
    s=(centers-center)@axis;bins=np.linspace(s.min(),s.max(),7);controls=[]
    for i in range(6):
        chosen=centers[(s>=bins[i])&(s<=bins[i+1])]
        if len(chosen):controls.append(chosen.mean(0))
    path,_=centripetal_curve(np.array(controls),max_step_m=.002,chord_tolerance_m=.00005)
    hit=scene.compute_closest_points(o3d.core.Tensor(path.astype(np.float32)),nthreads=1)
    skin=hit['points'].numpy().astype(float);normal=hit['primitive_normals'].numpy().astype(float)
    radius=.00045 if name=='venosum'else .00065
    # A short projection offset keeps the illustrative cord outside the host;
    # no capsule holes or new hilar solid disks are constructed.
    path=skin+normal*(radius*1.1);frames=rotation_minimizing_frames(path);angles=np.arange(12)*2*np.pi/12
    taper=.55+.45*np.sin(np.linspace(0,np.pi,len(path)))
    rings=path[:,None,:]+radius*taper[:,None,None]*(frames[:,:,0][:,None,:]*np.cos(angles)[None,:,None]+frames[:,:,1][:,None,:]*np.sin(angles)[None,:,None])
    vv=rings.reshape(-1,3);ff=loft_faces(len(path),12)
    # Solid fibrous cords have genuine terminal faces, unlike blood lumina.
    vv=np.r_[vv,path[[0,-1]]];first=len(vv)-2;last=len(vv)-1;end=(len(path)-1)*12
    caps=np.array([[first,(i+1)%12,i]for i in range(12)]+[[last,end+i,end+(i+1)%12]for i in range(12)],int)
    f=np.r_[ff,caps];offset=len(points);points.extend(vv.tolist());sid=name;surfaces[sid]=f+offset;cells[name]=[dict(surface=sid,side=1,neighbor=None)]
    evidence[name]=dict(hostGeometrySha256=m.geometry_sha,surfaceTerm=g['term'],radiusParameterMetres=radius,dimensionStatus='Illustrative cord calibre; anatomical fissure course follows actual source surface',path=path.tolist(),attachmentExtentProvisional=True)
ids={name:'TA98LIVER:male:ligament:'+name for name in cells}
result=pack(np.array(points),surfaces,cells,ids,m.frame,out/'packed',origin=[0,1.25,0],max_pairs=2000000,max_tests=16000000,checkpoint_directory=out/'checkpoints')
atomic_json(out/'evidence.json',dict(passed=result['passed'],courses=evidence,elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(passed=result['passed'],failures=result['checks']['failures'],bytes=result['total_bytes'],elapsedSeconds=time.monotonic()-start)),flush=True)
