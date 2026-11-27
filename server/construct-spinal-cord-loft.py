"""Conforming continuous grey/white core, clipped by the native cord envelope.

Continuous lofts avoid the disconnected pockets produced by undersampling a
very small canal in a volumetric grid. Anatomical profiles vary with level.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import manifold3d as md
from scipy.interpolate import CubicSpline
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.boolean_cells import BooleanCells

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic()
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True);assert not any(out.iterdir())
z=np.load(base/'female-surface-v1/surface.npz');v,inv=np.unique(z['positions'],axis=0,return_inverse=True);f=inv[z['faces']]
domain=np.load(base/'female-domain-v1/domain.npz');sample_y=domain['sample_y'];dims=domain['sample_dimensions']
right=np.array([[0.,.115],[.18,.15],[.25,.35],[.37,.53],[.52,.54],[.64,.41],
    [.61,.24],[.43,.065],[.41,-.13],[.43,-.33],[.50,-.54],[.57,-.77],
    [.40,-.63],[.27,-.39],[.20,-.17],[0.,-.115]])
controls=np.r_[right,right[-2:0:-1]*[-1.,1.],right[:1]]
curve=CubicSpline(np.arange(len(controls)),controls,bc_type='periodic')
profile=curve(np.linspace(0,len(controls)-1,64,endpoint=False))
ys=np.r_[v[:,1].min()-.003,np.linspace(sample_y[0],sample_y[-1],72),v[:,1].max()+.003]
grey=[];canal=[]
for y in ys:
    cx,cz,rx,rz=[np.interp(y,sample_y,dims[:,i])for i in range(4)]
    scale=np.interp(y,[1.148,1.158,1.178,1.198,1.24,1.42,1.47,1.53,1.565],[1.16,1.16,1.08,.95,.69,.69,.86,.86,.81])
    p=profile.copy();p[:,0]*=scale;p[:,1]=np.where(p[:,1]>0,p[:,1]*scale,p[:,1])
    grey.append(np.c_[cx+p[:,0]*rx,np.full(len(p),y),cz+p[:,1]*rz])
    angle=np.linspace(0,2*np.pi,12,endpoint=False);canal.append(np.c_[cx+.00008*np.cos(angle),np.full(len(angle),y),cz+.00008*np.sin(angle)])

def loft(rings):
    p=np.concatenate(rings);n=len(rings[0]);faces=[]
    for k in range(len(rings)-1):
        for i in range(n):
            j=(i+1)%n;a=k*n+i;b=k*n+j;c=(k+1)*n+j;d=(k+1)*n+i
            faces.extend([[a,b,c],[a,c,d]])
    # Concave grey horns use constrained polygon triangulation, never a
    # centre fan spanning the notches between horns.
    for k,side in [(0,-1),(len(rings)-1,1)]:
        xy=rings[k][:,[0,2]];poly_area=np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1))
        # Manifold's triangulator requires a CCW exterior contour. A clockwise
        # exterior must be reversed before triangulating, not triangle by triangle.
        order=np.arange(n) if poly_area>0 else np.arange(n-1,-1,-1)
        cap=order[np.asarray(md.triangulate([xy[order]],epsilon=0.,allow_convex=False),int)]
        t=xy[cap];cap_area=(t[:,1,0]-t[:,0,0])*(t[:,2,1]-t[:,0,1])-(t[:,1,1]-t[:,0,1])*(t[:,2,0]-t[:,0,0])
        reverse=cap_area*poly_area*side<0;cap[reverse]=cap[reverse][:,[0,2,1]]
        faces.extend(cap+k*n)
    ff=np.asarray(faces,int);t=topology(Mesh(p,ff,'spinal-profile-construction'))
    if t['oriented_triangle_integral_m3']<0:ff=ff[:,[0,2,1]];t=topology(Mesh(p,ff,'spinal-profile-construction'))
    if not(t['closed_oriented_manifold'] and len(t['components'])==1 and t['components'][0]['closed_orientable_genus']==0):
        atomic_npz(out/'rejected-loft.npz',positions=p,faces=ff)
        atomic_json(out/'rejected-loft.json',t)
        raise ValueError('Loft boundary gate: '+str({k:t[k]for k in ['open_edges','nonmanifold_edges','inconsistent_winding_edges','zero_area_faces','duplicate_faces']}))
    return p,ff

gv,gf=loft(grey);cv,cf=loft(canal)
b=BooleanCells(tolerance_m=2e-12);whole=b.source(v,f,'native-cord-envelope');g=b.source(gv,gf,'continuous-variable-grey-profile');c=b.source(cv,cf,'central-canal-void')
solids={'grey':(whole^g)-c,'white':whole-g}
for name,s in solids.items():assert str(s.status())=='Error.NoError' and not s.is_empty(),name
b.save(solids,out/'boolean-source.npz');result=b.assemble(solids)
arrays={'positions':result['positions'], 'grey_classifier_positions':gv,'grey_classifier_faces':gf,'canal_classifier_positions':cv,'canal_classifier_faces':cf}
for k,fs in result['surfaces'].items():arrays[k]=fs
atomic_npz(out/'cells.npz',**arrays)
atomic_json(out/'cells.json',dict(cells=result['cells'],surfaces={k:len(fs)for k,fs in result['surfaces'].items()},
    construction=result['evidence'],source_surface_sha256=digest_file(base/'female-surface-v1/surface.npz'),
    boundary_classification='Continuous reference-guided regional grey/white profiles; microscopic laminae and nuclei remain unbuilt.',
    canal_representation='Continuous 0.08 mm radius lumen, removed from tissue; not a solid cylinder',
    positions=len(result['positions']),anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(cells=list(result['cells']),surfaces={k:len(fs)for k,fs in result['surfaces'].items()},
    volumes_mL={k:float(s.volume())*1e6 for k,s in solids.items()},elapsed_seconds=time.monotonic()-start)),flush=True)
