"""An independently generated human-liver microscopic teaching unit.

Primary model: human Histology Guide MH-126b. Geometry depicts one-cell-thick
radiating/anastomosing plates, intervening sinusoidal spaces, peripheral
portal tracts and a central vein. Dimensions are illustrative construction
parameters; this is not patient-specific histology or a copied surface.
"""
from pathlib import Path
from collections import defaultdict
import sys,json,time,math
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.tissue import pack

start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'microscopic-unit-v151';out.mkdir(exist_ok=False)
points=[];surfaces={};cells=defaultdict(list);names=['lobule','arteries','portal-venules','central-vein','bile-ductules','perivascular-sheath']
for name in names:cells[name]=[]
def surface(name,vertices,faces,owners):
    offset=len(points);points.extend(np.asarray(vertices,float).tolist());sid=name
    surfaces[sid]=np.asarray(faces,int)+offset
    for owner,side,neighbor in owners:cells[owner].append(dict(surface=sid,side=side,neighbor=neighbor))
def tube(label,center,ri,ro,height,sheath=None,segments=20):
    angle=np.arange(segments)*2*np.pi/segments;center=np.array(center,float)
    fields=[]
    for r in [ri,ro]+([sheath]if sheath else[]):
        fields.append(np.array([center+[r*np.cos(a),y,r*np.sin(a)]for y in [-height/2,height/2]for a in angle]))
    outer=[]
    for i in range(segments):j=(i+1)%segments;outer.extend([[i,j,segments+j],[i,segments+j,segments+i]])
    # X-Z ring ordering points inward in this Y-up convention; reverse.
    outer=np.array(outer,int)[:,[0,2,1]]
    for k,radius in enumerate([ri,ro]+([sheath]if sheath else[])):
        if k==0:owners=[(label,-1,None)]
        elif k==1:owners=[(label,1,'perivascular-sheath'if sheath else None)]+([('perivascular-sheath',-1,label)]if sheath else[])
        else:owners=[('perivascular-sheath',1,None)]
        surface(f'{label}:{len(surfaces)}:radial',fields[k],outer,owners)
    for k,owner in [(0,label)]+([(1,'perivascular-sheath')]if sheath else[]):
        vv=np.r_[fields[k],fields[k+1]];n=segments*2;rims=[]
        for side in [0,1]:
            for i in range(segments):
                j=(i+1)%segments;a=side*segments+i;b=side*segments+j;c=b+n;d=a+n
                rims.extend([[a,c,b],[a,d,c]]if side==0 else[[a,b,c],[a,c,d]])
        surface(f'{owner}:{len(surfaces)}:open-lumen-annuli',vv,np.array(rims,int),[(owner,1,None)])

# Millimetre construction coordinates are converted once to metres.
tube('central-vein',[0,0,0],.025,.030,.60)
for i in range(6):
    a=i*np.pi/3;radial=np.array([np.cos(a),0,np.sin(a)]);tangent=np.array([-np.sin(a),0,np.cos(a)])
    center=radial*.58
    tube('portal-venules',center,.034,.042,.60,sheath=.046)
    tube('arteries',center+tangent*.082,.010,.018,.60,sheath=.022)
    tube('bile-ductules',center-tangent*.082,.009,.017,.60,sheath=.021)

def bar(a,b,width,height):
    a,b=np.array(a,float),np.array(b,float);delta=b-a;length=np.linalg.norm(delta)
    # Cuboid length along X; its Y dimension is the plate height.
    angle=-math.degrees(math.atan2(delta[2],delta[0]));center=(a+b)/2
    return md.Manifold.cube([length,height,width],center=True).rotate([0,angle,0]).translate(center)
pieces=[];count=24
for i in range(count):
    a=i*2*np.pi/count;direction=np.array([np.cos(a),0,np.sin(a)])
    pieces.append(bar(direction*.043,direction*.48,.014,.60))
    # Alternating anastomoses retain spaces between adjacent plates. They are
    # not a solid lobular envelope or a collagen septum.
    for radius in [.18,.33]:
        if (i+(radius>.2))%2==0:
            b=(i+1)*2*np.pi/count;pieces.append(bar(direction*radius,np.array([np.cos(b),0,np.sin(b)])*radius,.014,.60))
lobule=md.Manifold.batch_boolean(pieces,md.OpType.Add);mesh=lobule.to_mesh64()
surface('hepatocyte-plate-network',np.asarray(mesh.vert_properties)[:,:3],np.asarray(mesh.tri_verts),[('lobule',1,None)])
v=np.array(points)*.001;origin=np.zeros(3)
# Shared positional controls are identified only when exactly equal.
unique,inv=np.unique(v,axis=0,return_inverse=True);surfaces={sid:inv[f]for sid,f in surfaces.items()}
ids={name:'TA98LIVER:micro:'+name for name in names}
words={}
for sid,f in surfaces.items():
    q=unique[f];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
    words[sid]=np.repeat(np.rint(n*32767).astype('<i2')[:,None],3,axis=1)
result=pack(unique,surfaces,cells,ids,'liver:illustrative-human-lobule:metres',out/'packed',origin=origin,corner_words=words,
    max_pairs=2000000,max_tests=16000000,checkpoint_directory=out/'checkpoints')
atomic_json(out/'evidence.json',dict(passed=result['passed'],source='https://histologyguide.com/slideview/MH-126b-liver/15-slide-1.html',
    sourceType='Human histological architecture; independently generated illustrative three-dimensional unit',
    dimensionalCalibration='Illustrative: approximately 1.3 mm across; not patient measurement',
    features=names,oneCellThickPlateWidthMicrometres=14,noThickInterlobularSepta=True,
    sinusoidalSpacesOpen=True,vascularAndBiliaryLuminaOpen=True,newTextures=0,bytes=result['total_bytes'],elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(passed=result['passed'],failures=result['checks']['failures'],bytes=result['total_bytes'],elapsedSeconds=time.monotonic()-start)),flush=True)
