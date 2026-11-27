"""Deterministic actual-triangle surface/uncapped-section diagnostics."""
from pathlib import Path
import hashlib, json, time, sys
import numpy as np
from PIL import Image, ImageDraw

BASE=Path('.local-models/geometry-workbench/circulatory-completion/lateral-circumflex-femoral-native-v1/source-ownership-proposal-v1')

def clip(poly,normal,offset):
    out=[]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        fa=normal@a-offset; fb=normal@b-offset
        if fa>=0: out.append(a)
        if (fa>=0)!=(fb>=0):out.append(a+(b-a)*fa/(fa-fb))
    return np.asarray(out).reshape(-1,3)

def main():
    began=time.monotonic(); entries=[]; side=sys.argv[1] if len(sys.argv)>1 else 'left'
    offset=2 if side=='right' else 0
    for file,color in [(f'context-{offset}.npz',[170,170,150]),(f'context-{offset+1}.npz',[135,45,45]),(f'mesh-{offset}.npz',[180,65,60]),(f'mesh-{offset+1}.npz',[245,150,40])]:
        z=np.load(BASE/file,allow_pickle=False);entries.append((file,z['vertices'][z['faces']],np.array(color)))
    branch=np.concatenate([entries[-1][1].reshape(-1,3),entries[-2][1].reshape(-1,3)])
    low,high=branch.min(0)-.006,branch.max(0)+.006; center=(low+high)/2
    images=[]
    for name,direction,section in [('front',[0,0,1],False),('posterior',[0,0,-1],False),('lateral',[1,0,0],False),('root-section',[0,0,1],True)]:
        forward=np.array(direction,float);up=np.array([0.,1,0]);right=np.cross(up,forward);basis=np.array([right,up,forward]);corners=np.array(np.meshgrid(*zip(low,high))).reshape(3,-1).T;proj=(corners-center)@basis.T;scale=580/max(np.ptp(proj[:,:2],axis=0));polygons=[]
        for file,triangles,color in entries:
            for triangle in triangles:
                poly=triangle.copy()
                for axis in range(3):
                    n=np.eye(3)[axis];poly=clip(poly,n,low[axis])
                    if not len(poly):break
                    poly=clip(poly,-n,-high[axis])
                    if not len(poly):break
                if section and len(poly):poly=clip(poly,np.array([0.,1,0]),.884)
                if len(poly)<3:continue
                normal=np.cross(triangle[1]-triangle[0],triangle[2]-triangle[0]);normal/=max(np.linalg.norm(normal),1e-30)
                brightness=.45+.55*abs(normal@np.array([.3,.5,.81240384]));rgb=tuple(np.rint(color*brightness).astype(int));q=(poly-center)@basis.T;xy=np.c_[320+q[:,0]*scale,320-q[:,1]*scale];polygons.append((float(q[:,2].mean()),xy,rgb))
        image=Image.new('RGB',(640,640),(18,26,35));draw=ImageDraw.Draw(image)
        for depth,xy,rgb in sorted(polygons,key=lambda x:x[0]):draw.polygon([tuple(x)for x in xy],fill=rgb)
        path=BASE/((side+'-' if side=='right' else '')+name+'.png');image.save(path);images.append(dict(name=side+'-'+name,path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),section=section,actual_triangle_clip=True,inspected=False))
    proof=dict(schema='literal-lcfa-diagnostic-render/v1',views=images,source_mesh_sha256={file:hashlib.sha256((BASE/file).read_bytes()).hexdigest()for file,_,_ in entries},renderer='PIL actual triangle painter projection; diagnostic only, not photorealism or patency proof',uncapped_section=True,elapsed_seconds=time.monotonic()-began,anatomical_acceptance=False)
    (BASE/((side+'-' if side=='right' else '')+'render-proof.json')).write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))

if __name__=='__main__':main()
