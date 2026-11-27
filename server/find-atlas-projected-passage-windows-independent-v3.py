"""Find empty through-hole projection windows, then verify source handle cycles."""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.ndimage import label,distance_transform_edt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.passage_cycles import handle_cycle
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'

def main():
    began=time.monotonic();out=BASE/'atlas-projected-passage-windows-independent-v3';out.mkdir(exist_ok=False)
    inv=json.loads((BASE/'source-inventory.json').read_bytes());vs=[];fs=[]
    for pid in ['ZA:Atlas (C1):Bone-2','ZA:Atlas (C1):Cartilage']:
        row=inv['sourceBindings'][pid];assert digest_file(row['mesh'])==row['sha256'];z=np.load(row['mesh']);fs.append(z['faces']+sum(len(x)for x in vs));vs.append(z['vertices'])
    vertices=np.vstack(vs);v,ix=np.unique(vertices,axis=0,return_inverse=True);f=ix[np.vstack(fs)];t=v[f][:,:,[0,2]];spacing=.00025
    lo=v[:,[0,2]].min(0)-.002;hi=v[:,[0,2]].max(0)+.002;xx=np.arange(lo[0],hi[0]+spacing,spacing);zz=np.arange(lo[1],hi[1]+spacing,spacing);covered=np.zeros((len(xx),len(zz)),bool)
    for tri in t:
        a,b,c=tri;e=b-a;g=c-a;det=e[0]*g[1]-e[1]*g[0]
        if abs(det)<1e-16:continue
        ii=np.flatnonzero((xx>=tri[:,0].min()-1e-12)&(xx<=tri[:,0].max()+1e-12));jj=np.flatnonzero((zz>=tri[:,1].min()-1e-12)&(zz<=tri[:,1].max()+1e-12))
        if not len(ii)or not len(jj):continue
        x,z=np.meshgrid(xx[ii],zz[jj],indexing='ij');q=np.stack([x-a[0],z-a[1]],-1);u=(q[:,:,0]*g[1]-q[:,:,1]*g[0])/det;w=(e[0]*q[:,:,1]-e[1]*q[:,:,0])/det;covered[np.ix_(ii,jj)]|=(u>=-1e-12)&(w>=-1e-12)&(u+w<=1+1e-12)
    ll,count=label(~covered);border=set(np.r_[ll[0],ll[-1],ll[:,0],ll[:,-1]].tolist());rows=[];y=float(np.load(BASE/'l-literal-source-stations.npz')['centers'][90,1])
    for k in range(1,count+1):
        if k in border:continue
        mask=ll==k;dist=distance_transform_edt(mask)*spacing;i,j=np.unravel_index(dist.argmax(),dist.shape);seed=np.array([xx[i],y,zz[j]]);row=dict(component=k,gridSamples=int(mask.sum()),sourceAxis=[0,1,0],seed_m=seed,approximateRasterInradius_m=float(dist.max()),gridSpacing_m=spacing)
        try:
            cyc,proof=handle_cycle(v,f,seed,axis=(0,1,0),half_length=.025,radius=.018 if abs(seed[0])>.008 else .04,
                excluded_centers=[[-seed[0],seed[1],seed[2]],[0.,seed[1],-.03]]if abs(seed[0])>.008 else[],maximum_starts=12,wall_seconds=4.)
            row.update(sourceVertexIds=cyc,sourceCycleWorldPoints_m=v[cyc],sourceCycleCentroid_m=v[cyc].mean(0),proof=proof)
        except(ValueError,TimeoutError)as ex:row['failure']=dict(kind=type(ex).__name__,reason=str(ex))
        rows.append(row)
    atomic_npz(out/'projection-grid.npz',x_m=xx,z_m=zz,covered_by_source_triangle=covered,empty_components=ll)
    atomic_json(out/'evidence.json',dict(schema='atlas-native-projected-through-hole-windows/v3',sourceInventorySha256=digest_file(BASE/'source-inventory.json'),sourceCatalogueSha256=inv['sourceCatalogueSha256'],rows=rows,
        exactSourcePositionAndFaceAnalysis=True,noGeometryOrPartitionChanged=True,
        limitation='Raster clearance only proposes a source empty-window seed. Nonseparating source-edge cycle is a separate topological proof; neither establishes an ostial boundary or full-calibre passage.',elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(evidence=str(out/'evidence.json'),rows=[dict(component=r['component'],samples=r['gridSamples'],seed_mm=(r['seed_m']*1000).tolist(),inradius_mm=r['approximateRasterInradius_m']*1000,proof=r.get('proof'),failure=r.get('failure'))for r in rows],elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
