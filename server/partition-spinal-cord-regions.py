"""Disjoint macroscopic horn and funiculus territories in the native cord.

Classifier volumes are construction data. Only their intersections with the
already closed grey/white cells become anatomy. No additional overlapping coat.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import manifold3d as md
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.boolean_cells import BooleanCells

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic()
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
source=base/'female-core-continuous-v3';b,stored=BooleanCells.load(source/'boolean-source.npz')
# A fresh process starts Manifold's original-ID allocator at zero. Reserve the
# loaded provenance range before adding classifier/derived primitives.
md.Manifold.reserve_ids(max(b.primitives)+1)
grey=b.derived_source(stored['grey'],'reviewed-grey-core',retain_parents=True)
white=b.derived_source(stored['white'],'reviewed-white-core',retain_parents=True)
domain=np.load(base/'female-domain-v1/domain.npz');sample_y=domain['sample_y'];dims=domain['sample_dimensions']
ys=np.r_[1.13,np.linspace(sample_y[0],sample_y[-1],72),1.58]
xs=np.linspace(-.02,.02,17)

def mask(name,factor,wedge,above):
    rings=[]
    for y in ys:
        cx,cz,rx,rz=[np.interp(y,sample_y,dims[:,i])for i in range(4)]
        height=cz+rz*(factor+wedge*np.abs((xs-cx)/rx))
        if above:poly=np.r_[np.c_[xs,height],[[xs[-1],.06],[xs[0],.06]]]
        else:poly=np.r_[[[xs[0],-.16],[xs[-1],-.16]],np.c_[xs[::-1],height[::-1]]]
        rings.append(np.c_[poly[:,0],np.full(len(poly),y),poly[:,1]])
    p=np.concatenate(rings);n=len(rings[0]);ff=[]
    for k in range(len(rings)-1):
        for i in range(n):
            j=(i+1)%n;a=k*n+i;c=(k+1)*n+j
            ff.extend([[a,k*n+j,c],[a,c,(k+1)*n+i]])
    for k,side in [(0,-1),(len(rings)-1,1)]:
        xy=rings[k][:,[0,2]];area=np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1));order=np.arange(n)if area>0 else np.arange(n-1,-1,-1)
        cap=order[np.asarray(md.triangulate([xy[order]],epsilon=0.,allow_convex=False),int)]
        q=xy[cap];cross=(q[:,1,0]-q[:,0,0])*(q[:,2,1]-q[:,0,1])-(q[:,1,1]-q[:,0,1])*(q[:,2,0]-q[:,0,0]);flip=cross*area*side<0;cap[flip]=cap[flip][:,[0,2,1]];ff.extend(cap+k*n)
    ff=np.asarray(ff,int);t=topology(Mesh(p,ff,'spinal-region-classifier'))
    assert t['closed_oriented_manifold']
    if t['oriented_triangle_integral_m3']<0:ff=ff[:,[0,2,1]]
    return b.source(p,ff,name)

ga=mask('anterior-horn-boundary',.15,0.,True)
gp=mask('posterior-horn-boundary',-.18,0.,False)
wp=mask('posterior-funiculus-boundary',-.16,-.9,False)
wa=mask('anterior-funiculus-boundary',.13,.9,True)
macro={'grey-anterior':grey^ga,'grey-posterior':grey^gp,'grey-intermediate':grey-ga-gp,
    'white-anterior':white^wa,'white-posterior':white^wp,'white-lateral':white-wa-wp}
# Native female left-root coordinates are positive X. The fixed midsagittal
# plane is shared by every tissue and is not a camera transform.
left=b.source(np.array([[0.,1.1,-.2],[.03,1.1,-.2],[.03,1.6,-.2],[0.,1.6,-.2],[0.,1.1,.1],[.03,1.1,.1],[.03,1.6,.1],[0.,1.6,.1]]),
    np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]]),'fixed-midsagittal-classifier')
solids={}
for name,s in macro.items():
    solids[name+'-left']=s^left;solids[name+'-right']=s-left
for name,s in solids.items():assert str(s.status())=='Error.NoError'and not s.is_empty(),name
before=sum(float(s.volume())for s in stored.values());after=sum(float(s.volume())for s in solids.values())
assert abs(before-after)<1e-12
b.save(solids,out/'boolean-source.npz');result=b.assemble(solids)
atomic_npz(out/'cells.npz',positions=result['positions'],**result['surfaces'])
atomic_json(out/'cells.json',dict(cells=result['cells'],surfaces={k:len(fs)for k,fs in result['surfaces'].items()},
    construction=result['evidence'],source_sha256=digest_file(source/'cells.npz'),
    boundary_classification='Macroscopic regional horn/funiculus partitions; cytoarchitectural nuclei, laminae and detailed tracts are not implied.',
    source_envelope_preserved=True,canal_representation='Continuous lumen remains a void',
    volumes_mL={k:float(s.volume())*1e6 for k,s in solids.items()},whole_volume_difference_m3=after-before,
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(parts=len(solids),surfaces=len(result['surfaces']),volume_difference_m3=after-before,elapsed_seconds=time.monotonic()-start)),flush=True)
