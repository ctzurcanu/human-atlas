"""Bounded source-edge handle discovery; source coordinates remain literal."""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.algorithms.passage_cycles import handle_cycle

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'

def main():
    began=time.monotonic();level=sys.argv[1];assert level in ['C1','C2','C6']
    out=BASE/('actual-handle-cycles-independent-v2-'+level);out.mkdir(exist_ok=False)
    inv=json.loads((BASE/'source-inventory.json').read_bytes());oldpath=BASE/('actual-aperture-discovery-independent-v1-'+level)/'evidence.json';old=json.loads(oldpath.read_bytes());assert old['analysisAssemblyTopology']['boundaryEdges']==0 and old['analysisAssemblyTopology']['nonmanifoldEdges']==0
    vv=[];ff=[]
    for pid in old['analysisAssemblyTopology']['materialPieces']:
        row=inv['sourceBindings'][pid];z=np.load(row['mesh']);ff.append(z['faces']+sum(len(x)for x in vv));vv.append(z['vertices'])
    vertices=np.vstack(vv);v,ix=np.unique(vertices,axis=0,return_inverse=True);f=ix[np.vstack(ff)];rows=[]
    for side in 'lr':
        st=np.load(BASE/(side+'-literal-source-stations.npz'));c=st['centers'];station={'C1':90,'C2':70,'C6':25}[level];seed=c[station].copy();method='Native artery point is only a conduit-search seed; topological cycle must independently establish a passage.'
        if level=='C1':
            contours=old['allSourceContours'];sign=1 if side=='l'else-1
            island=max((q for q in contours if sign*q['center_m'][0]>.02),key=lambda q:sign*q['center_m'][0]);main=max(contours,key=lambda q:q['area_m2']);a=np.array(island['worldPoints_m']);b=np.array(main['worldPoints_m']);dist=np.linalg.norm(a[:,None]-b[None],axis=2);ii,jj=np.unravel_index(dist.argmin(),dist.shape);seed=(a[ii]+b[jj])/2
            method='Midpoint of the closest literal probe-contour vertices between the outer transverse-process island and main atlas mass; a search seed only, not an ostium or accepted axis.'
        excluded=[np.array([-seed[0],seed[1],seed[2]]),np.array([0.,seed[1],-.03])]
        row=dict(side=side,sourceSearchSeed_m=seed,searchSeedMethod=method,sourceArteryStationForProbe=station,otherPassageAxesExcluded=excluded)
        try:
            loop,proof=handle_cycle(v,f,seed,axis=(0,1,0),half_length=.025,radius=.018,excluded_centers=excluded,maximum_starts=12,wall_seconds=4.)
            p=v[loop];row.update(physicalSourceVertexIds=loop,worldPoints_m=p,cycleCentroid_m=p.mean(0),cycleBounds_m=[p.min(0),p.max(0)],proof=proof)
        except (ValueError,TimeoutError)as ex:row['failure']=dict(kind=type(ex).__name__,reason=str(ex))
        rows.append(row)
    atomic_json(out/'evidence.json',dict(schema='native-vertebral-nonseparating-handle-reference/v2',level=level,priorLiteralSectionSha256=digest_file(oldpath),sourceInventorySha256=digest_file(BASE/'source-inventory.json'),rows=rows,
        exactOriginalPositionsAndAllFacesOnly=True,noGeometryOrPartition=True,anatomicalAcceptance=False,
        limitation='A shortest source-edge nonseparating cycle is a throat reference, not necessarily the inferior entry or superior exit lip. Empty aperture and physiological full-calibre passage require explicit section and context review.',elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(evidence=str(out/'evidence.json'),rows=[dict(side=r['side'],center_mm=[float(x*1000)for x in r['cycleCentroid_m']]if'cycleCentroid_m'in r else None,proof=r.get('proof'),failure=r.get('failure'))for r in rows],elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
