"""Read-only vertebral regional boundary semantics and literal context witnesses."""
from pathlib import Path
import json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.distance import barycentric
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion'
OUT=BASE/'vertebral-boundary-semantics-independent-v1'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME='model:local-ta98:assembled-metres'

def segment_hits(a,b,tree):
    ids=tree.query_bounds(np.minimum(a,b),np.maximum(a,b));t=tree.triangles[ids]
    if not len(t):return []
    n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);da=np.einsum('ij,ij->i',n,a-t[:,0]);db=np.einsum('ij,ij->i',n,b-t[:,0]);take=da*db<0
    ids=ids[take];t=t[take];mu=-da[take]/(db[take]-da[take]);p=a+mu[:,None]*(b-a);w,valid=barycentric(p,t);keep=valid&np.all(w>=-1e-10,axis=1)
    return [dict(sourceFace=int(i),segmentParameter=float(u),point_m=x.tolist(),barycentric=ww.tolist())for i,u,x,ww in zip(ids[keep],mu[keep],p[keep],w[keep])]

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False)
    proof=BASE/'arterial-source-opportunities-v1/context-qualified-scan-v2/proof.json';old=json.loads(proof.read_bytes());prior={q['part_id']:q for q in old['nativeCourseMeshes']}
    sha=digest_file(CAT);atlas=json.loads(CAT.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();meshes={};trees={};bindings={}
    ids=['ZA:Vertebral artery.l','ZA:Vertebral artery.r','ZA:Vertebra C6:Bone-2','ZA:Axis (C2):Bone-2','ZA:Atlas (C1):Bone-2','ZA:Occipital bone:Bone-6','ZA:Spinal dura','LOCAL:male:posterior_atlanto_occipital_membrane:695']
    for pid in ids:
        m=reader.read_part(atlas,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha);meshes[pid]=m;trees[pid]=BVH(m)
        bindings[pid]=dict(geometrySha256=m.geometry_sha,triangles=len(m.faces),bounds_m=[m.vertices.min(0).tolist(),m.vertices.max(0).tolist()],sourceLicense=parts[pid]['license'])
        if pid in prior:assert m.geometry_sha==prior[pid]['geometry_sha256'];bindings[pid]['existingOpportunityGeometryExactlyReused']=True
    rows=[]
    for side in 'lr':
        pid='ZA:Vertebral artery.'+side;m=meshes[pid];v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];q=recover_station_rings(v,f);assert len(q)==1
        rr=list(q[0]['rings'])
        if v[rr[0]].mean(0)[1]>v[rr[-1]].mean(0)[1]:rr=rr[::-1]
        rr[0]=next(loop for loop in boundary_loops(f)if set(loop)==set(rr[0]));rings,_=align_station_correspondence(v,f,rr);c=np.array([v[r].mean(0)for r in rings])
        contexts=[]
        for target in ids[2:]:
            mm=meshes[target];selected=np.flatnonzero((c[:,1]>=mm.vertices[:,1].min()-.005)&(c[:,1]<=mm.vertices[:,1].max()+.005))
            witnesses=[]
            for i in selected:
                w=trees[target].closest_point(c[i]);witnesses.append(dict(station=int(i),center_m=c[i].tolist(),sourceSurfaceWitness=w))
            hits=[]
            for i,(a,b)in enumerate(zip(c[:-1],c[1:])):
                for h in segment_hits(a,b,trees[target]):hits.append(dict(sourceStationInterval=[i,i+1],**h))
            contexts.append(dict(partId=target,centerlineSurfaceCrossings=hits,nearLevelStationWitnesses=witnesses,
                neitherNearestDistanceNorWholeBoneHeightDefinesAnOstium=True))
        rows.append(dict(side=side,sourcePartId=pid,sourceStations=len(rings),stationCenters_m=c.tolist(),context=contexts,
            prevertebralEndNotChosen=True,atlanticStartNotChosen=True,intracranialStartNotChosen=True))
    refs=[dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11366551/',doi='10.31616/asj.2024.0012',study='2024 primary dissection of40 vertebral arteries from20 cadavers',
        finding='Study atlantic/V3 includes a vertical part from C2 transverse-foramen exit to C1 transverse-foramen exit, a horizontal atlas-groove part and an exit part reaching the dural entry. C1exit is an internal sub-boundary under this convention.'),
      dict(url='https://doi.org/10.1016/j.jocn.2008.05.033',study='2009 primary dissection of7 cadavers and CT comparison in24 individuals',
        finding='Uses an alternative classical cervical part through C6–C1 and a third part from the C1 foramen to dura. Do not silently mix this convention with C2-based atlantic/V3.'),
      dict(url='https://pubmed.ncbi.nlm.nih.gov/31756706/',doi='10.3171/2019.9.JNS191394',study='6 preserved adult cadaveric dural-entry sites with histology',
        finding='Dural connective tissue and arterial adventitial connective tissue have no clear intervening plane. Actual dural entry is a tissue interface; a bony foramen-magnum plane does not identify that interface.'),
      dict(url='https://pubmed.ncbi.nlm.nih.gov/25031813/',doi='10.3344/kjp.2014.27.3.266',study='CTA review of920 arteries in460 patients',
        finding='C6 was the entry level in855 arteries, about93%; other levels occur. A generic C6 horizontal plane is not a substitute for actual entry into the modeled transverse-foramen aperture.')]
    atomic_json(OUT/'evidence.json',dict(schema='vertebral-native-boundary-semantics-independent/v1',sourceCatalogueSha256=sha,existingOpportunityProofSha256=digest_file(proof),sourceBindings=bindings,sides=rows,primaryReferences=refs,
        recommendedRecordedConvention='Use C2exit as cervical-to-atlantic boundary if adopting the2024 atlantic/V3 convention; retain C1exit separately as atlantic internal anatomical landmark. Dural entry is required for V3/V4 under either convention.',
        localPrivateMembraneRejectedAsBoundaryGuide='Its unresolved-frame Y1.70254–1.71546m lies wholly above this native artery maximumY1.69089m.',
        noGeometryOrCatalogueChanged=True,noPartitionChosen=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),sides=[dict(side=r['side'],stations=r['sourceStations'],contextHits=[dict(partId=x['partId'],hits=len(x['centerlineSurfaceCrossings']),intervals=[h['sourceStationInterval']for h in x['centerlineSurfaceCrossings']])for x in r['context']])for r in rows],elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
