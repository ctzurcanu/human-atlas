"""Saved baseline against every actual retained local physical context mesh."""
from pathlib import Path
import json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
OUT=BASE/'left-whole-context-baseline-v24'

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False)
    contract=json.loads((BASE/'course-contract-v23/contract.json').read_bytes());runtime=json.loads(Path(contract['fixedCurrentContext']).read_bytes())
    spec=next(r for r in contract['regions']if r['side']=='l'and r['kind']=='lateral')
    source=Mesh.load(contract['sourceBindings'][spec['sourcePartId']]['mesh']);ids=np.array(spec['exactSourceFaceOwnership'],int)
    region=Mesh(source.vertices,source.faces[ids],source.frame);tree=BVH(region);box=region.bounds;rows=[]
    for q in runtime['rows'][0]['rows']:
        if not q['retained']or q['id']==spec['sourcePartId']:continue
        binding=contract['sourceBindings'][q['id']];assert digest_file(binding['mesh'])==binding['meshSha256'];m=Mesh.load(binding['mesh']);bb=m.bounds
        if np.any(bb[1]<box[0])or np.any(bb[0]>box[1]):
            rows.append(dict(partId=q['id'],system=q['system'],broadPhaseSeparate=True,gate=None));continue
        data,gate=pair_contacts(region,m,0.,first_bvh=tree,max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget']
        rows.append(dict(partId=q['id'],system=q['system'],gate=gate,sourceFacePairs=data['pairs'].tolist(),relations=data['relation'].tolist(),broadPhaseSeparate=False))
    result=dict(schema='left-lateral-foot-complete-current-context-baseline/v24',contractSha256=digest_file(BASE/'course-contract-v23/contract.json'),sourceCatalogueSha256=contract['sourceCatalogueSha256'],
        rows=rows,sourceRegionFaces=len(ids),physicalMeshCount=len(rows),anatomicalAcceptance=False,noChanges=True,elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'evidence.json',result)
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),physicalMeshes=len(rows),elapsedSeconds=result['elapsedSeconds'],
        contacts=[dict(id=r['partId'],cross=r['gate']['crossing_face_pairs'],touch=r['gate']['touching_face_pairs'],overlap=r['gate']['coplanar_overlap_face_pairs'])for r in rows if r['gate']and r['gate']['returned_triangle_pairs']])))

if __name__=='__main__':main()
