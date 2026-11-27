"""Complete actual normal-loader context audit of unchanged vertebral sources."""
from pathlib import Path
import json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'
OUT=BASE/'complete-current-context-v3'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    inventory=json.loads((BASE/'source-inventory.json').read_bytes());runtime=json.loads((BASE/'actual-runtime-context-v2.json').read_bytes());atlas=json.loads(CAT.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    source={s:Mesh.load(inventory['sourceBindings']['ZA:Vertebral artery.'+s]['mesh'])for s in ['l','r']};trees={s:BVH(m)for s,m in source.items()};rows=[]
    for index,p in enumerate(runtime['rows']):
        if not p['normalAvailable']:continue
        pid=p['id'];record=parts[pid];m=reader.read_part(atlas,record,frame='model:local-ta98:assembled-metres',catalogue_path=CAT,catalogue_sha=EXPECTED)
        row=dict(partId=pid,system=p['system'],normalRuntimeAvailable=True,recordSha256=identity(record),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'],wholeSourceBounds_m=m.bounds.tolist(),wholeSourceTriangles=len(m.faces),checks=[])
        if pid.startswith('ZA:Vertebral artery.'):
            row['role']='source artery; whole-source self audit already complete';rows.append(row);continue
        tri=m.triangles
        for s,artery in source.items():
            lo,hi=artery.bounds;faceIds=np.flatnonzero(np.all(tri.max(1)>=lo,axis=1)&np.all(tri.min(1)<=hi,axis=1))
            if not len(faceIds):row['checks'].append(dict(side=s,completeSourceTriangleAABBsDisjoint=True,crossingPairs=0));continue
            clipped=Mesh(m.vertices,m.faces[faceIds],m.frame,m.attributes,dict(m.source,actualSourceFaceIds=faceIds.tolist(),broadPhaseOnly=True));path=OUT/(identity(pid+':'+s)[:20]+'.npz');clipped.save(path)
            details,gate=pair_contacts(artery,clipped,0.,first_bvh=trees[s],max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget']
            row['checks'].append(dict(side=s,actualCorridorSnapshot=str(path),snapshotSha256=digest_file(path),sourceFaceIds=faceIds.tolist(),gate=gate,
                defects=[dict(arterySourceFace=int(a),fixedSourceFace=int(faceIds[b]),relation=str(rel))for(a,b),rel in zip(details['pairs'],details['relation'])if rel in ['crossing','coplanar-overlap']],crossingPairs=gate['crossing_face_pairs']))
        rows.append(row)
        if len(rows)%100==0:atomic_json(OUT/('progress-'+str(len(rows))+'.json'),dict(processed=len(rows),elapsedSeconds=time.monotonic()-began))
    result=dict(schema='vertebral-whole-current-source-actual-context/v3',sourceCatalogueSha256=EXPECTED,sourceInventorySha256=digest_file(BASE/'source-inventory.json'),runtimeContextSha256=digest_file(BASE/'actual-runtime-context-v2.json'),
        actualNormalSourcesDecoded=len(rows),allSystemsIncluded=True,wholeSourceCourseAndSelfAuditReused=True,allActualContactsComplete=True,rows=rows,
        noGeometryOrVisibilityChanged=True,contactDoesNotEstablishAnatomicalInterface=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',result);assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),actualSources=len(rows),elapsedSeconds=result['elapsedSeconds'],crossings=[dict(id=q['partId'],system=q['system'],side=z['side'],pairs=z['crossingPairs'])for q in rows for z in q['checks']if z['crossingPairs']])))

if __name__=='__main__':main()
