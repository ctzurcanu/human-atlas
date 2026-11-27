"""Install the user's explicitly requested provisional bilateral refinement.

Preparation: python -m geometry_tools.commands.install_tibial_source_preview
Apply the prepared, replayed local revision: add --apply.
This never installs the rejected fitted courses, changes remote deployments,
or claims anatomical acceptance of the remaining recipient contacts.
"""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from geometry_tools.common import ROOT
from geometry_tools.course_io import sha,load_arrays,write_json
from geometry_tools.install import catalogue_diff,transaction
from geometry_tools.model import AtlasReader

B=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
SOURCE=B/'bilateral-tibial-shaft-v85'
BASE=B/'atlas-cardiovascular-review-v31/atlas.json'
OUT=B/'atlas-cardiovascular-review-v32'
ROUTE=ROOT/'app/local-ta98-layer-catalogue.ts'


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--apply',action='store_true')
    args=parser.parse_args();started=time.monotonic()
    if args.apply:
        m=json.loads((OUT/'merge-report.json').read_text())
        for path,digest in [(BASE,m['baselineSha256']),(OUT/'atlas.json',m['catalogueSha256']),
                            (OUT/'geometry.bin',m['bankSha256']),(ROUTE,m['routeSourceSha256'])]:
            if sha(path)!=digest:raise ValueError('Prepared installation changed: '+str(path))
        old=ROUTE.read_bytes();needle=b'atlas-cardiovascular-review-v31/atlas.json'
        if old.count(needle)!=3:raise ValueError('Unexpected active local routes')
        new=old.replace(needle,b'atlas-cardiovascular-review-v32/atlas.json')
        result=transaction({ROUTE:new},OUT/'rollback-route',expected={str(ROUTE):sha(ROUTE)})
        m.update(installed=True,transaction=result,installedRouteSha256=sha(ROUTE))
        write_json(OUT/'merge-report.json',m)
        print(json.dumps(dict(installed=True,changedParts=len(m['changedParts']),
                             catalogueSha256=m['catalogueSha256'],elapsedSeconds=time.monotonic()-started)))
        return
    if OUT.exists():raise ValueError('Prepared local revision already exists')
    r=json.loads((SOURCE/'report.json').read_text())
    replay=json.loads((SOURCE/'independent-replay.json').read_text())
    review=json.loads((SOURCE/'visual-review-v95.json').read_text())
    gate_path=B/'tibial-install-source-gate-v396.json';gate=json.loads(gate_path.read_text())
    if r['sourceProfileSha256']!=sha(BASE) or replay['candidateReportSha256']!=sha(SOURCE/'report.json'):
        raise ValueError('Source/refinement replay is stale')
    if not gate['allSelfChecksPassed'] or r['negativeTransportedFaceOrientations'] or not r['topologyPreserved']:
        raise ValueError('Source-preserving refinement failed an intrinsic geometry check')
    if not replay['originalPosedControlsRetained'] or replay['unmatchedMovingMouths']:
        raise ValueError('Original controls or source mouths were not retained')
    if review['candidateReportSha256']!=sha(SOURCE/'report.json'):raise ValueError('Stale rendered review')
    for view in review['views']:
        if sha(view['path'])!=view['sha256']:raise ValueError('Reviewed image changed')
    before=json.loads(BASE.read_text());after=copy.deepcopy(before);by={p['id']:p for p in after['parts']}
    bindings={x['file']:x['sha256']for x in r['candidateArrayBindings']}
    gates={x['partId']:x for x in gate['rows']};bank=bytearray();changed=[];arrays={};chunk=len(after['chunks'])
    def pack(value,dtype):
        while len(bank)%4:bank.append(0)
        offset=len(bank);bank.extend(np.asarray(value,dtype=dtype).tobytes());return offset
    for row in r['rows']:
        pid=row['partId'];path=SOURCE/row['geometryFile']
        if sha(path)!=bindings[path.name] or sha(path)!=gates[pid]['sha256'] or by[pid]!=row['sourceRecord']:
            raise ValueError('Refinement source record changed')
        a=load_arrays(path);arrays[pid]=a;p=copy.deepcopy(by[pid])
        p.update(chunk=chunk,positions=pack(a['local_positions'],'<f4'),
                 normals=pack(a['normals'],'<i2'),normalSigns=pack(a['normalSigns'],'i1'),
                 indices=pack(a['faces'],'<u4'),vertexCount=len(a['vertices']),indexCount=int(a['faces'].size),
                 bounds=[a['vertices'].min(0).tolist(),a['vertices'].max(0).tolist()])
        if 'uvs'in a:p['uvs']=pack(a['uvs'],'<f4')
        elif 'uvs'in p:raise ValueError('Refinement unexpectedly removed source UVs')
        by[pid]=p;changed.append(pid)
    OUT.mkdir();(OUT/'geometry.bin').write_bytes(bank)
    bankinfo=dict(url='/local-models/'+str((OUT/'geometry.bin').relative_to(ROOT/'.local-models')),
                  bytes=len(bank),sha256=sha(OUT/'geometry.bin'))
    after['chunks'].append(bankinfo)
    # This map is indexed by the CURRENT chunk and returns its canonical source
    # identity. Preserve every used slot; discard only stale unused tail slots.
    source_map=before.get('sourceChunkIndices',list(range(chunk)))
    if len(source_map)<chunk:raise ValueError('Incomplete existing canonical bank map')
    canonical=before.get('sourceChunkCount',max(source_map)+1)
    after['sourceChunkIndices']=source_map[:chunk]+[canonical]
    after['sourceChunkCount']=canonical+1
    after['parts']=[by[p['id']]for p in before['parts']]
    if isinstance(after.get('triangles'),int):after['triangles']+=r['newTriangles']
    limits=['Provisional source-preserving refinement installed at the user request; not anatomical acceptance.',
            'Existing muscle, tendon, bone, nerve and other recipient contacts remain for later polish.',
            'Rejected independent-profile fitting fields and recipient muscle trials are not installed.']
    after['biodigitalReview']['tibialSourceRefinement']=dict(partIds=changed,sourceCandidate='bilateral-tibial-shaft-v85',
        installedForLocalPreview=True,provisional=True,anatomicalAcceptance=False,newTerms=0,remainingLimits=limits)
    catalogue_diff(before,after,changed)
    for key in ['concepts','materials','ta98TermBindings','defaultHidden','sourceCardiacIntegration']:
        if key in before and before[key]!=after[key]:raise ValueError('Unrelated metadata changed: '+key)
    if any(after['sourceChunkIndices'][p['chunk']]!=source_map[p['chunk']]for p in before['parts']if p['id']not in changed):
        raise ValueError('Unchanged carrier canonical bank identity changed')
    write_json(OUT/'atlas.json',after)
    reader=AtlasReader();checks=[]
    for pid in changed:
        a=arrays[pid];mesh=reader.read_part(after,by[pid],frame='atlas-world-metres',catalogue_path=OUT/'atlas.json')
        used,inv=np.unique(a['faces'],return_inverse=True)
        if not np.array_equal(mesh.vertices,a['vertices'][used]) or not np.array_equal(mesh.faces,inv.reshape(-1,3)):
            raise ValueError('Emitted bank coordinates/faces failed exact replay')
        for name in ['normals','normalSigns','uvs']:
            if name in a and not np.array_equal(mesh.attributes[name],a[name][used]):raise ValueError('Stored attribute mismatch: '+name)
        checks.append(dict(partId=pid,positionsFacesAndAttributesLiteral=True))
    write_json(OUT/'packed-bank-replay.json',dict(rows=checks,allOtherRecordsLiteral=True,duplicatePartsAdded=0))
    report=dict(schema='provisional-bilateral-tibial-source-install/v32',baselineSha256=sha(BASE),
        catalogueSha256=sha(OUT/'atlas.json'),bankSha256=sha(OUT/'geometry.bin'),
        candidateReportSha256=sha(SOURCE/'report.json'),selfGateSha256=sha(gate_path),
        priorReviewSha256=sha(SOURCE/'visual-review-v95.json'),routeSourceSha256=sha(ROUTE),
        changedParts=changed,sourceProfilesPreserved=6130,newProfiles=r['newProfiles'],
        newTriangles=r['newTriangles'],duplicatePartsAdded=0,skeletonSkinFeetHeartAndOtherPartsLiteral=True,
        userDecision='Install what you have now and defer tibial polish; hand fingers next.',
        remainingLimits=limits,provisional=True,installed=False,anatomicalAcceptance=False,
        elapsedSeconds=time.monotonic()-started)
    write_json(OUT/'merge-report.json',report)
    print(json.dumps({k:v for k,v in report.items()if k!='changedParts'}))


if __name__=='__main__':main()
