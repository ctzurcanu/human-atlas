"""Expand only the necessary native basilar dependency; inventory full authority."""
from pathlib import Path
import json,sqlite3,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file,identity
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'
OUT=BASE/'necessary-basilar-source-closure-v5'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    authority=ROOT/'.local-models/ta98-authority.sqlite';db=sqlite3.connect(authority);db.row_factory=sqlite3.Row
    terms=[dict(q)for q in db.execute('WITH RECURSIVE t(id) AS (VALUES(?) UNION ALL SELECT a.id FROM ta98 a JOIN t ON a.parent_id=t.id) SELECT a.* FROM ta98 a JOIN t ON a.id=t.id ORDER BY a.id',('A12.2.07.081',))];assert len(terms)==40
    censuspath=ROOT/'reports/ta98-blood-circulation-inventory-installed-v12.json';census=json.loads(censuspath.read_bytes());ct={q['id']:q for q in census['terms']};atlas=json.loads(CAT.read_bytes());parts={q['id']:q for q in atlas['parts']}
    for q in terms:
        current=ct.get(q['id'],{});q.update(existingSourcePartIds=current.get('active_exact_source_part_ids',[]),representationStatus=current.get('representation_status','missing-direct-representation'),
            plannedRepresentationClass='named macroscopic arterial conduit/branch/network, with no new wall or microscopic inference',geometryNotConstructed=True)
    atomic_json(OUT/'authority-inventory.json',dict(schema='complete-basilar-applicable-descendants-before-source-correction/v5',sourceCatalogueSha256=EXPECTED,authoritySha256=digest_file(authority),censusSha256=digest_file(censuspath),terms=terms))
    # Required checked-in authority report is written before any candidate.
    report=ROOT/'reports/mesh-audit/REVIEW.md';heading='## Necessary proximal basilar dependency: complete authority inventory v5'
    assert heading not in report.read_text()
    lines=['',heading,'',
        'The literal vertebral/basilar junction crosses the fixed occipital source: each vertebral terminal ring has two actual ring-edge/triangle hits, and 11 of 12 points of both vertebral terminals and the basilar proximal ring are bracketed by posterior/anterior occipital source faces. The occipital source is not a certified closed solid; direct crossing and slice witnesses remain explicit. The basilar distal ring has zero such AP-bracketed points. Root authorized proximal basilar dependency expansion; source skeleton and installed heart remain fixed. No named branch is inferred from the existing unbranched carrier.',
        '', '| TA98 | Parent | Structure | Current direct source/reference IDs | Status |','|---|---|---|---|---|']
    for q in terms:lines.append('| '+q['id']+' | '+str(q['parent_id'])+' | '+q['name_en']+' | '+(', '.join(q['existingSourcePartIds'])or'—')+' | '+q['representationStatus']+' |')
    lines+=['', 'All 40 applicable basilar descendants are recorded, including its posterior-cerebral and labyrinthine subtrees. Only the proximal basilar source may join the moving closure as required by actual interfaces. Existing AICA/pontine/PICA/anterior-spinal source origins are protected outside explicit support or added as specifically measured dependency roots; distant arterial courses do not move by generic containment. Named branches without independent geometry remain missing. New regions must cover the repaired whole vertebral source once, with shared literal boundaries and no duplicated whole carrier faces.', '']
    with report.open('a')as f:f.write('\n'.join(lines))
    reader=AtlasReader();bindings={};courses=[]
    ids=sorted({pid for q in terms for pid in q['existingSourcePartIds']if pid.startswith('ZA:')})
    ids+=['ZA:Anterior inferior cerebellar artery.l','ZA:Anterior inferior cerebellar artery.r','ZA:Medial pontine branches of basilar artery.l','ZA:Medial pontine branches of basilar artery.r','ZA:Lateral pontine branches of basilar artery.l','ZA:Lateral pontine branches of basilar artery.r']
    for pid in dict.fromkeys(ids):
        p=parts[pid];m=reader.read_part(atlas,p,frame='model:local-ta98:assembled-metres',catalogue_path=CAT,catalogue_sha=EXPECTED);path=OUT/(identity(pid)[:20]+'.npz');m.save(path)
        bindings[pid]=dict(partId=pid,mesh=str(path),sha256=digest_file(path),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,recordSha256=identity(p),sourceBankSha256=m.source['bank_sha256'],record=p)
        if pid!='ZA:Basilar artery':continue
        v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];rr=recover_station_rings(v,f);assert len(rr)==1;rings=list(rr[0]['rings'])
        if v[rings[0]].mean(0)[1]>v[rings[-1]].mean(0)[1]:rings=rings[::-1]
        rings,phase=align_station_correspondence(v,f,rings);rings=np.array(rings);sp=OUT/'basilar-literal-source-stations.npz';atomic_npz(sp,physical_vertices=v,physical_faces=f,source_to_physical=inv,rings=rings,centers=v[rings].mean(1),source_faces=m.faces)
        courses.append(dict(partId=pid,stations=len(rings),stationSnapshot=str(sp),stationSnapshotSha256=digest_file(sp),sourceHasNoActualBranchSockets=True))
    result=dict(schema='necessary-vertebrobasilar-proximal-source-moving-closure/v5',sourceCatalogueSha256=EXPECTED,authorityInventorySha256=digest_file(OUT/'authority-inventory.json'),sourceBindings=bindings,courses=courses,
        minimumMovingOwners=['ZA:Vertebral artery.l','ZA:Vertebral artery.r','ZA:Basilar artery'],
        conditionalSpecificDependencies=['ZA:Anterior inferior cerebellar artery.l','ZA:Anterior inferior cerebellar artery.r','ZA:Lateral pontine branches of basilar artery.l','ZA:Lateral pontine branches of basilar artery.r','ZA:Medial pontine branches of basilar artery.l','ZA:Medial pontine branches of basilar artery.r','ZA:Posterior inferior cerebellar artery.l','ZA:Posterior inferior cerebellar artery.r','ZA:Anterior spinal artery.l','ZA:Anterior spinal artery.r'],
        priorPortWitnessSha256=digest_file(BASE/'native-terminal-bone-port-witnesses.json'),bothBasilarMouthWitnessSha256=digest_file(BASE/'basilar-both-terminal-bone-witnesses.json'),
        distalBasilarAndGenuineOutsideSupportOriginsProtected=True,skeletonAndInstalledHeartFixed=True,noCandidateGeometryOrPartitionYet=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'contract.json',result);print(json.dumps(dict(contract=str(OUT/'contract.json'),terms=len(terms),nativeSources=len(bindings),courses=courses,elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
