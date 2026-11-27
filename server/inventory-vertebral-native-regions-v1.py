"""Full vertebral-subtree inventory and literal source bindings, before partition."""
from pathlib import Path
import json,sqlite3,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file,identity
from geometry_tools.model import AtlasReader
from geometry_tools.contact import self_contact_summary
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    atlas=json.loads(CAT.read_bytes());parts={q['id']:q for q in atlas['parts']};censuspath=ROOT/'reports/ta98-blood-circulation-inventory-installed-v12.json';census=json.loads(censuspath.read_bytes());ct={q['id']:q for q in census['terms']}
    authority=ROOT/'.local-models/ta98-authority.sqlite';db=sqlite3.connect(authority);db.row_factory=sqlite3.Row
    terms=[dict(q)for q in db.execute('WITH RECURSIVE t(id) AS (VALUES(?) UNION ALL SELECT a.id FROM ta98 a JOIN t ON a.parent_id=t.id) SELECT a.* FROM ta98 a JOIN t ON a.id=t.id ORDER BY a.id',('A12.2.08.002',))]
    assert len(terms)==17
    regional={'A12.2.08.003','A12.2.08.004','A12.2.08.009','A12.2.08.010'}
    for row in terms:
        q=ct[row['id']];row.update(existingSourcePartIds=q['active_exact_source_part_ids'],representationStatus=q['representation_status'],
            representationKind='contiguous regional tissue of the one native whole artery'if row['id']in regional else'named macroscopic arterial conduit or branch/network',
            plannedLiteralPartition=row['id']in regional,openingsOrWallLayersOrMicroscopicSolidsNotInferred=True)
    atomic_json(OUT/'authority-inventory.json',dict(schema='vertebral-complete-applicable-ta98-subtree/v1',authoritySha256=digest_file(authority),censusSha256=digest_file(censuspath),
        sourceCatalogueSha256=EXPECTED,terms=terms,total=17,currentSourceTerms=3,missingDirectTerms=14,proposedRegionalTerms=sorted(regional),noGeometryChanged=True))
    priorpath=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1/context-qualified-scan-v2/proof.json';prior=json.loads(priorpath.read_bytes());earlier={q['part_id']:q for q in prior['nativeCourseMeshes']}
    reader=AtlasReader();bindings={};courses=[]
    selected=['ZA:Vertebral artery.l','ZA:Vertebral artery.r','ZA:Basilar artery','ZA:Left subclavian artery','ZA:Right subclavian artery','ZA:Spinal dura']
    selected +=[p['id']for p in atlas['parts']if p['system']=='skeletal'and(p['id'].startswith(('ZA:Atlas (C1):','ZA:Axis (C2):','ZA:Vertebra C3','ZA:Vertebra C4','ZA:Vertebra C5','ZA:Vertebra C6','ZA:Vertebra C7','ZA:Occipital bone:')))]
    selected+=sorted({pid for q in terms for pid in q['existingSourcePartIds']if pid.startswith('ZA:')})
    for pid in dict.fromkeys(selected):
        p=parts[pid];mesh=reader.read_part(atlas,p,frame='model:local-ta98:assembled-metres',catalogue_path=CAT,catalogue_sha=EXPECTED);path=OUT/(identity(pid)[:20]+'.npz');mesh.save(path)
        row=dict(partId=pid,mesh=str(path),sha256=digest_file(path),geometrySha256=mesh.geometry_sha,attributesSha256=mesh.attributes_sha,
            recordSha256=identity(p),sourceBankSha256=mesh.source['bank_sha256'],sourceBankPath=mesh.source['bank_path'],bounds_m=mesh.bounds.tolist(),record=p)
        bindings[pid]=row
        if pid not in earlier or not pid.startswith('ZA:Vertebral artery.'):continue
        old=earlier[pid];assert mesh.geometry_sha==old['geometry_sha256']and mesh.attributes_sha==old['attributes_sha256']and mesh.source['bank_sha256']==old['bank_sha256']
        v,first,inv=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True);f=inv[mesh.faces];rr=recover_station_rings(v,f);assert len(rr)==1
        rings=list(rr[0]['rings'])
        if v[rings[0]].mean(0)[1]>v[rings[-1]].mean(0)[1]:rings=rings[::-1]
        rings[0]=next(r for r in boundary_loops(f)if set(r)==set(rings[0]));rings,phase=align_station_correspondence(v,f,rings);rings=np.array(rings)
        side=pid[-1];sp=OUT/(side+'-literal-source-stations.npz');atomic_npz(sp,physical_vertices=v,physical_faces=f,source_to_physical=inv,rings=rings,centers=v[rings].mean(1),source_faces=mesh.faces)
        selfgate,_=self_contact_summary(mesh,max_pairs=250000,max_tests=2000000);assert selfgate['complete_within_pair_budget']
        courses.append(dict(partId=pid,stationSnapshot=str(sp),stationSnapshotSha256=digest_file(sp),sourceStations=len(rings),sourceWholeTopologyReusedFromProof=old['topology'],
            exactEarlierGeometryAndAttributesAndBankReplayed=True,wholeSourceSelfContact=selfgate,sourceOnlyTwoMouthUnbranched=True,
            sourceInterfaces=['Native subclavian sidewall-origin cast','C6 transverse-foramen entry','C2 transverse-foramen exit','C1 transverse-foramen/groove course','Actual dural entry (not inferred from foramen-magnum plane)','Basilar union cast'],
            regionalBoundariesNotYetChosen=True))
    result=dict(schema='vertebral-native-source-regional-opportunity/v1',sourceCatalogueSha256=EXPECTED,priorSourceProof=str(priorpath),priorSourceProofSha256=digest_file(priorpath),
        inventorySha256=digest_file(OUT/'authority-inventory.json'),sourceBindings=bindings,courses=courses,
        wholeOrganAssembly='Two bilateral one-component outer-source cylinders; region ownership will cover each original source face exactly once, preserving original whole concept as the union.',
        fixedOwners='Every other actual anatomical mesh, including whole skeleton/foramina, installed heart, nerves, dura, skin, subclavian and basilar recipient source anatomy',
        noSourceTransformOrPartition=True,noAtlasChanges=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'source-inventory.json',result);assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(path=str(OUT/'source-inventory.json'),terms=len(terms),sources=len(bindings),elapsedSeconds=result['elapsedSeconds'],courses=[dict(id=q['partId'],stations=q['sourceStations'],selfCrossings=q['wholeSourceSelfContact']['crossing_face_pairs'])for q in courses])))

if __name__=='__main__':main()
