"""Add a non-owning canal boundary control without another geometry bank.

Only literal, audited index subranges are referenced. Grey matter retains every
physical face and the void remains empty. The view is hidden except when chosen.
"""
from pathlib import Path
import argparse,copy,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('candidate');ap.add_argument('--boundary',required=True)
args=ap.parse_args();started=time.monotonic();folder=Path(args.candidate).resolve();path=folder/'atlas.json'
boundary=Path(args.boundary).resolve();proof=json.loads(boundary.read_text());a=json.loads(path.read_text())
assert proof['passed']and proof['source_catalogue_sha256']==digest_file(path)
marker=a['finiteSpinalIntegration'];assert marker['bankSha256']==proof['bank_sha256']
assert not marker.get('referencePartIds')and not any(g['term']=='A14.1.02.019'for g in marker['termGroups'])
by_id={p['id']:p for p in a['parts']};cid='TA98CORD:female:term:A14.1.02.019';refs=[];bank_data={}
for number,row in enumerate(proof['face_runs']):
    owner=by_id[row['part_id']];chunk=owner['chunk']
    if chunk not in bank_data:
        url=a['chunks'][chunk]['url'];assert url.startswith('/local-models/')
        bank=ROOT/'.local-models'/url.removeprefix('/local-models/')
        assert digest_file(bank)==proof['bank_sha256'];bank_data[chunk]=bank.read_bytes()
    raw=bank_data[chunk];indices=np.ndarray((row['indexCount']//3,3),'<u4',buffer=raw,offset=row['indices'])
    v=np.ndarray((row['vertexCount'],3),'<f4',buffer=raw,offset=row['positions'])
    used=v[np.unique(indices)].astype(float)+row['geometryOrigin']
    start=(row['indices']-owner['indices'])//4;assert start>=0 and start+row['indexCount']<=owner['indexCount']
    pid=f'TA98CORD:female:canal-boundary:{number}'
    part=dict(id=pid,name='Central canal',conceptId=cid,system='nervous',chunk=chunk,
        **{k:row[k]for k in ['indices','indexCount','positions','normals','normalSigns','vertexCount','geometryOrigin']},
        bounds=[used.min(0).tolist(),used.max(0).tolist()],material=owner['material'],ta98Term='A14.1.02.019',
        sourceId='spinal_cord:central-canal-boundary',defaultHidden=True,pairedTissue=False,sectionSolid=False,
        spinalBoundaryReference=dict(kind='lumen-boundary',ownerPartId=owner['id'],sourceBankSha256=proof['bank_sha256'],
            sourceIndexStart=start,sourceIndexCount=row['indexCount'],physicalOwnerRetained=True),
        coverageLimitation='Actual boundary of the modelled canal void. It is not a solid tube or an ependymal tissue layer.')
    a['parts'].append(part);refs.append(pid)
a['concepts'].append(dict(id=cid,name='Central canal',elements=refs,componentElements=refs,ta98Term='A14.1.02.019',
    ta98Kind='exact',partOf=marker['wholeConceptId'],
    coverageLimitation='Boundary-owned canal void; native terminal ventricle and cranial continuation remain unfinished.'))
marker['referencePartIds']=refs;marker['termGroups'].append(dict(term='A14.1.02.019',conceptId=cid,partIds=refs))
marker['representedCordTerms']+=1;marker['remainingCordTerms']-=1
marker['canalBoundaryEvidence']=dict(bankSha256=proof['bank_sha256'],proofSha256=digest_file(boundary),
    physicalGeometryAdded=False,referenceKind='lumen-boundary',independentSolidAdded=False)
marker['remaining'][0]=str(marker['remainingCordTerms'])+' cord terms remain unrepresented'
a['defaultHidden']=list(dict.fromkeys([*a.get('defaultHidden',[]),*refs]));a['partsCount']=len(a['parts'])
# These are selectable views of existing faces, not new anatomical triangles.
a['spinalPhysicalTriangleCount']=sum(p['indexCount']//3 for p in a['parts']if not p.get('suppressed')and not p.get('spinalBoundaryReference'))
atomic_json(path,a,overwrite=True)
inventory_path=folder/'inventory.json';inventory=json.loads(inventory_path.read_text())
for row in inventory['terms']:
    if row['authority']['id']=='A14.1.02.019':row['candidate_part_ids']=refs
atomic_json(inventory_path,inventory,overwrite=True)
job_path=folder/'job.json';job=json.loads(job_path.read_text());p=job['parameters']
p['require_closed']=False
# The only open pieces are non-owning boundary references. All actual tissue
# cells remain explicitly listed in closed_parts and retain that strict gate.
assert set(p['closed_parts'])==set(marker['partIds']+marker['coveringPartIds'])
p['changed_parts']=sorted(set(p['changed_parts'])|set(refs));p['assembly']['part_ids']+=refs
p['evidence_files'].append(str(boundary));p['assembly']['remaining_limits']=marker['remaining']
atomic_json(job_path,job,overwrite=True)
atomic_json(folder/'canal-control.json',dict(passed=True,reference_parts=len(refs),new_geometry_bytes=0,
    bank_sha256=proof['bank_sha256'],new_finite_tissues=0,canal_is_void=True,anatomicallyAccepted=False,
    candidate_catalogue_sha256=digest_file(path),elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(passed=True,cord_terms=marker['representedCordTerms'],missing=marker['remainingCordTerms'],
    reference_parts=len(refs),new_geometry_bytes=0)),flush=True)
