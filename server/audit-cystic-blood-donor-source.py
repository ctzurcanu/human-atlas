"""Bounded, source-bound cystic blood-vessel donor inspection; never installs."""
from pathlib import Path
import json, sqlite3, time, sys
from collections import Counter, defaultdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, digest_file, atomic_json
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology

WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
SOURCE = ROOT / 'public/models/atlas-hra-female.json'
TARGET = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '9cf8b30456620f7b7bd6f52311ba818fc94f935beee24e6e07cb953d14106f64'
TERMS = {'A12.2.12.031':'HRA:VH_F_cystic_artery', 'A12.3.12.012':'HRA:VH_F_cystic_vein'}

def boundaries(mesh):
    vertices, inv = np.unique(mesh.vertices, axis=0, return_inverse=True)
    faces = inv[mesh.faces]
    edges = Counter(tuple(sorted((int(a),int(b)))) for f in faces for a,b in zip(f,np.roll(f,-1)))
    border = [e for e,n in edges.items() if n == 1]
    adj = defaultdict(set)
    for a,b in border: adj[a].add(b);adj[b].add(a)
    loops=[]
    while adj:
        first=min(adj);sequence=[first];previous=None;current=first
        while True:
            choices=adj[current]-({previous} if previous is not None else set())
            if not choices: break
            nxt=min(choices)
            if nxt==first: break
            sequence.append(nxt);previous,current=current,nxt
            if len(sequence)>len(border): raise ValueError('Invalid source boundary walk')
        points=vertices[sequence]
        loops.append(dict(weldedVertexIds=sequence,coordinates_m=points.tolist(),center_m=points.mean(0).tolist(),closed=first in adj[current],degrees=[len(adj[i]) for i in sequence]))
        for i in sequence: del adj[i]
    return loops

def main():
    began=time.monotonic();assert digest_file(TARGET)==EXPECTED
    assert not WORK.exists();WORK.mkdir(parents=True)
    source=json.loads(SOURCE.read_text());target=json.loads(TARGET.read_text());reader=AtlasReader()
    donor_ids=list(TERMS.values())
    source_ids=donor_ids+['HRA:VH_F_gallbladder','HRA:VH_F_cystic_duct','HRA:VH_F_right_hepatic_artery','HRA:VH_F_hepatic_portal_vein','HRA:VH_F_right_branch_of_portal_vein']
    targets=[p for p in target['parts'] if p['id'].startswith('TA98PART:ZA:Gallbladder:GALLWALL:') or p['id'] in ['TA98HEPATIC:male:artery','TA98HEPATIC:male:portal','TA98HEPATIC:male:region:A12.2.12.030','TA98HEPATIC:male:region:A12.3.12.001']]
    rows=[]
    for label,atlas,path,parts in [('donor',source,SOURCE,[next(p for p in source['parts'] if p['id']==i) for i in source_ids]),('fixed-target',target,TARGET,targets)]:
        for index,part in enumerate(parts):
            mesh=reader.read_part(atlas,part,frame='HRA:united-female-v1.10:original-metres' if label=='donor' else 'model:local-ta98:assembled-metres',catalogue_path=path,catalogue_sha=digest_file(path))
            key=label+'-'+str(index);mesh.save(WORK/(key+'.npz'))
            row=dict(id=part['id'],name=part['name'],role=label,path=str((WORK/(key+'.npz')).relative_to(ROOT)),geometrySha256=mesh.geometry_sha,attributesSha256=mesh.attributes_sha,bounds_m=mesh.bounds.tolist(),sourceBankSha256=mesh.source['bank_sha256'],topology=topology(mesh))
            if part['id'] in donor_ids: row.update(boundaryProfiles=boundaries(mesh),license=part.get('license','CC BY 4.0'),provenance=part.get('provenance'),interpretation='Original named vessel envelope; finite vessel wall and communicating recipient attachment not yet established')
            rows.append(row)
    authority=ROOT/'.local-models/ta98-authority.sqlite'
    with sqlite3.connect('file:'+str(authority)+'?mode=ro',uri=True) as connection:
        connection.row_factory=sqlite3.Row;terms={r['id']:dict(r) for r in connection.execute('SELECT * FROM ta98')}
    descendants=set(TERMS)
    while True:
        more={code for code,row in terms.items() if row['parent_id'] in descendants}
        if more<=descendants:break
        descendants|=more
    inventory=dict(schema='ta98-cystic-blood-donor-inventory/v1',authorityPath=str(authority.relative_to(ROOT)),authoritySha256=digest_file(authority),sourceCatalogue=str(SOURCE.relative_to(ROOT)),sourceCatalogueSha256=digest_file(SOURCE),targetCatalogue=str(TARGET.relative_to(ROOT)),targetCatalogueSha256=EXPECTED,terms=[dict(**terms[code],representationType='3D vessel envelope',sourcePart=TERMS.get(code),currentMaleDirectParts=[],plannedRepresentation='Use original named donor geometry only after topology, source anatomy, positive similarity registration and fixed-context review') for code in sorted(descendants)],parentTerms=[terms[code] for code in sorted({terms[c]['parent_id'] for c in descendants})],allApplicableDescendantsIncluded=True,sourceProvenance=source.get('provenance'),installed=False,anatomicalAcceptance=False)
    atomic_json(ROOT/'reports/ta98-cystic-blood-donor-inventory-v1.json',inventory)
    assembly=dict(schema='cystic-blood-donor-assembly/v1',sourceFrame='HRA:united-female-v1.10:original-metres',targetFrame='model:local-ta98:assembled-metres',movingClosure=donor_ids,donorRegistrationContext=source_ids[2:],fixedTargetParts=[p['id'] for p in targets],skeletonFixed=True,organsFixed=True,interfaces=[dict(owner=donor_ids[0],sourceParent='HRA:VH_F_right_hepatic_artery',targetParent='TA98HEPATIC:male:artery',sourceRecipient='HRA:VH_F_gallbladder',targetRecipient='male gallbladder wall',status='pending actual boundary and wall-entry witnesses'),dict(owner=donor_ids[1],sourceParent='original source portal system; endpoint must be measured',targetParent='portal system or hepatic-bed termination only if source supports it',status='pending original source endpoint classification')],prohibited=['No biliary-lumen entry','No organs or skeleton moved','No lumen plugged by a disk','No invented enlarged venous destination','No installation before native render/cut and context checks'],rows=rows,installed=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(WORK/'source-audit.json',assembly)
    print(json.dumps(dict(work=str(WORK),sourceRows=len(source_ids),fixedRows=len(targets),donors=[{k:r[k] for k in ['id','topology','boundaryProfiles']} for r in rows if r['id'] in donor_ids],elapsedSeconds=assembly['elapsedSeconds'])))

if __name__=='__main__':main()
