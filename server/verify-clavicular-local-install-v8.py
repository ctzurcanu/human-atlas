"""Verify the installed four-record change, served assets and preserved models."""
import hashlib,json,sys,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
W=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracoacromial-clavicular-branch-v1'
D=W/'bilateral-packed-exterior-v8';snap=json.loads((D/'pre-install-snapshot.json').read_text())
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
J=lambda p:hashlib.sha256(json.dumps(p,sort_keys=True,separators=(',',':')).encode()).hexdigest()
master=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';current=json.loads(master.read_text());parts={p['id']:p for p in current['parts']};new=set(snap['new_parts']);changed=set(snap['changed_existing']);old=snap['part_record_hashes']
assert parts.keys()-old.keys()==new and not old.keys()-parts.keys()
untouched=[pid for pid,h in old.items()if pid not in changed];assert all(J(parts[pid])==old[pid]for pid in untouched)
assert H(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json')==snap['female_catalogue_sha256']
assert H(ROOT/'.local-models/organ-models/spinal-cord/female-body.json')==snap['female_body_sha256']
bodypath=ROOT/'.local-models/organ-models/spinal-cord/male-body.json';body=json.loads(bodypath.read_text());served={p['id']:p for p in body['parts']};reader=AtlasReader();rows=[]
for pid in sorted(new|changed):
 assert pid in served and pid not in set(body.get('defaultHidden',[]))
 a=reader.read_part(current,parts[pid],frame='model:local-ta98:assembled-metres',catalogue_path=master)
 b=reader.read_part(body,served[pid],frame='model:local-ta98:assembled-metres',catalogue_path=bodypath)
 assert a.geometry_sha==b.geometry_sha and a.attributes_sha==b.attributes_sha
 url='http://localhost:3016'+body['chunks'][served[pid]['chunk']]['url']
 with urllib.request.urlopen(url,timeout=8)as response:assert response.status==200;prefix=response.read(64)
 assert prefix
 rows.append(dict(part_id=pid,geometry_sha256=a.geometry_sha,attributes_sha256=a.attributes_sha,http_status=200,vertex_count=len(a.vertices),triangles=len(a.faces)))
for pid in new:
 c=next(c for c in body['concepts']if c['id']==pid);assert c['ta98Term']=='A12.2.09.008'and c['ta98Kind']=='exact'
proof=dict(complete=True,passed=True,current_master_sha256=H(master),current_body_sha256=H(bodypath),part_count=len(parts),unchanged_existing_parts=len(untouched),changed_existing_parts=sorted(changed),new_parts=sorted(new),new_TA98_terms=['A12.2.09.008'],skeleton_basilar_pons_and_all_unrelated_records_literal=True,female_catalogue_and_served_body_literal=True,decoded_current_served_geometry_and_attributes_literal=True,http_assets_verified=rows,anatomical_acceptance=False)
(D/'installed-verification.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
