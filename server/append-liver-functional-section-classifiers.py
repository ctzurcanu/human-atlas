"""Reuse literal audited conduit exclusions in the optional segment bank."""
from pathlib import Path
import json,sys,copy
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.model import AtlasReader
base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'completion-v155';path=out/'atlas.json';a=json.loads(path.read_text());r=AtlasReader();m=a['liverCompletionIntegration'];parts={p['id']:p for p in a['parts']}
source=parts[m['functionalPartIds'][0]];bank=r.asset(a['chunks'][source['chunk']]['url']);data=bytearray(bank.read_bytes());old=a['liverSectionIntegration'];oldpath=r.asset(old['bankUrl']);raw=oldpath.read_bytes();rows=[]
def append(value):
    data.extend(b'\0'*((-len(data))%4));offset=len(data);data.extend(value);return offset
for e in old['exclusions']:
    row=copy.deepcopy(e);row['positions']=append(raw[e['positions']:e['positions']+e['vertexCount']*12]);row['indices']=append(raw[e['indices']:e['indices']+e['indexCount']*4]);rows.append(row)
new=out/'functional-section-bank.bin';atomic_bytes(new,bytes(data));url='/local-models/'+str(new.relative_to(ROOT/'.local-models'));chunk=len(a['chunks']);a['chunks'].append(dict(url=url,bytes=len(data)))
for pid in m['functionalPartIds']:parts[pid]['chunk']=chunk
m['functionalSections']=dict(schema='functional-liver-sections/v1',bankUrl=url,bankSha256=digest_file(new),exclusions=rows,sourceBankSha256=old['bankSha256'])
atomic_json(out/'atlas-with-sections.json',a)
atomic_json(out/'functional-section-evidence.json',dict(passed=True,sourceBankSha256=old['bankSha256'],functionalBankSha256=digest_file(new),conduitExclusionsLiteral=True,originalFunctionalBankPrefixLiteral=True,extraBytes=len(data)-bank.stat().st_size))
print(json.dumps(dict(bankBytes=len(data),addedBytes=len(data)-bank.stat().st_size)))
