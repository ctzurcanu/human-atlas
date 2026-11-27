"""Merge two independently replayed, same-source branch candidates once.

This changes no geometry/banks, current catalogue or database. Right bank index
is remapped after the left bank; every other old record remains literal.
"""
import argparse,copy,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def run(left,right,out):
    start=time.monotonic();master=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';source_sha=H(master);base=read(master);old={p['id']:p for p in base['parts']};n=len(base['chunks']);out.mkdir(exist_ok=False)
    candidates=[];proofs=[];changed=[];newconcepts=[];newbindings=[];banks=[];tri_delta=0
    for side,directory in [('left',left),('right',right)]:
        replay=read(directory/'independent-compiled-replay.json');pack=read(directory/'packing-proof.json');proposal=directory/'proposed-atlas.json';candidate=read(proposal)
        assert replay['passed'] and replay['complete'] and pack['passed'] and pack['complete']
        assert replay['source_catalogue_sha256']==pack['source_catalogue_sha256']==source_sha
        assert H(proposal)==replay['proposed_catalogue_sha256']==pack['proposed_catalogue_sha256']
        suffix='l'if side=='left'else'r';parent=f'ZA:Thoraco-acromial artery.{suffix}';child=f'TA98PART:ZA:Thoraco-acromial artery.{suffix}:A12.2.09.008';parts={p['id']:p for p in candidate['parts']}
        assert parts.keys()-old.keys()=={child} and not old.keys()-parts.keys()
        assert all(parts[k]==v for k,v in old.items()if k!=parent)
        assert candidate['chunks'][:-1]==base['chunks'] and len(candidate['chunks'])==n+1
        assert candidate['materials']==base['materials'] and candidate['concepts'][:-1]==base['concepts'] and candidate['ta98TermBindings'][:-1]==base['ta98TermBindings']
        assert candidate['concepts'][-1]['id']==child and candidate['concepts'][-1]['ta98Term']=='A12.2.09.008'
        changed.extend([parent,child]);newconcepts.append(candidate['concepts'][-1]);newbindings.append(candidate['ta98TermBindings'][-1]);banks.append(candidate['chunks'][-1]);tri_delta+=candidate['triangles']-base['triangles']
        copies={pid:copy.deepcopy(parts[pid])for pid in [parent,child]}
        for rec in copies.values():assert rec['chunk']==n;rec['chunk']=n+len(candidates)
        candidates.append(copies);proofs.append(dict(side=side,directory=str(directory),proposed_catalogue_sha256=H(proposal),independent_proof=str(directory/'independent-compiled-replay.json'),independent_proof_sha256=H(directory/'independent-compiled-replay.json'),bank_sha256=pack['bank_sha256'],bank_bytes=pack['bank_bytes']))
    after=copy.deepcopy(base)
    for records in candidates:
        for pid,rec in records.items():
            if pid in old:after['parts'][next(i for i,x in enumerate(after['parts'])if x['id']==pid)]=rec
            else:after['parts'].append(rec)
    after['chunks'].extend(banks);after['concepts'].extend(newconcepts);after['ta98TermBindings'].extend(newbindings);after['partsCount']=len(after['parts']);after['triangles']=base['triangles']+tri_delta
    assert len(after['parts'])==len(base['parts'])+2 and len(set(changed))==4
    proposal=out/'proposed-atlas.json';dump(proposal,after)
    inventory=left/'registry-authority-inventory.json';iv=read(inventory);assert iv['schema']=='ta98-geometry-inventory/v1' and len(iv['terms'])==7
    (out/'registry-authority-inventory.json').write_bytes(inventory.read_bytes());inventory=out/'registry-authority-inventory.json'
    assembly=dict(schema='bilateral-independent-clavicular-open-ostia/v1',whole_organ='Both complete native thoraco-acromial exteriors, each with one independently constructed illustrative clavicular main/muscular sidecourse and an open parent ostium.',adjoining_interfaces=['Both native parent proximal/distal profiles remain literal and open.','Each side has one shared parent/child exterior ostium with literal packed seam coordinates/normals.','Main and muscular distal courses remain open continuations; no lumen disk, finite wall or complete perfusion claim.','Left and right recipient/bone/vein/nerve context is reviewed independently; skeleton and all unrelated geometry fixed.'],source_catalogue_sha256=source_sha,source_candidates=proofs,sourceOuterSurfaceOnly=True,anatomical_acceptance=False,installed=False)
    dump(out/'assembly.json',assembly);evidence=[str(out/'assembly.json')]
    for d in [left,right]:evidence.extend(str(d/p)for p in ['packing-proof.json','native-tree-proof.json','independent-compiled-replay.json','compiled-merged-audit.json','diagnostic-connected-blood-classifier-audit.json'])
    job=dict(kind='catalogue-candidate',parameters=dict(model='local-ta98',proposed_catalogue=str(proposal),proposed_catalogue_sha256=H(proposal),expected_catalogue_sha256=source_sha,changed_parts=changed,require_closed=False,closed_parts=[],inventory=str(inventory),inventory_sha256=H(inventory),assembly=assembly,evidence=evidence,max_pairs=600000,max_tests=2000000));dump(out/'catalogue-job.json',job)
    proof=dict(schema='bilateral-same-source-clavicular-merge/v1',passed=True,complete=True,source_catalogue_sha256=source_sha,proposed_catalogue_sha256=H(proposal),changed_parts=changed,new_parts=[newconcepts[0]['id'],newconcepts[1]['id']],new_TA98_terms=['A12.2.09.008'],unchanged_existing_parts=len(old)-2,original_materials_and_other_records_literal=True,existing_chunks_literal=True,geometry_or_bank_rewrite=False,banks=proofs,bank_bytes=sum(x['bank_bytes']for x in proofs),right_bank_index=n+1,sourceOuterSurfaceOnly=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start);dump(out/'packing-proof.json',proof);print(json.dumps(proof))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--left',type=Path,required=True);p.add_argument('--right',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.left.resolve(),a.right.resolve(),a.output.resolve())
