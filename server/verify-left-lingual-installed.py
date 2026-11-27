"""Independent read-only postinstall check; run only after a successful apply.

No installer helpers are used. This compares immutable source/proposal bytes,
every physical record, actual ordinary-body decoding, and the real cardiac
bundle loader without a skip or a browser camera interaction.
"""
from pathlib import Path
import argparse,json,subprocess,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader

SOURCE=ROOT/'.local-models/geometry-workbench/jobs/b9c8595514cdc1b8671263ec425dbf7bbad8f059c484638ae4b24211c0b2d9dc/prepared-2/proposal.json'
SOURCE_SHA='9cf8b30456620f7b7bd6f52311ba818fc94f935beee24e6e07cb953d14106f64'
STAGE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-lingual-install-proposal-v1'
PROPOSED_SHA='978f5a0a25d505f1fcfcca9d2b4cd9fd07023ae22e1363f15131dcd7d0e750b1'
FEMALE_SHA='8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'
BUNDLES={'male':'375bed821531d0bb6c85f5b330362161c9e50dcef2af40f670292844a03e3040','female':'f88dd44e18df5fe343d3a9724a68ef4d76d17803b7f0b388f9eb26d5772bded5'}

def main():
    start=time.monotonic();parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=STAGE/'independent-after-install.json');args=parser.parse_args()
    current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';proposal=STAGE/'proposed-atlas.json'
    assert digest_file(SOURCE)==SOURCE_SHA and digest_file(proposal)==PROPOSED_SHA
    assert digest_file(current)==PROPOSED_SHA,'Apply has not completed, or the installed master differs from the reviewed proposal'
    old=json.loads(SOURCE.read_bytes());expected=json.loads(proposal.read_bytes());actual=json.loads(current.read_bytes());assert actual==expected
    job=json.loads((STAGE/'catalogue-candidate-job.json').read_text())['parameters'];allowed=set(job['changed_parts']);assert len(allowed)==7
    bp={p['id']:p for p in old['parts']};ap={p['id']:p for p in actual['parts']};assert len(bp)==len(old['parts'])and len(ap)==len(actual['parts']);changed={pid for pid in bp.keys()|ap.keys()if bp.get(pid)!=ap.get(pid)};assert changed==allowed
    new=allowed-bp.keys();assert len(new)==6 and allowed&bp.keys()=={'ZA:External carotid artery.l'}
    assert all(ap[pid]==part for pid,part in bp.items()if pid not in allowed)
    skeleton={pid for pid,p in bp.items()if p.get('system')=='skeletal'};assert skeleton and all(ap[pid]==bp[pid]for pid in skeleton)
    heart_proof=json.loads((ROOT/'.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/preview-proof.json').read_text());cardiac={r['part']for r in heart_proof['physical_tissues']};assert len(cardiac)==33
    cardiac.update(pid for key,value in old['finiteCardiacIntegration'].items()if key.endswith('PartIds')for pid in value)
    assert all(ap[pid]==bp[pid]for pid in cardiac)
    assert actual['chunks'][:len(old['chunks'])]==old['chunks']and len(actual['chunks'])==len(old['chunks'])+1
    female=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';assert digest_file(female)==FEMALE_SHA
    for sex,sha in BUNDLES.items():
        path=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/cardiac-label-sites.json';assert digest_file(path)==sha;bundle=json.loads(path.read_text());assert bundle['integration']['labeledTermCount']==113
    bodypath=ROOT/'.local-models/organ-models/spinal-cord/male-body.json';body=json.loads(bodypath.read_bytes());mapping=body['sourceChunkIndices'];assert len(mapping)==len(body['chunks'])and len(set(mapping))==len(mapping);assert body['sourceChunkCount']==len(actual['chunks']);rows={p['id']:p for p in body['parts']};assert new<=rows.keys();reader=AtlasReader();decoded=[]
    for pid in sorted(allowed):
        row=rows[pid];master=dict(row,chunk=mapping[row['chunk']]);assert master==ap[pid],('Ordinary body dropped or changed a new owner',pid)
        ordinary=reader.read_part(body,row,frame='model:local-ta98:assembled-metres',catalogue_path=bodypath);canonical=reader.read_part(actual,ap[pid],frame=ordinary.frame,catalogue_path=current)
        assert ordinary.geometry_sha==canonical.geometry_sha and ordinary.attributes_sha==canonical.attributes_sha
        assert ordinary.source['bank_sha256']==canonical.source['bank_sha256'];decoded.append(dict(part_id=pid,geometry_sha256=ordinary.geometry_sha,bank_sha256=ordinary.source['bank_sha256']))
    terms={f'A12.2.05.{n:03d}'for n in range(15,20)};assert {rows[pid].get('ta98Term')for pid in new}==terms
    for term in terms:
        concepts=[c for c in body['concepts']if c.get('ta98Term')==term and c.get('ta98Kind')=='exact'];assert concepts and any(set(c['elements'])&new for c in concepts),(term,'No ordinary-body direct control')
    # Use the actual application loader, rather than counting bundle records.
    javascript="""
import assert from 'node:assert/strict';import fs from 'node:fs';import {build} from 'esbuild';
const built=await build({entryPoints:['app/cardiac-label-sites.ts'],bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const result={};for(const sex of ['male','female']){
 const raw=JSON.parse(fs.readFileSync(`.local-models/organ-models/spinal-cord/${sex}-body.json`));
 const bundle=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/cardiac-label-sites.json`));
 assert.ok(raw.sourceChunkIndices);const loaded=api.withCardiacLabelBundle(raw,bundle);assert.ok(api.cardiacLabelRegistry(loaded));
 result[sex]=api.cardiacLabeledTermParts(loaded).size;assert.equal(result[sex],113);
}process.stdout.write(JSON.stringify(result));
"""
    node_start=time.monotonic();node=subprocess.run(['node','--input-type=module','--eval',javascript],cwd=ROOT,capture_output=True,text=True,timeout=20,check=True);controls=json.loads(node.stdout);node_elapsed=time.monotonic()-node_start
    report=dict(schema='independent-installed-left-lingual-family/v1',passed=True,source_catalogue_sha256=SOURCE_SHA,installed_catalogue_sha256=PROPOSED_SHA,exact_changed_records=sorted(changed),new_owners=6,canonical_TA98_terms=sorted(terms),unchanged_source_records=len(bp)-1,unchanged_skeleton_records=len(skeleton),unchanged_cardiac_records=len(cardiac),unchanged_female_catalogue_sha256=FEMALE_SHA,unchanged_cardiac_bundle_sha256=BUNDLES,ordinary_body_actual_bank_decoding=decoded,ordinary_body_cardiac_loader_controls=controls,cardiac_loader_skips=0,node_timeout_seconds=20,node_elapsed_seconds=node_elapsed,anatomical_completion=False,elapsed_seconds=time.monotonic()-start)
    atomic_json(args.output,report);print(json.dumps({k:v for k,v in report.items()if k not in ['exact_changed_records','ordinary_body_actual_bank_decoding']},sort_keys=True))

if __name__=='__main__':main()
