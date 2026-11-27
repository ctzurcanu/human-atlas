"""Rebase only the seven reviewed right-lingual records onto live local state.

No job is run and no installed catalogue is written. New unrelated parts,
left lingual owners, heart and skeleton remain literal in the resulting file.
Any changed/new geometry near this family requires a new context review.
"""
from pathlib import Path
import argparse,copy,json,sys,time,gzip
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_bytes,digest_file
T=time.monotonic();args=argparse.ArgumentParser();args.add_argument('--output',default='male-right-lingual-install-proposal-v2');args=args.parse_args();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';BASE=B/'male-left-lingual-install-proposal-v1/proposed-atlas.json';OLD=B/'male-right-lingual-install-proposal-v1';CURRENT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(BASE)=='978f5a0a25d505f1fcfcca9d2b4cd9fd07023ae22e1363f15131dcd7d0e750b1';assert digest_file(OLD/'proposed-atlas.json')=='736b56fbf4dca7253c54f2c234cc5c0ab4da550f4073e7bda9078b8860ae42b4';sha=digest_file(CURRENT);base=json.loads(BASE.read_bytes());old=json.loads((OLD/'proposed-atlas.json').read_bytes());live=json.loads(CURRENT.read_bytes());job=json.loads((OLD/'catalogue-candidate-job.json').read_text());parameters=job['parameters'];owned=set(parameters['changed_parts']);assert len(owned)==7;bp={p['id']:p for p in base['parts']};op={p['id']:p for p in old['parts']};lp={p['id']:p for p in live['parts']};assert all(lp.get(pid)==bp.get(pid)for pid in owned);assert live['chunks'][:len(base['chunks'])]==base['chunks'];assert len(old['chunks'])==len(base['chunks'])+1
context=json.loads((B/'male-right-lingual-finite-network-v2/nearby-context-emitted-v1/evidence.json').read_text());assert context['passed'];assert all(lp[r['part_id']]==r['record']for r in context['checked']);low,high=np.array(context['expanded_bounds_m']);unrelated_changed={pid for pid in bp.keys()|lp.keys()if bp.get(pid)!=lp.get(pid)};assert not unrelated_changed&owned
for pid in unrelated_changed:
    part=lp.get(pid)
    if part is None or part.get('suppressed'):continue
    bounds=part.get('bounds');assert bounds is not None,('Cannot prove new context is distant',pid)
    assert not(np.all(np.array(bounds)[1]>=low)and np.all(np.array(bounds)[0]<=high)),('Changed/new nearby anatomy requires refreshed context check',pid)
leftjob=json.loads((B/'male-left-lingual-install-proposal-v1/catalogue-candidate-job.json').read_text());left_ids=leftjob['parameters']['changed_parts'];assert all(lp[pid]==bp[pid]for pid in left_ids)
after=copy.deepcopy(live);chunk=len(after['chunks']);bank=copy.deepcopy(old['chunks'][-1]);assert all(c['url']!=bank['url']for c in after['chunks']);after['chunks'].append(bank)
for pid in owned:
    part=copy.deepcopy(op[pid]);part['chunk']=chunk
    if pid in lp:after['parts'][after['parts'].index(lp[pid])]=part
    else:after['parts'].append(part)
bc={c['id']:c for c in base['concepts']};lc={c['id']:c for c in live['concepts']};newconcepts=[c for c in old['concepts']if c['id']not in bc];assert len(newconcepts)==5 and all(c['id']not in lc for c in newconcepts);after['concepts'].extend(copy.deepcopy(newconcepts));ap={p['id']:p for p in after['parts']};assert {pid for pid in lp.keys()|ap.keys()if lp.get(pid)!=ap.get(pid)}==owned;assert all(ap[pid]==part for pid,part in lp.items()if pid not in owned)
out=B/args.output;out.mkdir(exist_ok=False);atomic_bytes(out/'source-atlas.json.gz',gzip.compress(CURRENT.read_bytes(),mtime=0));proposal=out/'proposed-atlas.json';atomic_json(proposal,after);inventory=json.loads(Path(parameters['inventory']).read_text());inventory.update(source_catalogue=str(CURRENT),source_catalogue_sha256=sha);inv=out/'source-bound-inventory.json';atomic_json(inv,inventory)
proof=dict(schema='right-lingual-seven-record-source-preserving-rebase/v1',base_catalogue_sha256=digest_file(BASE),current_source_sha256=sha,proposed_sha256=digest_file(proposal),changed_right_records=sorted(owned),all_live_unrelated_records_literal=True,installed_left_records_literal=True,preserved_unrelated_changed_or_new_records=sorted(unrelated_changed),native_nearby_context_records_literal=len(context['checked']),all_new_context_geometry_outside_family_expanded_bounds=True,source_bank_literal=True,only_new_bank_chunk_index_rebound=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-T);atomic_json(out/'current-source-rebase-proof.json',proof)
parameters.update(expected_catalogue_sha256=sha,proposed_catalogue=str(proposal),inventory=str(inv));parameters['evidence'].update(current_source_rebase=str(out/'current-source-rebase-proof.json'),source_bound_inventory=str(inv));parameters['assembly']['preserved_unrelated_current_revision']=sha
review_path=B/'male-right-lingual-finite-network-v2/actual-native-render-review.json';review=json.loads(review_path.read_text());assert review['accepted_for_local_first_version'];parameters['evidence']['actual-native-render-review']=str(review_path)
for ri,render in enumerate(review['reviews']):
    assert digest_file(Path(render['result']))==render['result_sha256'];parameters['evidence']['native-render-result-'+str(ri)]=render['result']
    for vi,image in enumerate(render['images']):
        assert digest_file(Path(image['path']))==image['sha256'];parameters['evidence'][f'native-render-{ri}-view-{vi}']=image['path']
atomic_json(out/'catalogue-candidate-job.json',job);assert digest_file(CURRENT)==sha;print(json.dumps(dict(source_sha256=sha,proposed=str(proposal),proposed_sha256=digest_file(proposal),job=str(out/'catalogue-candidate-job.json'),preserved_unrelated_records=len(unrelated_changed),elapsed_seconds=time.monotonic()-T)))
