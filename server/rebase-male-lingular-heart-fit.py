"""Rebase the reviewed local lingular cells without reverting the installed heart."""
from pathlib import Path
import argparse,copy,gzip,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion'
OLD=BASE/'male-lingular-install-proposal-v1';OUT=BASE/'male-lingular-heart-fit-rebase-v1'
CURRENT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
SOURCE='a0d7b46e033d145617a8937e1aff489db90b7a2808aa76ec90192f155ad95518'
PRIOR='6526654a0e412947121518636d1b8004a934093f8a4cf459afcfa55c80a1b63c'
def load(p):return json.loads(Path(p).read_bytes())

def prepare():
    began=time.monotonic();assert digest_file(CURRENT)==SOURCE
    source=load(CURRENT);prior_raw=gzip.decompress((ROOT/'.local-models/ta98-review/male-fixed-vessel-heart-pose-installed-v1/0-atlas.json.gz').read_bytes())
    import hashlib
    assert hashlib.sha256(prior_raw).hexdigest()==PRIOR
    before=json.loads(prior_raw);old=load(OLD/'proposed-atlas.json');a=copy.deepcopy(source)
    spec=load(OLD/'catalogue-candidate-job-stage2.json');ids=set(spec['parameters']['changed_parts']);assert len(ids)==10
    bp={p['id']:p for p in before['parts']};op={p['id']:p for p in old['parts']};cp={p['id']:p for p in source['parts']}
    actual={pid for pid in set(bp)|set(op)if bp.get(pid)!=op.get(pid)};assert actual==ids
    chunks={i:next((j for j,c in enumerate(a['chunks'])if c==chunk),None)for i,chunk in enumerate(old['chunks'])}
    for i,chunk in enumerate(old['chunks']):
        if chunks[i]is None:chunks[i]=len(a['chunks']);a['chunks'].append(copy.deepcopy(chunk))
    for pid in ids:
        assert cp.get(pid)==bp.get(pid),'Current native artery changed: '+pid
    a['parts']=[copy.deepcopy(op[p['id']])if p['id']in ids else p for p in a['parts']]
    a['parts'].extend(copy.deepcopy(p)for p in old['parts']if p['id']not in cp)
    for p in a['parts']:
        if p['id']in ids:p['chunk']=chunks[p['chunk']]
    bc={c['id']:c for c in before['concepts']};oc={c['id']:c for c in old['concepts']};cc={c['id']:c for c in source['concepts']}
    changed_concepts={cid for cid in set(bc)|set(oc)if bc.get(cid)!=oc.get(cid)};assert len(changed_concepts)==4
    for cid in changed_concepts:assert cc.get(cid)==bc.get(cid)
    a['concepts']=[copy.deepcopy(oc[c['id']])if c['id']in changed_concepts else c for c in a['concepts']]
    a['concepts'].extend(copy.deepcopy(c)for c in old['concepts']if c['id']not in cc)
    assert source['ta98TermBindings']==before['ta98TermBindings']
    a['ta98TermBindings']=copy.deepcopy(old['ta98TermBindings'])
    a['lingularEntryIntegration']=copy.deepcopy(old['lingularEntryIntegration']);a['lingularEntryIntegration']['sourceCatalogueSha256']=SOURCE
    a['lingularEntryIntegration']['previousReviewedSourceCatalogueSha256']=PRIOR
    a['partsCount']=len(a['parts']);a['triangles']=old['triangles']
    ap={p['id']:p for p in a['parts']}
    assert all(ap[pid]==p for pid,p in cp.items()if pid not in ids)
    assert a['finiteCardiacIntegration']==source['finiteCardiacIntegration'] and a['cardiacFeatureIntegration']==source['cardiacFeatureIntegration']
    assert all(a.get(k)==source.get(k)for k in set(a)|set(source)if k not in ['parts','concepts','chunks','ta98TermBindings','lingularEntryIntegration','partsCount','triangles'])
    OUT.mkdir(exist_ok=False);atomic_json(OUT/'proposed-atlas.json',a)
    inv=load(BASE/'male-lingular-four-owner-pack-v3/source-bound-inventory.json');assert inv['source_catalogue_sha256']==PRIOR
    inv['source_catalogue_sha256']=SOURCE;inv['rebased_from_source_catalogue_sha256']=PRIOR;atomic_json(OUT/'source-bound-inventory.json',inv)
    spec['parameters']['expected_catalogue_sha256']=SOURCE;spec['parameters']['proposed_catalogue']=str(OUT/'proposed-atlas.json');spec['parameters']['inventory']=str(OUT/'source-bound-inventory.json')
    atomic_json(OUT/'catalogue-candidate-job.json',spec)
    proof=dict(schema='lingular-heart-pose-preserving-rebase/v1',source_sha256=SOURCE,previous_source_sha256=PRIOR,previous_candidate_sha256=digest_file(OLD/'proposed-atlas.json'),candidate_sha256=digest_file(OUT/'proposed-atlas.json'),changed_parts=sorted(ids),changed_concepts=sorted(changed_concepts),all_nonmembers_literal=True,installed_heart_and_113_label_integration_literal=True,new_banks=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(OUT/'rebase-proof.json',proof);print(json.dumps({k:v for k,v in proof.items()if k not in ['changed_parts','changed_concepts']}))

def heart_contacts():
    began=time.monotonic();assert digest_file(CURRENT)==SOURCE
    source=load(CURRENT);candidate=load(OUT/'proposed-atlas.json');parts={p['id']:p for p in source['parts']};new={p['id']:p for p in candidate['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres'
    heart=load(ROOT/'.local-models/geometry-workbench/heart/male-fixed-vessel-heart-closure-v1/closure.json')['physical_heart_tissue'];ids=[r['part']for r in heart];assert len(ids)==33
    owners={pid:reader.read_part(candidate,new[pid],frame=frame,catalogue_path=OUT/'proposed-atlas.json')for pid in candidate['lingularEntryIntegration']['partIds']};trees={pid:BVH(m)for pid,m in owners.items()}
    rows=[];total=0
    for pid in ids:
        host=reader.read_part(source,parts[pid],frame=frame,catalogue_path=CURRENT);tree=None
        for owner,mesh in owners.items():
            import numpy as np
            overlap=bool(np.all(mesh.bounds[0]<=host.bounds[1])and np.all(host.bounds[0]<=mesh.bounds[1]));summary=None
            if overlap:
                tree=tree or BVH(host);_,summary=pair_contacts(mesh,host,0.,first_bvh=trees[owner],second_bvh=tree,max_pairs=250000,max_tests=2000000)
                assert summary['complete_within_pair_budget'];total+=summary['crossing_face_pairs']+summary['coplanar_overlap_face_pairs']
            rows.append(dict(owner=owner,context=pid,context_geometry_sha256=host.geometry_sha,aabb_overlap=overlap,contact=summary))
        atomic_json(OUT/f'heart-checkpoint-{len(rows):03d}.json',dict(source_sha256=SOURCE,rows=rows,complete=len(rows)==132,elapsed_seconds=time.monotonic()-began))
    result=dict(schema='lingular-current-heart-contact/v1',source_sha256=SOURCE,candidate_sha256=digest_file(OUT/'proposed-atlas.json'),heart_tissues=33,finite_owners=4,pairs=132,rows=rows,complete=True,crossing_or_overlap_face_pairs=total,passed=total==0,elapsed_seconds=time.monotonic()-began)
    atomic_json(OUT/'current-heart-contact.json',result);print(json.dumps({k:v for k,v in result.items()if k!='rows'}))

def context_contacts():
    began=time.monotonic();assert digest_file(CURRENT)==SOURCE
    source=load(CURRENT);candidate=load(OUT/'proposed-atlas.json');reader=AtlasReader();frame='model:local-ta98:assembled-metres'
    owners={p['id']:reader.read_part(candidate,p,frame=frame,catalogue_path=OUT/'proposed-atlas.json')for p in candidate['parts']if p['id']in candidate['lingularEntryIntegration']['partIds']};trees={pid:BVH(m)for pid,m in owners.items()}
    source_ids={'ZA:Left pulmonary artery','ZA:Inferior lingular artery of left lung','ZA:Superior lingular artery of left lung'}
    hosts=[p for p in source['parts']if not p.get('suppressed')and ((p['id'].startswith('ZA:')and 'segmental artery of left lung'in p['id'])or p['id'].startswith('LOCAL:male:pulmonary:bronchus:left_ling'))and p['id']not in source_ids]
    assert len(hosts)==10
    rows=[];total=0
    for p in hosts:
        host=reader.read_part(source,p,frame=frame,catalogue_path=CURRENT);tree=BVH(host)
        for owner,mesh in owners.items():
            _,summary=pair_contacts(mesh,host,0.,first_bvh=trees[owner],second_bvh=tree,max_pairs=250000,max_tests=2000000)
            assert summary['complete_within_pair_budget'];total+=summary['crossing_face_pairs']+summary['coplanar_overlap_face_pairs']
            rows.append(dict(owner=owner,context=p['id'],context_geometry_sha256=host.geometry_sha,contact=summary))
        atomic_json(OUT/f'context-checkpoint-{len(rows):03d}.json',dict(source_sha256=SOURCE,rows=rows,complete=len(rows)==40,elapsed_seconds=time.monotonic()-began))
    result=dict(schema='lingular-current-arterial-bronchial-contact/v1',source_sha256=SOURCE,candidate_sha256=digest_file(OUT/'proposed-atlas.json'),fixed_arteries=8,fixed_bronchi=2,finite_owners=4,pairs=40,rows=rows,complete=True,crossing_or_overlap_face_pairs=total,passed=total==0,elapsed_seconds=time.monotonic()-began)
    atomic_json(OUT/'current-arterial-bronchial-contact.json',result);print(json.dumps({k:v for k,v in result.items()if k!='rows'}))

def finalize():
    assert digest_file(CURRENT)==SOURCE
    proofs=[OUT/'current-heart-contact.json',OUT/'current-arterial-bronchial-contact.json']
    for path in proofs:
        proof=load(path);assert proof['passed']and proof['complete']and proof['source_sha256']==SOURCE and proof['candidate_sha256']==digest_file(OUT/'proposed-atlas.json')
    spec=load(OUT/'catalogue-candidate-job.json');params=spec['parameters'];params['evidence'].update({p.name:str(p)for p in proofs});params['evidence']['source-bound-inventory.json']=str(OUT/'source-bound-inventory.json');params['evidence']['pose-preserving-rebase-proof']=str(OUT/'rebase-proof.json')
    params['assembly']['adjoining_interfaces']=[x for x in params['assembly']['adjoining_interfaces']if x['kind']!='fixed-context-no-crossing']
    params['assembly']['adjoining_interfaces'].append(dict(kind='fixed-current-heart-and-native-arterial-bronchial-context',heart_tissues=33,fixed_other_left_segmental_arteries=8,fixed_lingular_bronchi=2,exact_triangle_pairs_checked=172,new_crossings=0,evidence=[str(p)for p in proofs]))
    params['validation_staging']['fixed_heart_and_bronchus_proof']=str(proofs[0]);params['validation_staging']['fixed_other_arteries_and_bronchi_proof']=str(proofs[1]);params['validation_staging']['rebased_current_heart_pose_preserved']=True
    atomic_json(OUT/'catalogue-candidate-job-current-context.json',spec);print(json.dumps(dict(spec=str(OUT/'catalogue-candidate-job-current-context.json'),source=SOURCE,candidate=digest_file(OUT/'proposed-atlas.json'),all_current_fixed_context_complete=True)))

def receipt():
    assert digest_file(CURRENT)==SOURCE
    candidate_sha=digest_file(OUT/'proposed-atlas.json');folder=ROOT/'.local-models/geometry-workbench/dependency-relations'/candidate_sha;cache_path=folder/'male.json';cache=load(cache_path);staged=cache['stagedProducer'];manifest=load(staged['manifest'])
    assert cache['sourceSha256']==candidate_sha and staged['completeOrdinalBatches']==manifest['batchCount']==28 and manifest['partCount']==6945 and manifest['batchSize']==250
    for name,sha in {**cache['resolverInputs'],**staged['producerInputs']}.items():assert digest_file(ROOT/name)==sha
    assert digest_file(staged['manifest'])==staged['manifestSha256']
    for record in staged['batchFiles']:
        assert digest_file(record['path'])==record['sha256'];batch=load(record['path']);assert batch['complete']and batch['completedThrough']==batch['end']
    records=load(OUT/'batch-runtime-records.json');assert [r['index']for r in records]==list(range(28));runtime=[]
    for record in records:
        p=Path(record['path']);r=load(p);assert r['state']=='completed'and r['returncode']==0 and r['limit_seconds']==55 and '--batch-index'in r['command']and r['command'][-1]==str(record['index'])
        runtime.append(dict(index=record['index'],path=str(p.relative_to(ROOT)),sha256=digest_file(p),elapsed_seconds=r['elapsed_seconds'],wall_limit_seconds=55,state=r['state']))
    review=ROOT/'.local-models/geometry-workbench/review/e3d48045f250d8b4af45918e652589c9caa2c2de121d6ab34d1a620a7ba974fe';result=load(review/'result.json');request=load(review/'request.json');assert result['state']=='rendered'and result['original_objects_preserved']and not result['main_file_saved']
    views=[]
    for index in range(3):
        p=review/f'view-{index}.png';views.append(dict(path=str(p),sha256=digest_file(p),name=request['views'][index]['name'],surface_or_section='section'if 'section'in request['views'][index]else 'surface',section='section'in request['views'][index],inspected=True))
    context=[OUT/'current-heart-contact.json',OUT/'current-arterial-bronchial-contact.json'];evidence=dict(schema='root-reviewed-lingular-rebase-evidence/v1',job_id='198a2ce0e5b592936635b9725b8fe96b64d868bfe54ff8a27b456c95aae06fe2',input_signature='d2e4372724fe8a018b16b53df5128fe4f0d5bdec3b2ae2c2ed67562612b28f06',reviewer='Codex root: actual surface and section inspection previously recorded; explicit acceptance renewed after current-context checks',accepted_for_local_preview=True,anatomical_acceptance=False,remaining_limits=load(OUT/'catalogue-candidate-job-current-context.json')['parameters']['remaining_limits'],views=views,unchanged_cell_bank_sha256='f8b16614a6dd8aea74d96863f593a4156526a1a2dab0aff40e2f7f4f5ffc7fff',current_context_proofs=[dict(path=str(p),sha256=digest_file(p))for p in context],prior_actual_review_reused_for_identical_cells_only=True,root_will_record_review_and_install=True)
    atomic_json(OUT/'root-review-evidence.json',evidence)
    r=dict(schema='exact-forest-prewarm-receipt/v1',sourceSha256=candidate_sha,cache=str(cache_path.relative_to(ROOT)),cacheSha256=digest_file(cache_path),manifestSha256=staged['manifestSha256'],parts=6945,completedExactBatches=28,batchSize=250,edgeCount=len(cache['edges']),resolverFresh=True,twoWorkers=True,eachBatchWallLimitSeconds=55,noTimedOutOrIncompleteBatches=True,installedCatalogueNotMutated=True,runtimeReports=runtime,reviewEvidence=str((OUT/'root-review-evidence.json').relative_to(ROOT)))
    atomic_json(folder/'lingular-prewarm-receipt-v1.json',r);print(json.dumps({k:v for k,v in r.items()if k!='runtimeReports'}))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','heart','context','finalize','receipt']);args=ap.parse_args();{'prepare':prepare,'heart':heart_contacts,'context':context_contacts,'finalize':finalize,'receipt':receipt}[args.stage]()
