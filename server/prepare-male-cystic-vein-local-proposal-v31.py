"""Prepare only the passed venous family on the explicitly installed revision."""
from pathlib import Path
import sys,json,time,copy,hashlib,sqlite3
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion'
WORK=BASE/'male-cystic-blood-donor-v1';OUT=BASE/'male-cystic-vein-install-proposal-v3'
CURRENT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='5f03a758aeea29e282be5a5dc077a663eeacc1c5e4141acf57b57e630ca3642a'
FRAME='model:local-ta98:assembled-metres'


def main():
    began=time.monotonic();assert digest_file(CURRENT)==EXPECTED;assert not OUT.exists();OUT.mkdir()
    before=json.loads(CURRENT.read_text());after=copy.deepcopy(before);parts={p['id']:p for p in before['parts']};reader=AtlasReader()
    stage=WORK/'vein-recipient-bank-v29/vein';manifest=json.loads((stage/'manifest.json').read_text());assert manifest['emittedGeometryChecksPassed']
    source=json.loads((WORK/'vein-recipient-partition-v28/vein/source-records.json').read_text());records=json.loads((stage/'records.json').read_text());parent_id=manifest['parentId'];branch_id=manifest['branchId'];bank=Path(manifest['bank']);assert digest_file(bank)==manifest['bankSha256']
    for pid,record in source['records'].items():assert parts[pid]==record,'Native hepatic selector changed since geometry checks: '+pid
    host=reader.read_part(before,parts[parent_id],frame=FRAME,catalogue_path=CURRENT);assert host.geometry_sha==source['parent_geometry_sha256'] and host.source['bank_sha256']==source['parent_bank_sha256']
    gb=[]
    for number in range(4,20):
        saved=Mesh.load(WORK/('fixed-target-'+str(number)+'.npz'));pid=saved.source['record']['id'];current=reader.read_part(before,parts[pid],frame=FRAME,catalogue_path=CURRENT);assert current.geometry_sha==saved.geometry_sha and current.attributes_sha==saved.attributes_sha;gb.append(dict(id=pid,geometry_sha256=current.geometry_sha,bank_sha256=current.source['bank_sha256']))
    chunk=len(after['chunks']);url='/local-models/'+str(bank.relative_to(ROOT/'.local-models'));after['chunks'].append(dict(url=url,bytes=bank.stat().st_size))
    parent=records[parent_id];n=parent['vertexCount'];raw=bank.read_bytes();signature=hashlib.sha256(b''.join(raw[parent[k]:parent[k]+n*width]for k,width in [('positions',12),('normals',6),('normalSigns',1)])+json.dumps([parent['geometryOrigin'],1.]).encode()).hexdigest()
    for pid,p in records.items():
        p['chunk']=chunk
        if pid!=branch_id:p['hepaticAttributeSignature']=signature
    branch=records[branch_id];branch['provenance']['detail']='Original 18-control proximal profile and separately transported source sections; actual native right-portal venous blood-space junction. Finite wall thickness is illustrative.';branch['groups']=list(dict.fromkeys(branch.get('groups',[])+['Gallbladder','Gallbladder venous drainage']));branch['ta98Kind']='exact';assert branch_id not in parts
    after['parts']=[copy.deepcopy(records[p['id']])if p['id']in records else p for p in after['parts']];after['parts'].append(copy.deepcopy(branch))
    limitation='Local macroscopic cystic vein, transported from its own HRA source sections without section shrinkage. A connected smooth transverse corridor clears the fixed installed cystic artery and checked biliary/native vascular tissues. The vein communicates with actual native right-portal blood space and enters the same gallbladder fundus wall region with a 0.300 mm terminal-annulus shift, without breaching its innermost surface. First eight portal sections remain literal. Finite wall thickness and microscopic terminal perfusion are illustrative; sampled maximum station turn is 31.14 degrees and requires actual emitted-course review. The fixed pre-existing cystic artery retains 603 crossing face pairs with the whole hepatic-biliary mesh; no new vein/context crossing is detected.'
    concept=dict(id=branch_id,name='Cystic vein',elements=[branch_id],partOf='TA98HEPATIC:male:region:A12.3.12.002',ta98Term='A12.3.12.012',ta98Kind='exact',coverageLimitation=limitation);after['concepts'].append(concept)
    network_concept=next(c for c in after['concepts']if c['id']=='HEPATIC:male:network:portal');assert network_concept['elements']==[parent_id];network_concept['elements'].append(branch_id)
    after['ta98TermBindings'].append(dict(term='A12.3.12.012',conceptId=branch_id,partIds=[branch_id]))
    after['hepaticConnectionIntegration']['banks']['portal']=dict(url=url,sha256=digest_file(bank));after['hepaticConnectionIntegration']['termGroups'].append(dict(term='A12.3.12.012',conceptId=branch_id,parts=[branch_id]))
    # Only this obsolete metadata sentence changes on the protected artery;
    # every physical arterial record and its bank remain literally unchanged.
    obsolete='The separately constructed cystic vein is held outside this candidate because it crosses this artery; the native portal tree remains literal.'
    corrected='The cystic vein is represented separately and joins actual right-portal blood space. The fixed cystic artery retains 603 crossing face pairs with the whole hepatic-biliary mesh; microscopic perfusion is unresolved.'
    artery_concept=next(c for c in after['concepts']if c['id']=='TA98CYSTICBLOOD:male:artery');assert obsolete in artery_concept['coverageLimitation'];artery_concept['coverageLimitation']=artery_concept['coverageLimitation'].replace(obsolete,corrected)
    artery_limits=after['cysticArterialIntegration']['remainingLimits'];assert len(artery_limits)==1 and obsolete in artery_limits[0];after['cysticArterialIntegration']['remainingLimits']=[artery_limits[0].replace(obsolete,corrected)]
    assert all(after['cysticArterialIntegration'][k]==v for k,v in before['cysticArterialIntegration'].items()if k!='remainingLimits')
    after['cysticVenousIntegration']=dict(schema='male-cystic-vein/v1',sourceCatalogueSha256=EXPECTED,previousGeometrySourceCatalogueSha256=source['catalogue_sha256'],parentPartId=parent_id,partIds=[branch_id],regionPartIds=manifest['records'][2:],recipient='ZA:Gallbladder',bankSha256=digest_file(bank),sourceDonor='HRA:VH_F_cystic_vein',sourceLicense='CC BY 4.0',sourceSectionDiameterTolerance=.15,sourceSectionDiameterMaximumRelativeError=.04495321282431752,originalProximalProfileControls=18,remainingLimits=[limitation],anatomicallyAccepted=False)
    changed=sorted(records);ap={p['id']:p for p in after['parts']};assert sorted(pid for pid in parts.keys()|ap.keys()if parts.get(pid)!=ap.get(pid))==changed
    assert all(ap[pid]==p for pid,p in parts.items()if pid not in records)
    after['partsCount']=len(after['parts']);after['triangles']=int(before['triangles']+sum(p.get('triangles',0)for p in after['parts'])-sum(p.get('triangles',0)for p in before['parts']))
    assert all(after.get(k)==before.get(k)for k in set(before)|set(after)if k not in ['parts','concepts','chunks','partsCount','triangles','ta98TermBindings','hepaticConnectionIntegration','cysticVenousIntegration','cysticArterialIntegration'])
    proposal=OUT/'proposed-atlas.json';atomic_json(proposal,after)
    # Layout binds the exact already audited packed buffers, including its
    # appended regional index views, without rewriting the geometry bank.
    layout=OUT/'audited-finite-layout.json';atomic_json(layout,dict(schema='shared-finite-tissues/v1',frame=FRAME,bank=str(bank),bank_sha256=digest_file(bank),origin=parent['geometryOrigin'],parts={pid:records[pid]for pid in [parent_id,branch_id]},passed=True,checks=json.loads((stage/'checks.json').read_text())))
    authority=ROOT/'.local-models/ta98-authority.sqlite';corrections=ROOT/'geometry_tools/ta98_authority.py';con=sqlite3.connect(authority);con.row_factory=sqlite3.Row;term=dict(con.execute('SELECT * FROM ta98 WHERE id=?',('A12.3.12.012',)).fetchone());assert not con.execute('SELECT id FROM ta98 WHERE parent_id=?',(term['id'],)).fetchall()
    inventory=dict(schema='ta98-geometry-inventory/v1',model='local-ta98',model_sex='male',root_term=term['id'],authority_path=str(authority),authority_sha256=digest_file(authority),authority_corrections_path=str(corrections),authority_corrections_sha256=digest_file(corrections),source_catalogue=str(CURRENT),source_catalogue_sha256=EXPECTED,term_count=1,applicable_term_count=1,checked_in_report=str(ROOT/'reports/ta98-cystic-blood-donor-inventory-v1.json'),automatic_mapping_is_not_representation_acceptance=True,anatomical_acceptance=False,terms=[dict(authority=term,applicable_to_model=True,representation_class='tissue-region',source_meshes=[dict(part_id='HRA:VH_F_cystic_vein',source_catalogue=str(ROOT/'public/models/atlas-hra-female.json'),geometry_sha256=Mesh.load(WORK/'donor-1.npz').geometry_sha,role='Original source-specific course and 18-control proximal profile'),dict(part_id=parent_id,geometry_sha256=host.geometry_sha,role='Actual native parent wall and blood-space host; not an alias of the new cystic branch')],missing_distinct_representation=False,distinct_representation_verified=False,planned_representation=dict(part_ids=[branch_id],representation_class='finite-vessel-wall',anatomical_scope=limitation))]);inv=OUT/'source-bound-inventory.json';atomic_json(inv,inventory)
    params={}
    for name,pid,peer in [('first',parent_id,branch_id),('second',branch_id,parent_id)]:
        params[name+'_faces']=[i for g in records[pid]['tissueDrawGroups']if g.get('neighbor')==peer for i in range(g['start']//3,(g['start']+g['count'])//3)];assert len(params[name+'_faces'])==36
    interface=dict(first=parent_id,second=branch_id,kind='shared-interface',parameters=params)
    proof=dict(schema='male-cystic-vein-source-rebase/v21',sourceCatalogueSha256=EXPECTED,oldGeometrySourceCatalogueSha256=source['catalogue_sha256'],proposedCatalogueSha256=digest_file(proposal),changedParts=changed,newPhysicalTA98Terms=[term['id']],allUnrelatedRecordsAndBanksLiteral=True,nativeArterialTreeLiteral=True,installedLingualFamiliesLiteral=True,onlyArterialConceptMetadataChanged=True,obsoleteHeldVeinSentenceRemoved=True,heart33And113LabelIntegrationLiteral=True,affectedNamedRegionIds=manifest['records'][2:],actualNativeParentGeometrySha256=host.geometry_sha,actualNativeParentBankSha256=host.source['bank_sha256'],fixedGallbladderSources=gb,anatomicalAcceptance=False,installed=False,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'source-rebase-proof.json',proof)
    evidence=dict(source_bound_inventory=str(inv),checked_in_authority_inventory=str(ROOT/'reports/ta98-cystic-blood-donor-inventory-v1.json'),source_rebase_proof=str(OUT/'source-rebase-proof.json'),source_actual_parent_sheets=str(WORK/'native-parent-envelopes-v9/proof.json'),continuous_profile_and_patency=str(WORK/'vein-recipient-course-v27/gates.json'),ancestry_partition=str(WORK/'vein-recipient-partition-v28/vein/proof.json'),actual_packed_bank=str(stage/'manifest.json'),emitted_checks=str(stage/'checks.json'),packed_selectors_and_context=str(WORK/'vein-recipient-emitted-v30/proof.json'),branch_self_replay=str(WORK/'vein-recipient-course-v27/branch-self-replay.json'),new_vein_delta_and_inherited_defect=str(WORK/'vein-recipient-emitted-v30/delta-and-inherited-defect.json'),new_parent_support_context=str(WORK/'vein-recipient-emitted-v30/new-parent-support-context.json'),actual_emitted_longitudinal_review=str(WORK/'vein-recipient-course-v27/true-longitudinal-review.json'))
    assembly=dict(whole_organ=dict(name='Native male portal tree and source-specific cystic vein',authority_root=term['id'],part_ids=changed,scope=limitation),adjoining_interfaces=[dict(kind='shared-annular-finite-tissue-interface',parent=parent_id,child=branch_id,triangles=36),dict(kind='actual-native-right-portal-blood-space-mouth',evidence=str(WORK/'vein-recipient-course-v27/gates.json')),dict(kind='gallbladder-wall-recipient',innermostCrossings=0,evidence=str(WORK/'vein-recipient-emitted-v30/proof.json')),dict(kind='native-upstream-annular-interface',sourceNativeTissueFacesLiteral=True)],fixed_skeleton=True,fixed_heart=True,fixed_gallbladder=True,fixed_arterial_tree=True,fixed_all_unrelated_geometry=True,anatomical_acceptance=False)
    spec=dict(model='local-ta98',expected_catalogue_sha256=EXPECTED,proposed_catalogue=str(proposal),changed_parts=changed,closed_parts=[parent_id,branch_id],require_closed=False,inventory=str(inv),assembly=assembly,interface_checks=[interface],evidence=evidence,remaining_limits=[limitation],tissue_audit_checkpoints=str(stage/'validation-checkpoints'),tissue_audit_layouts=[str(layout)],max_pairs=1000000,max_tests=12000000)
    atomic_json(OUT/'catalogue-candidate-job.json',dict(kind='catalogue-candidate',parameters=spec));assert digest_file(CURRENT)==EXPECTED
    print(json.dumps(dict(output=str(OUT),source=EXPECTED,candidate=digest_file(proposal),changedParts=len(changed),newPhysicalTerms=1,nativeArterialLiteral=True,elapsedSeconds=time.monotonic()-began)))


if __name__=='__main__':main()
