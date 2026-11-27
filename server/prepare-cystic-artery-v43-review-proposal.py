"""Queue literal emitted reviews and prepare an inert source-bound proposal."""
from pathlib import Path
import sys,time,json,copy,hashlib,sqlite3
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.review import prepare

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion';B=BASE/'male-cystic-blood-donor-v1';C=B/'artery-final-source-frame-v43';OUT=BASE/'male-cystic-artery-correction-proposal-v43';CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345';FRAME='model:local-ta98:assembled-metres'

def crop(mesh,low,high):
    t=mesh.triangles;keep=np.all(t.max(1)>=low,axis=1)&np.all(t.min(1)<=high,axis=1)
    return Mesh(mesh.vertices,mesh.faces[keep],mesh.frame,mesh.attributes,dict(mesh.source,diagnosticOnlyLiteralTriangleSubset=True,retainedActualFaceIds=np.flatnonzero(keep).tolist()))

def queue(parent,branch,current,atlas):
    sec=np.load(C/'candidate/sections.npz');centers=sec['centers'];reader=AtlasReader();vein=reader.read_part(atlas,current['TA98CYSTICBLOOD:male:vein'],frame=FRAME,catalogue_path=CAT);duct=reader.read_part(atlas,current['TA98HEPATIC:male:biliary'],frame=FRAME,catalogue_path=CAT)
    gbid='ZA:Gallbladder';gb=reader.read_part(atlas,current[gbid],frame=FRAME,catalogue_path=CAT)
    out=OUT/'literal-review-inputs';out.mkdir();jobs=[]
    def entry(m,name,color):
        path=out/(name.replace(' ','-')+'.npz');m.save(path);return dict(mesh=str(path),name=name,color=color)
    def submit(name,entries,views):
        request=prepare(None,dict(frame=FRAME,meshes=entries,views=views,resolution=1100,wall_seconds=55));jobs.append(dict(name=name,**request))
    # First review uses the complete branch and only literal surrounding faces
    # in its ROI, so an unrelated distant hepatic branch cannot shrink framing.
    low=branch.triangles.min((0,1))-.004;high=branch.triangles.max((0,1))+.004
    course_entries=[entry(branch,'Full emitted corrected artery',[.88,.12,.08,1]),entry(crop(duct,low,high),'Actual duct around full artery',[.66,.46,.23,.4]),entry(crop(vein,low,high),'Fixed vein around full artery',[.16,.24,.85,1])]
    tangent=centers[18]-centers[12];tangent/=np.linalg.norm(tangent);normal=np.cross(tangent,np.array([1.,0.,0.]));normal/=np.linalg.norm(normal)
    submit('whole-source-course-and-duct-clearance',course_entries,[dict(name='complete-artery-surface-front',direction=[0,0,1]),dict(name='complete-artery-surface-oblique',direction=[1,.25,1]),dict(name='actual-course-longitudinal-cut',direction=normal.tolist(),section=dict(normal=normal.tolist(),point=centers[15].tolist())),dict(name='actual-course-longitudinal-opposite-cut',direction=(-normal).tolist(),section=dict(normal=(-normal).tolist(),point=centers[15].tolist()))])
    # A true local longitudinal plane contains the donor-root direction and
    # native host's local principal direction; its points are actual bank data.
    root=(centers[0]+centers[7])/2;forward=centers[7]-centers[0];forward/=np.linalg.norm(forward);blood=Mesh.load(B/'native-parent-envelopes-v9/artery-blood.npz');distance=np.linalg.norm(blood.vertices-centers[0],axis=1);local=blood.vertices[distance<.004];assert len(local)>20
    hostaxis=np.linalg.svd(local-local.mean(0),full_matrices=False)[2][0];normal=np.cross(forward,hostaxis)
    if np.linalg.norm(normal)<.1:normal=np.cross(forward,np.array([0.,0.,1.]))
    normal/=np.linalg.norm(normal);low=root-.007;high=root+.007
    root_entries=[entry(crop(parent,low,high),'Literal native and proximal arterial root',[.72,.12,.08,1]),entry(crop(branch,low,high),'Literal recipient root of corrected artery',[1,.32,.12,1])]
    submit('actual-native-hepatic-communicating-root',root_entries,[dict(name='native-hepatic-root-longitudinal-1',direction=normal.tolist(),section=dict(normal=normal.tolist(),point=centers[0].tolist())),dict(name='native-hepatic-root-longitudinal-2',direction=(-normal).tolist(),section=dict(normal=(-normal).tolist(),point=centers[0].tolist())),dict(name='shared-annulus-surface-oblique',direction=[1,.25,1])])
    tip=centers[-1];low=tip-.007;high=tip+.007;axis=centers[-1]-centers[-4];axis/=np.linalg.norm(axis);normal=np.cross(axis,np.array([0.,1.,0.]));normal/=np.linalg.norm(normal)
    tip_entries=[entry(crop(branch,low,high),'Literal corrected arterial gallbladder entry',[.88,.12,.08,1]),entry(crop(gb,low,high),'Actual gallbladder entry wall',[.48,.67,.24,.45]),entry(crop(vein,low,high),'Fixed vein at gallbladder entry',[.16,.24,.85,1]),entry(crop(duct,low,high),'Actual nearby biliary wall at entry',[.68,.45,.23,.5])]
    # Some fixed context can be absent from this small ROI; it is omitted,
    # never represented by a synthetic point or substituted geometry.
    tip_entries=[e for e in tip_entries if len(Mesh.load(Path(e['mesh'])).faces)]
    submit('gallbladder-entry-and-biliary-separation',tip_entries,[dict(name='GB-entry-surface-oblique',direction=[1,.25,1]),dict(name='GB-entry-local-longitudinal-cut',direction=normal.tolist(),section=dict(normal=normal.tolist(),point=tip.tolist())),dict(name='GB-entry-opposite-local-longitudinal-cut',direction=(-normal).tolist(),section=dict(normal=(-normal).tolist(),point=tip.tolist()))])
    atomic_json(OUT/'native-review-queue.json',dict(sourceCatalogueSha256=EXPECTED,emittedBankSha256=digest_file(C/'candidate/bank/artery/cells-and-regions.bin'),allReviewCoordinatesAndEffectiveNormalsLiteral=True,diagnosticContextCropsAreTriangleSubsetsOnly=True,diagnosticCutsUncapped=True,jobs=jobs,installed=False));return jobs

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED;assert not OUT.exists();OUT.mkdir();before=json.loads(CAT.read_text());current={p['id']:p for p in before['parts']};after=copy.deepcopy(before);stage=C/'candidate/bank/artery';records=json.loads((stage/'records.json').read_text());manifest=json.loads((stage/'manifest.json').read_text());proof=json.loads((C/'emitted-final-proof.json').read_text());bank=stage/'cells-and-regions.bin';assert digest_file(bank)==proof['bankSha256']==manifest['bankSha256']and proof['wholeEmittedGatesPassed'];assert len(records)==12
    for pid,p in records.items():
        assert pid in current
        for k,v in current[pid].items():
            if k not in {'geometryScale','geometryInstance','sourceOffset','uvs','colors','positions','normals','normalSigns','indices','indexCount','vertexCount','geometryOrigin','triangles','bounds','tissueDrawGroups','pairedTissue','sectionSolid'}:assert p.get(k)==v,(pid,k)
    parentid='TA98HEPATIC:male:artery';branchid='TA98CYSTICBLOOD:male:artery';chunk=len(after['chunks']);url='/local-models/'+str(bank.relative_to(ROOT/'.local-models'));after['chunks'].append(dict(url=url,bytes=bank.stat().st_size));parent=records[parentid];raw=bank.read_bytes();n=parent['vertexCount'];signature=hashlib.sha256(b''.join(raw[parent[k]:parent[k]+n*w]for k,w in [('positions',12),('normals',6),('normalSigns',1)])+json.dumps([parent['geometryOrigin'],1.]).encode()).hexdigest()
    for pid,p in records.items():
        p['chunk']=chunk
        if pid!=branchid:p['hepaticAttributeSignature']=signature
    after['parts']=[copy.deepcopy(records[p['id']])if p['id']in records else p for p in after['parts']];after['partsCount']=len(after['parts']);after['triangles']=int(before['triangles']+sum(p.get('triangles',0)for p in after['parts'])-sum(p.get('triangles',0)for p in before['parts']));after['hepaticConnectionIntegration']['banks']['artery']=dict(url=url,sha256=digest_file(bank));after['cysticArterialIntegration']['bankSha256']=digest_file(bank)
    changed=sorted(records);ap={p['id']:p for p in after['parts']};assert sorted(pid for pid in current if current[pid]!=ap[pid])==changed;assert all(ap[pid]==p for pid,p in current.items()if pid not in records)
    allowed={'parts','chunks','partsCount','triangles','hepaticConnectionIntegration','cysticArterialIntegration'};assert all(after.get(k)==before.get(k)for k in set(after)|set(before)if k not in allowed)
    proposal=OUT/'proposed-atlas.json';atomic_json(proposal,after);atomic_json(OUT/'candidate-records.json',records);reader=AtlasReader();parent_mesh=reader.read_part(after,ap[parentid],frame=FRAME,catalogue_path=proposal);branch_mesh=reader.read_part(after,ap[branchid],frame=FRAME,catalogue_path=proposal);reviews=queue(parent_mesh,branch_mesh,current,before)
    layout=OUT/'audited-finite-layout.json';atomic_json(layout,dict(schema='shared-finite-tissues/v1',frame=FRAME,bank=str(bank),bank_sha256=digest_file(bank),origin=parent['geometryOrigin'],parts={pid:records[pid]for pid in [parentid,branchid]},passed=True,checks=json.loads((stage/'checks.json').read_text())))
    authority=ROOT/'.local-models/ta98-authority.sqlite';correction=ROOT/'geometry_tools/ta98_authority.py';con=sqlite3.connect(authority);con.row_factory=sqlite3.Row;term=dict(con.execute('SELECT * FROM ta98 WHERE id=?',('A12.2.12.031',)).fetchone());assert not con.execute('SELECT id FROM ta98 WHERE parent_id=?',(term['id'],)).fetchall();con.close()
    scope='Correct existing source-specific cystic artery course using proper translations of unchanged paired source profiles. Root rings 0–7, native shared annulus and terminal gallbladder entry remain literal. Actual global hepatic-biliary outer/fluid envelope occupation and whole-wall crossings are zero; fixed cystic vein and other checked native vascular/GB tissues remain clear. Ten regional arterial selectors retain the same oriented source triangles and normal words. Wall thickness and microscopic perfusion remain illustrative. Native emitted review is queued and not yet accepted.'
    oldbranch=reader.read_part(before,current[branchid],frame=FRAME,catalogue_path=CAT)
    inventory=dict(schema='ta98-geometry-inventory/v1',model='local-ta98',model_sex='male',root_term=term['id'],authority_path=str(authority),authority_sha256=digest_file(authority),authority_corrections_path=str(correction),authority_corrections_sha256=digest_file(correction),source_catalogue=str(CAT),source_catalogue_sha256=EXPECTED,term_count=1,applicable_term_count=1,checked_in_report=str(ROOT/'reports/ta98-cystic-blood-donor-inventory-v1.json'),automatic_mapping_is_not_representation_acceptance=True,anatomical_acceptance=False,terms=[dict(authority=term,applicable_to_model=True,representation_class='tissue-region',source_meshes=[dict(part_id=branchid,source_catalogue=str(CAT),geometry_sha256=oldbranch.geometry_sha,role='Existing source-bound paired wall; its source profile dimensions, root and tip are retained')],missing_distinct_representation=False,distinct_representation_verified=False,planned_representation=dict(part_ids=[branchid],representation_class='finite-vessel-wall',anatomical_scope=scope))]);inv=OUT/'source-bound-inventory.json';atomic_json(inv,inventory)
    params={}
    for name,pid,peer in [('first',parentid,branchid),('second',branchid,parentid)]:params[name+'_faces']=[i for g in records[pid]['tissueDrawGroups']if g.get('neighbor')==peer for i in range(g['start']//3,(g['start']+g['count'])//3)];assert len(params[name+'_faces'])==36
    stale=[]
    for c in before.get('concepts',[]):
        if '603'in c.get('coverageLimitation',''):stale.append(dict(kind='concept',id=c['id'],coverageLimitation=c['coverageLimitation']))
    for k in ['cysticArterialIntegration','cysticVenousIntegration']:
        for text in before.get(k,{}).get('remainingLimits',[]):
            if '603'in text:stale.append(dict(kind='integration',id=k,text=text))
    rebase=dict(sourceCatalogueSha256=EXPECTED,sourceRegistryRevision='691a562b61ef319d9e72a14f0f0d6c8211d45222e2b1088bc3656ac6b6582d37',proposedCatalogueSha256=digest_file(proposal),changedParts=changed,physicalOwners=2,regionalControls=10,newTA98Terms=0,allUnrelatedPhysicalRecordsAndBanksLiteral=True,allConceptsAndTermBindingsLiteral=True,allUnrelatedTopLevelMetadataLiteral=True,allowedTopLevelBankReferencesUpdated=['hepaticConnectionIntegration.banks.artery','cysticArterialIntegration.bankSha256'],protectedNativeParentAndRegionalTrianglesAndCornerWordsLiteral=True,protectedCysticBranchFacesAndWordsLiteral=proof['protectedBranchFacesAndCornerWordsLiteral'],staleInheritedDefectMetadataPreservedForRootReview=stale,onlyProspectiveFilesWritten=True,noCachePrepared=True,noGeometryJobCreated=True,noLiveInstall=True,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'source-rebase-proof.json',rebase)
    evidence=dict(source_bound_inventory=str(inv),checked_in_authority_inventory=str(ROOT/'reports/ta98-cystic-blood-donor-inventory-v1.json'),source_rebase_proof=str(OUT/'source-rebase-proof.json'),actual_packed_bank=str(stage/'manifest.json'),emitted_checks=str(stage/'checks.json'),continuous_profile_and_patency=str(C/'candidate/gates.json'),ancestry_partition=str(C/'candidate/owned/artery/proof.json'),complete_emitted_context_and_occupancy=str(C/'emitted-final-proof.json'),actual_global_duct_oracle=str(B/'biliary-actual-envelope-v39/proof.json'),fixed_vein_complete_occupancy=str(B/'cystic-complete-duct-occupancy-v40/proof.json'),native_review_queue=str(OUT/'native-review-queue.json'))
    assembly=dict(whole_organ=dict(name='Corrected native male hepatic arterial and cystic wall family',authority_root=term['id'],part_ids=changed,scope=scope),adjoining_interfaces=[dict(kind='shared-annular-finite-tissue-interface',parent=parentid,child=branchid,triangles=36),dict(kind='literal-native-right-hepatic-blood-space-root',sourceNativeTissueFacesLiteral=True),dict(kind='literal-terminal-gallbladder-wall-entry',innermostCrossings=0),dict(kind='fixed-duct-exterior-exclusion',fullGlobalOccupation_mm3=0)],fixed_skeleton=True,fixed_heart=True,fixed_gallbladder=True,fixed_portal_tree=True,fixed_cystic_vein=True,fixed_all_unrelated_geometry=True,anatomical_acceptance=False)
    spec=dict(model='local-ta98',expected_catalogue_sha256=EXPECTED,proposed_catalogue=str(proposal),changed_parts=changed,closed_parts=[parentid,branchid],require_closed=False,inventory=str(inv),assembly=assembly,interface_checks=[dict(first=parentid,second=branchid,kind='shared-interface',parameters=params)],evidence=evidence,remaining_limits=[scope],tissue_audit_checkpoints=str(stage/'validation-checkpoints'),tissue_audit_layouts=[str(layout)],max_pairs=1000000,max_tests=12000000)
    atomic_json(OUT/'install-spec.json',spec);atomic_json(OUT/'catalogue-candidate-job.json',dict(kind='catalogue-candidate',parameters=spec));assert digest_file(CAT)==EXPECTED
    atomic_json(OUT/'handoff.json',dict(sourceCatalogueSha256=EXPECTED,proposalSha256=digest_file(proposal),bankSha256=digest_file(bank),bank=str(bank),records=str(OUT/'candidate-records.json'),manifest=str(stage/'manifest.json'),checkpoints=str(stage/'validation-checkpoints'),layout=str(layout),installSpec=str(OUT/'install-spec.json'),jobSpec=str(OUT/'catalogue-candidate-job.json'),reviewJobs=[dict(name=r['name'],request_id=r['request_id'],folder=r['folder'])for r in reviews],finalAppCodeAndCacheBindingPendingRoot=True,installed=False,elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(output=str(OUT),proposalSha256=digest_file(proposal),bankSha256=digest_file(bank),changedParts=12,reviewJobs=[dict(name=r['name'],request_id=r['request_id'])for r in reviews],elapsedSeconds=time.monotonic()-began)))
if __name__=='__main__':main()
