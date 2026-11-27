"""Prepare source-verified declaration refinements; never edits the live DB."""
from pathlib import Path
import sys,time,json,copy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader
from geometry_tools.constraints import check_constraint
from geometry_tools.interface_transfer import _sheet

P=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-artery-correction-proposal-v43';OUT=P/'verified-interface-parameters-v1';CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345';PROPOSED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69';FRAME='model:local-ta98:assembled-metres'
PREFIXES=['099ecf5d','88f3253e','ba6d6026']

def drawfaces(record,peer):
    return [i for g in record.get('tissueDrawGroups',[])if g.get('neighbor')==peer for i in range(g['start']//3,(g['start']+g['count'])//3)]

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED and digest_file(P/'proposed-atlas.json')==PROPOSED;assert not OUT.exists();OUT.mkdir()
    before=json.loads(CAT.read_text());after=json.loads((P/'proposed-atlas.json').read_text());bp={p['id']:p for p in before['parts']};ap={p['id']:p for p in after['parts']};declarations=json.loads((P/'constraint-declarations-v1.json').read_text())['constraints'];rows=[r for r in declarations if any(r['id'].startswith(s)for s in PREFIXES)];assert len(rows)==3;reader=AtlasReader();proofs=[];refinements=[];bindings=[];reuse=[]
    for declaration in rows:
        interface=declaration['evidence']['forest_interface'];first,second=interface['from'],interface['to'];assert declaration['kind']=='shared-interface'and declaration['parameters']=={'co_move':False,'verification_pending':True}
        old=[reader.read_part(before,bp[pid],frame=FRAME,catalogue_path=CAT)for pid in [first,second]];new=[reader.read_part(after,ap[pid],frame=FRAME,catalogue_path=P/'proposed-atlas.json')for pid in [first,second]]
        oldparams={'co_move':False};newparams={};aliases=[]
        for role,peer,pid in [('first',second,first),('second',first,second)]:
            physicalpeer=bp[peer].get('hepaticRegionOf',peer);oldids=drawfaces(bp[pid],peer);newids=drawfaces(ap[pid],peer)
            if not oldids:
                assert physicalpeer!=peer and ap[peer].get('hepaticRegionOf')==physicalpeer,'Only an explicit source-bound hepatic regional view may resolve a physical neighbor'
                oldids=drawfaces(bp[pid],physicalpeer);newids=drawfaces(ap[pid],physicalpeer);assert oldids and newids
                # A regional alias is permitted only because the full annulus
                # and every stored normal word are exactly the same owner data.
                alias_index=0 if peer==first else 1;ownerold=reader.read_part(before,bp[physicalpeer],frame=FRAME,catalogue_path=CAT);ownernew=reader.read_part(after,ap[physicalpeer],frame=FRAME,catalogue_path=P/'proposed-atlas.json');owner_old_ids=drawfaces(bp[physicalpeer],pid);owner_new_ids=drawfaces(ap[physicalpeer],pid);region_old_ids=drawfaces(bp[peer],pid);region_new_ids=drawfaces(ap[peer],pid)
                assert _sheet(ownerold,owner_old_ids)==_sheet(old[alias_index],region_old_ids)and _sheet(ownernew,owner_new_ids)==_sheet(new[alias_index],region_new_ids)
                aliases.append(dict(endpoint=peer,declaredPhysicalNeighbor=physicalpeer,sourceRegionalAnnulusAndNormalWordsLiteralOwnerAlias=True,candidateRegionalAnnulusAndNormalWordsLiteralOwnerAlias=True,regionalSourceFaceIds=region_old_ids,regionalCandidateFaceIds=region_new_ids,physicalSourceFaceIds=owner_old_ids,physicalCandidateFaceIds=owner_new_ids))
            assert len(oldids)==len(newids)==36;oldparams[role+'_faces']=oldids;newparams[role+'_faces']=newids
        source_constraint=dict(kind='shared-interface',parameters=oldparams);source_check=check_constraint(source_constraint,*old);candidate_check=check_constraint(dict(kind='shared-interface',parameters=newparams),*new);assert source_check['passed']and candidate_check['passed']
        literal=[]
        for role,m,n in zip(['first','second'],old,new):
            same=_sheet(m,oldparams[role+'_faces'])==_sheet(n,newparams[role+'_faces']);assert same;literal.append(same)
        assert oldparams['first_faces']==newparams['first_faces']and oldparams['second_faces']==newparams['second_faces'],'This source reuse proof requires literally unchanged face indices'
        proof=dict(id=declaration['id'],first=first,second=second,firstMeshId=declaration['first_mesh_id'],secondMeshId=declaration['second_mesh_id'],pendingParametersSha256=identity(declaration['parameters']),verifiedSourceParameters=oldparams,verifiedSourceParametersSha256=identity(oldparams),candidateParameters=newparams,sourceCheck=source_check,candidateCheck=candidate_check,sourceGeometrySha256=[m.geometry_sha for m in old],candidateGeometrySha256=[m.geometry_sha for m in new],sourceAttributeSha256=[m.attributes_sha for m in old],candidateAttributeSha256=[m.attributes_sha for m in new],sourceAndCandidateOrientedFacetsAndCornerWordsLiteral=True,sourceAndCandidateFaceIndicesLiteral=True,physicalNeighborAliases=aliases,allCompletePhysicalReciprocalDrawGroupsCovered=True,requiresCandidateRebinding=False,reason='Refined source parameters also apply literally to the candidate; the constraint remains checked normally')
        proofs.append(proof)
    proofpath=OUT/'proof.json';atomic_json(proofpath,dict(schema='cystic-artery-source-bound-interface-refinement/v1',sourceCatalogueSha256=EXPECTED,proposedCatalogueSha256=PROPOSED,constraintDeclarationsSha256=digest_file(P/'constraint-declarations-v1.json'),rows=proofs,allThreeActualSourceAndCandidateChecksPass=True,noGeometryOrCatalogueChanges=True,noConstraintSkipped=True,noLiveDatabaseMutation=True,elapsedSeconds=time.monotonic()-began))
    for declaration,proof in zip(rows,proofs):
        evidence=dict(source=str(proofpath),source_sha256=digest_file(proofpath),construction_source=str(P/'audited-finite-layout.json'),construction_source_sha256=digest_file(P/'audited-finite-layout.json'),source_catalogue_sha256=EXPECTED,proposed_catalogue_sha256=PROPOSED,verification='Exact 36-triangle reciprocal finite annulus; original source and candidate coordinates and effective integer corner words verified. Source face indices also remain literal.')
        refine=dict(parameters=proof['verifiedSourceParameters'],evidence=evidence);path=OUT/('refine-'+declaration['id'][:8]+'.json');atomic_json(path,refine);refinements.append(dict(id=declaration['id'],spec=str(path),expected_pending_parameters_sha256=proof['pendingParametersSha256'],expected_refined_parameters_sha256=proof['verifiedSourceParametersSha256']))
        binding=dict(id=declaration['id'],first=proof['first'],second=proof['second'],expected_parameters_sha256=proof['verifiedSourceParametersSha256'],parameters=proof['candidateParameters'],evidence=evidence)
        # Two physical-parent directions meet the existing install validator's
        # complete exact-peer rule. The regional alias must instead reuse its
        # verified source parameters and is still checked by constraints.check.
        if not proof['physicalNeighborAliases']:
            for role,peer in [('first','second'),('second','first')]:assert binding['parameters'][role+'_faces']==drawfaces(ap[binding[role]],binding[peer])
        reuse.append(dict(id=declaration['id'],source_parameters_sha256=proof['verifiedSourceParametersSha256'],parameters=proof['verifiedSourceParameters'],candidate_check=proof['candidateCheck'],constraintRemainsCheckedNormally=True,noRebindingNeededBecauseAllSourceFaceIndicesRemainLiteral=True,regionalViewAliasesPhysicalParent=bool(proof['physicalNeighborAliases'])))
    atomic_json(OUT/'refinement-requests.json',dict(sourceCatalogueSha256=EXPECTED,refinements=refinements,liveDbEditsNotPerformed=True));atomic_json(OUT/'shared-interface-spec-patch.json',dict(shared_interface_rebindings=bindings,verified_existing_interface_parameter_reuse=reuse))
    assert digest_file(CAT)==EXPECTED and digest_file(P/'proposed-atlas.json')==PROPOSED
    print(json.dumps(dict(output=str(OUT),verifiedSourceInterfaces=3,sourceRefinementsPrepared=len(refinements),validExplicitPhysicalRebindings=len(bindings),verifiedUnchangedRegionalAliasConstraints=len(reuse),proposalUnchanged=True,noDbEdits=True,elapsedSeconds=time.monotonic()-began)))
if __name__=='__main__':main()
