"""Exact authoritative FMA correspondence to licensed cerebral source records."""
from pathlib import Path
import sys,time,json,sqlite3
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file

CAT=ROOT/'public/models/atlas.json';AUTH=ROOT/'.local-models/ta98-authority.sqlite';CENSUS=ROOT/'reports/ta98-blood-circulation-inventory-installed-v11.json'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1/bp3d-cerebral-term-opportunities-v8'

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False);atlas=json.loads(CAT.read_text());parts={p['id']:p for p in atlas['parts']};concepts={c['id']:c for c in atlas['concepts']};census=json.loads(CENSUS.read_text());status={r['id']:r for r in census['terms']}
    db=sqlite3.connect('file:'+str(AUTH)+'?mode=ro',uri=True);db.row_factory=sqlite3.Row
    terms=[dict(r)for r in db.execute("select id,name_en,name_la,parent_id,parent_name,fma_id,type_of_entity,immaterial from ta98 where id like 'A12.2.07.%' order by id")]
    rows=[]
    for t in terms:
        cid='FMA'+str(t['fma_id']);c=concepts.get(cid)
        if not c:continue
        ids=[i for i in c.get('elements',[])if i in parts and parts[i].get('system')=='arterial']
        if not ids:continue
        r=status.get(t['id'],{});direct=r.get('resolved_direct_part_ids',[])
        # Explicitly retain broad/source-aggregate cases as audit-only; they
        # cannot manufacture a new local branch from an unchanged carrier.
        specific=t['type_of_entity']=='material entity'and not t['immaterial']and len(ids)<=20 and t['id']not in ['A12.2.07.001','A12.2.07.080']
        rows.append(dict(authority=t,sourceConceptId=cid,originalSourceConceptName=c['name'],identityEvidence='Authoritative TA98 fma_id equals original official BodyParts3D concept ID exactly; source element membership remains explicit',sourceParts=[dict(id=i,originalName=parts[i]['name'],originalFmaId=parts[i]['conceptId'],triangles=parts[i]['indexCount']//3,record=parts[i])for i in ids],currentDirectPartIds=direct,currentlyMissingDirectRepresentation=not bool(direct),specificBoundedSourceFamilyCandidate=specific,notAnatomicallyAccepted=True,neededChecks=['Exact original source common-frame replay and fixed brain/skeleton placement','Whole-source topology and actual complete course','Actual parent and child root positions, shared mouths and lumen continuity','No duplicate active physical owners; exact source-face authority ownership','Native surface and sections before any local installation']))
    candidates=[r for r in rows if r['currentlyMissingDirectRepresentation']and r['specificBoundedSourceFamilyCandidate']]
    out=dict(schema='bp3d-cerebral-exact-authority-source-opportunities/v8',sourceCatalogueSha256=digest_file(CAT),authoritySha256=digest_file(AUTH),censusSha256=digest_file(CENSUS),sourceProvenance=atlas['provenance'],licence='BodyParts3D CC BY 4.0; ATTRIBUTION.md',authoritativeCerebralTerms=len(terms),exactFmaMatchedSourceFamilies=len(rows),missingSpecificCandidateFamilies=len(candidates),missingSpecificCandidateTermIds=[r['authority']['id']for r in candidates],rows=rows,missingTermsWithoutAnExactOriginalSourceFmaConcept=[t for t in terms if not status.get(t['id'],{}).get('resolved_direct_part_ids')and not any(r['authority']['id']==t['id']for r in rows)],noLooseNameMatching=True,noBroadCarrierAliasPromotion=True,noGlbDecoded=True,noGeometryOrLiveMetadataChanged=True,elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'proof.json',out)
    print(json.dumps(dict(output=str(OUT/'proof.json'),matchedFamilies=len(rows),missingSpecificCandidates=[dict(term=r['authority']['id'],name=r['authority']['name_en'],sourceConcept=r['sourceConceptId'],sourceParts=[p['id']for p in r['sourceParts']])for r in candidates],elapsedSeconds=out['elapsedSeconds'])))
if __name__=='__main__':main()
