"""Read-only exact terminology/source retrieval; no carrier is accepted here."""
from pathlib import Path
from collections import defaultdict,Counter
import sys,time,json,re,sqlite3
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file

OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1'
REPORT=ROOT/'reports/ta98-blood-circulation-inventory-installed-v11.json'
EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'

def normal(s):
    s=str(s).lower().replace('_',' ').replace('arteria ','').replace('arteriae ','')
    s=re.sub(r'\s*\((?:left|right)\)\s*$','',s)
    s=re.sub(r'\.(?:l|r)$','',s)
    s=re.sub(r'^(?:left|right)\s+|\s+(?:left|right)$','',s)
    s=re.sub(r'\s+',' ',s).strip(' .')
    return s

def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir();report=json.loads(REPORT.read_text());current=ROOT/report['source_catalogue'];assert digest_file(current)==EXPECTED;terms={r['id']:r for r in report['terms']};missing={c:r for c,r in terms.items()if c.startswith('A12.2.')and r['male_applicable']and not r['resolved_direct_part_ids']}
    def ancestors(code):
        out=[]
        while code in terms and code not in out:out.append(code);code=terms[code]['parent_id']
        return out
    excluded=set()
    for code,r in missing.items():
        names=' | '.join(terms[c]['name_en'].lower()for c in ancestors(code))
        if any(x in names for x in ['facial artery','lingual artery','coronary artery','pulmonary trunk','pulmonary artery','aortic root','ascending aorta','aortic arch']):excluded.add(code)
    scoped={c:r for c,r in missing.items()if c not in excluded}
    names=defaultdict(set)
    for c,r in scoped.items():
        for n in [r['name_en'],r['name_la']]:names[normal(n)].add(c)
    with sqlite3.connect('file:'+str(ROOT/'.local-models/ta98-authority.sqlite')+'?mode=ro',uri=True)as db:
        for code,syn,lang in db.execute('select id,synonym,lang from synonyms'):
            if code in scoped and lang.lower()in ['en','eng','english']:names[normal(syn)].add(code)
    paths=[current,ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json',ROOT/'public/models/atlas-z-anatomy.json',ROOT/'public/models/atlas-male-complete.json',ROOT/'public/models/atlas-hra-female.json',ROOT/'public/models/atlas.json']
    candidates=defaultdict(list);catalogues=[]
    generic={'muscular branches','cutaneous branches','sternal branches','glandular branches','articular branches','terminal branch','anterior branch','posterior branch','superficial branch','deep branch','dorsal branches','perforating branches'}
    for path in paths:
        if not path.exists():continue
        a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};exact=defaultdict(set)
        for p in parts.values():
            if p.get('ta98Term'):exact[p['ta98Term']].add(p['id'])
        for c in a.get('concepts',[]):
            code=c.get('ta98Term')
            if code and c.get('ta98Kind','exact')=='exact':
                for pid in c.get('elements',[]):
                    if pid in parts:exact[code].add(pid)
        for b in a.get('ta98TermBindings',[]):
            if b.get('kind','exact')=='exact':
                for pid in b.get('partIds',[]):
                    if pid in parts:exact[b['term']].add(pid)
        inverse=defaultdict(set)
        for c,ids in exact.items():
            for pid in ids:inverse[pid].add(c)
        catalogues.append(dict(path=str(path),sha256=digest_file(path),parts=len(parts),arterialParts=sum(p.get('system')=='arterial'for p in parts.values())))
        seen=set()
        for pid,p in parts.items():
            named=names.get(normal(p.get('name','')),set());mapped={c for c in inverse[pid]if c in scoped}
            for code in named|mapped:
                key=(code,pid)
                if key in seen:continue
                seen.add(key);other=sorted(inverse[pid]-{code});method='explicit-exact-source-term'if code in mapped else'strict-authority-name-or-synonym';ambiguous=normal(scoped[code]['name_en'])in generic
                candidates[code].append(dict(part_id=pid,name=p.get('name'),catalogue=str(path),catalogue_sha256=catalogues[-1]['sha256'],system=p.get('system'),license=p.get('license'),provenance=p.get('provenance'),suppressed=bool(p.get('suppressed')),defaultHidden=bool(p.get('defaultHidden')),triangles=p.get('triangles',p.get('indexCount',0)//3),bounds=p.get('bounds'),geometryInstance=p.get('geometryInstance'),sourceId=p.get('sourceId'),chunk_url=a['chunks'][p['chunk']]['url']if'chunk'in p else None,identity_evidence=method,other_exact_terms_using_same_source=other,generic_name_needs_parent_branch_identity=ambiguous,coarse_or_shared_carrier_not_distinct=len(other)>0 and method!='explicit-exact-source-term',not_accepted=True))
    rows=[]
    for code,r in scoped.items():
        row={k:r[k]for k in ['id','name_en','name_la','parent_id','parent_name','planned_representation_class','variant','bilaterality']};row['ancestry']=[dict(id=c,name=terms[c]['name_en'])for c in ancestors(code)];row['candidates']=candidates[code];rows.append(row)
    families=[]
    for b in report['branches']:
        if not b['root'].startswith('A12.2.')or b['root']=='A12.2.00.001':continue
        codes={c for c in scoped if b['root']in ancestors(c)};available=[c for c in codes if any(not p['generic_name_needs_parent_branch_identity']and not p['coarse_or_shared_carrier_not_distinct']for p in candidates[c])]
        if available:families.append(dict(root=b['root'],name=b['name'],missingInScope=len(codes),termsWithExactRetrievalCandidates=len(available),terms=sorted(available)))
    families.sort(key=lambda f:(-f['termsWithExactRetrievalCandidates'],f['missingInScope']))
    proof=dict(schema='missing-arterial-source-opportunities/v1',sourceCatalogueSha256=EXPECTED,censusSha256=digest_file(REPORT),missingMaleArterialTerms=len(missing),excludedFacialLingualHeartRootTerms=len(excluded),scopedMissingTerms=len(scoped),termsWithAnyCandidate=sum(bool(v)for v in candidates.values()),sourceCatalogues=catalogues,rankedFamilies=families,terms=rows,retrievalIsNotAnatomicalAcceptance=True,noNewDirectRepresentationClaimed=True,noGeometryCatalogueOrDbChanged=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'metadata-source-scan.json',proof)
    print(json.dumps(dict(scopedMissing=len(scoped),termsWithCandidates=len(candidates),topFamilies=families[:12],elapsedSeconds=proof['elapsedSeconds'])))
if __name__=='__main__':main()
