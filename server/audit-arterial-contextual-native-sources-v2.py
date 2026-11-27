"""Strict parent-qualified source retrieval; no anatomy/name fuzzing accepted."""
from pathlib import Path
from collections import defaultdict
import sys,time,json,re
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1';OUT=B/'context-qualified-scan-v2';EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'

def normal(s):
    s=str(s).lower().replace('_',' ');s=re.sub(r'\s*\((?:left|right)\)\s*$|\.(?:l|r)$','',s);s=re.sub(r'\s+',' ',s).strip(' .')
    return s.replace('arteries','artery').replace('branches','branch').replace('veins','vein')

def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir();prior=json.loads((B/'metadata-source-scan.json').read_text());terms=prior['terms'];reader=AtlasReader();hits=[];catalogues={};records={}
    for cat in prior['sourceCatalogues']:
        path=Path(cat['path']);a=json.loads(path.read_text());assert digest_file(path)==cat['sha256'];catalogues[str(path)]=a;records[str(path)]={p['id']:p for p in a['parts']}
        for p in a['parts']:
            if p.get('system')!='arterial':continue
            raw=normal(p['name']);variants={raw}
            if raw.startswith(('left ','right ')):variants.add(raw.split(' ',1)[1])
            for term in terms:
                if term['planned_representation_class']!='vessel-or-vascular-network':continue
                desired=normal(term['name_en']);parent=normal(term['parent_name']);expected={desired,desired+' of '+parent,desired+' of left '+parent,desired+' of right '+parent}
                generic=desired in {'anterior branch','posterior branch','superficial branch','deep branch','muscular branch','glandular branch','cutaneous branch','spinal branch','perforating branch','ascending branch','descending branch'}
                matched=variants&expected
                if not matched:continue
                if generic and matched&{desired}:continue
                if not generic and len(variants&{desired})==0 and len(matched)==0:continue
                hits.append(dict(term=term['id'],name=term['name_en'],parent=term['parent_name'],part_id=p['id'],part_name=p['name'],catalogue=str(path),catalogue_sha256=cat['sha256'],identity_evidence='Strict singular/plural and explicit source parent qualifier; still requires source review',explicit_parent_qualifier=any(' of 'in x for x in matched),license=p.get('license'),suppressed=bool(p.get('suppressed')),geometryInstance=p.get('geometryInstance'),sourceId=p.get('sourceId'),chunk_url=a['chunks'][p['chunk']]['url'],triangles=p.get('triangles',p.get('indexCount',0)//3),nativeCurrentCatalogue=str(path)==prior['sourceCatalogues'][0]['path']))
    unique={(r['term'],r['catalogue'],r['part_id']):r for r in hits};hits=list(unique.values());hits.sort(key=lambda r:(not r['nativeCurrentCatalogue'],r['term'],r['part_id']))
    # Independent native course families are worth a bounded topology audit
    # even when their detailed regions do not yet exist as separate objects.
    currentpath=Path(prior['sourceCatalogues'][0]['path']);assert digest_file(currentpath)==EXPECTED;current=catalogues[str(currentpath)];cp=records[str(currentpath)];chosen=['ZA:Internal carotid artery.l','ZA:Internal carotid artery.r','ZA:Vertebral artery.l','ZA:Vertebral artery.r','ZA:Internal thoracic artery.l','ZA:Internal thoracic artery.r','ZA:Posterior intercostal arteries.l','ZA:Posterior intercostal arteries.r']
    native=[]
    for pid in chosen:
        if pid not in cp:continue
        m=reader.read_part(current,cp[pid],frame='model:local-ta98:assembled-metres',catalogue_path=currentpath);q=topology(m);native.append(dict(part_id=pid,name=cp[pid]['name'],license=cp[pid].get('license'),geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,bank_sha256=m.source['bank_sha256'],bank_path=m.source.get('bank_path'),source_catalogue_sha256=EXPECTED,triangles=len(m.faces),bounds=[m.triangles.min((0,1)).tolist(),m.triangles.max((0,1)).tolist()],topology=q,not_a_distinct_region_until_bounded_by_actual_landmarks=True))
    result=dict(schema='strict-parent-qualified-arterial-source-scan/v2',sourceCatalogueSha256=EXPECTED,candidates=hits,termsWithContextualCandidates=len({r['term']for r in hits}),nativeCourseMeshes=native,noTermBoundToWholeCarrier=True,noGeometryOrCatalogueChanged=True,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'proof.json',result);print(json.dumps(dict(contextualCandidates=len(hits),termsWithCandidates=result['termsWithContextualCandidates'],candidateTerms=sorted({r['term']+': '+r['name']for r in hits}),nativeCourseSummary=[dict(id=r['part_id'],triangles=r['triangles'],components=len(r['topology']['components']),openEdges=r['topology']['open_edges'])for r in native],elapsedSeconds=result['elapsedSeconds'])))
if __name__=='__main__':main()
