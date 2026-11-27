"""Compact source availability only; no fitting, registration or mesh changes."""
from pathlib import Path
from collections import defaultdict,Counter
import json,re,sqlite3,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology
T=time.monotonic();census_path=ROOT/'reports/ta98-blood-circulation-inventory-installed-v7.json';census=json.loads(census_path.read_text());missing={r['id']:r for r in census['terms']if r['male_applicable']and r['representation_status']=='missing-direct-representation'};current_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';current=json.loads(current_path.read_text());live={p.get('ta98Term')for p in current['parts']if not p.get('suppressed')and p.get('ta98Kind')=='exact'};reader=AtlasReader()
def normalized(text):
    text=re.sub(r'\b(left|right)\b|\(([^)]*)\)',' ',text.lower());return re.sub(r'[^a-z0-9]+',' ',text).strip()
lookup=defaultdict(set)
for term,row in missing.items():lookup[normalized(row['name_en'])].add(term)
with sqlite3.connect('file:'+str(ROOT/'.local-models/ta98-authority.sqlite')+'?mode=ro',uri=True)as db:
    for term,name,lang in db.execute('SELECT id,synonym,lang FROM synonyms'):
        if term in missing and lang in ['en','EN','eng']and any(s in name.lower()for s in ['artery','arteries','vein','veins','sinus','plexus']):lookup[normalized(name)].add(term)
    fma={term:str(fid)for term,fid in db.execute('SELECT id,fma_id FROM ta98')}
catalogues=[('HRA',ROOT/'public/models/atlas-hra-female.json','CC BY 4.0'),('BodyParts3D',ROOT/'public/models/atlas.json','CC BY 4.0')];rows=[];counts=[];already=[];ambiguous=[]
for label,path,license_name in catalogues:
    a=json.loads(path.read_text());parts=[p for p in a['parts']if p.get('system')in ['arterial','venous']and p.get('indexCount',0)>=3];direct_fma=defaultdict(list)
    for term in missing:
        for p in parts:
            if str(p.get('conceptId','')).replace('FMA','')==fma.get(term):direct_fma[term].append(p['id'])
    matched=0
    for p in parts:
        key=normalized(p['name']);codes=lookup.get(key,set())
        if not codes:continue
        if len(codes)!=1:ambiguous.append(dict(source=label,part_id=p['id'],possible_terms=sorted(codes)));continue
        term=next(iter(codes));matched+=1
        if term in live:already.append(dict(term=term,name=missing[term]['name_en'],part_id=p['id'],reason='Now directly present after the v7 census; excluded'));continue
        # Named individual vascular source, never an entire mixed-system carrier.
        assert not any(s in key for s in ['arterial system','venous system','blood vessels','circulatory system'])
        m=reader.read_part(a,p,frame='source:'+label+':catalogue-metres',catalogue_path=path);t=topology(m);part=dict(source=label,catalogue=str(path.relative_to(ROOT)),catalogue_sha256=digest_file(path),source_sex=a.get('sex','male'),license=license_name,part_id=p['id'],name=p['name'],index_count=p['indexCount'],triangles=len(m.faces),geometry_sha256=m.geometry_sha,bank_path=m.source['bank_path'],bank_sha256=m.source['bank_sha256'],source_bounds_m=m.bounds.tolist(),connected_components=len(t['components']),open_edges=t['open_edges'],identity_match='Unambiguous normalized full English vessel name or authoritative English synonym; registration is not inferred',term=term,ta98_name=missing[term]['name_en'],ta98_parent=missing[term]['parent_id'],ta98_parent_name=missing[term]['parent_name'],wall_lumen_and_attachment_verified=False)
        rows.append(part)
    counts.append(dict(source=label,licensed_vascular_mesh_records=len(parts),strict_name_or_synonym_matches_before_live_exclusion=matched,direct_missing_FMA_matches=direct_fma,available_after_live_exclusion=sum(r['source']==label for r in rows)))
families=[]
for term in sorted({r['term']for r in rows}):
    source=[r for r in rows if r['term']==term];heart=term.startswith('A12.2.03.')or term.startswith('A12.3.01.');inflight=term=='A12.3.12.012'
    barriers=['Source is an outer vessel envelope; no finite wall, inner lumen or recipient mouth is established.','Male source registration and actual parent/recipient interfaces require a source-bound construction contract; fixed existing anatomy must remain unchanged.']
    if term=='A12.3.06.111':barriers+=['Fit independently to the fixed male optic nerve/retinal and superior-ophthalmic-vein interfaces. Do not infer a microscopic retinal capillary bed from this vein mesh.']
    if inflight:barriers+=['Cystic venous donor is already in the current gallbladder/portal task; keep it queued behind the unresolved native portal host interface.']
    if heart:barriers+=['Held outside the requested away-from-heart donor shortlist. Any future coronary/cardiac-vein work must preserve the fitted heart and existing fixed surface-vessel courses.']
    families.append(dict(term=term,name=source[0]['ta98_name'],mesh_count=len(source),triangles=sum(s['triangles']for s in source),source_parts=[s['part_id']for s in source],recommended_new_batch=not heart and not inflight,already_in_active_task=inflight,held_heart_surface_vessel=heart,contract_barriers=barriers))
out=ROOT/'.local-models/geometry-workbench/circulatory-completion/missing-male-cc-blood-donor-availability-v1.json';report=dict(schema='missing-male-blood-licensed-donor-availability/v1',census=str(census_path.relative_to(ROOT)),census_sha256=digest_file(census_path),male_missing_terms_in_v7=len(missing),licensing_evidence=dict(path='ATTRIBUTION.md',sha256=digest_file(ROOT/'ATTRIBUTION.md'),hra_lines='23–31',bodyparts3d_lines='33–37'),method='Local indexed mesh availability using only specific full vessel names, authoritative English synonyms and direct FMA identity; no fuzzy name or shared-carrier promotion.',catalogue_counts=counts,families=families,sources=rows,now_installed_excluded=already,ambiguous_names_excluded=ambiguous,limits=['Zero strict BodyParts3D matches is not proof that no reusable subbranch can be extracted from a source network.','HRA donor frame is unregistered to the male; no inherited anatomical acceptance.','The two heart-surface vessels are listed only as held availability evidence, not recommended next work.'],geometry_or_catalogues_modified=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T);atomic_json(out,report);print(json.dumps(dict(output=str(out),strict_available_terms=len(families),actual_meshes=len(rows),recommended_families=[f for f in families if f['recommended_new_batch']],held_or_inflight=[f['term']for f in families if not f['recommended_new_batch']],counts=counts,elapsed_seconds=report['elapsed_seconds'])))
