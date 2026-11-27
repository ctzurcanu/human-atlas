"""Bounded deterministic morphology/context evidence for selected source vessels.

No fuzzy distance, generic parent mapping, geometry writes or runtime bindings.
"""
import collections,hashlib,json,re,sqlite3,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.ta98_authority import corrected_terms
started=time.monotonic();work=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';output=work/'ta98-inventory-v4'
prior_path=work/'ta98-inventory-v3/inventory.json';prior=json.loads(prior_path.read_bytes());previous_rules=prior['reviewed_identity_rules']
if output.exists():
 if '--refine-initial' not in sys.argv or (output/'inventory.json').exists() or (output/'reviewed-context-rules-initial-morphology.json').exists():raise SystemExit('Refusing to overwrite v4 evidence')
 (output/'reviewed-context-rules.json').rename(output/'reviewed-context-rules-initial-morphology.json');(output/'rule-selection-summary.json').rename(output/'rule-selection-initial-morphology-summary.json')
else:output.mkdir()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
with sqlite3.connect('file:'+str(ROOT/'.local-models/ta98-authority.sqlite')+'?mode=ro',uri=True) as db:
 db.row_factory=sqlite3.Row;terms=corrected_terms({r['id']:dict(r) for r in db.execute('select * from ta98')});synonyms=[dict(r) for r in db.execute('select id,synonym,lang from synonyms')]
scope={r['id']:r for r in prior['terms']};permitted={c for c,r in scope.items() if r['male_applicable'] and r['planned_representation_class']=='vessel-or-vascular-network' and not c.startswith('A12.0.')}
def canon(s):
 s=s.lower().replace('’',"'").replace('-','');s=re.sub(r'\b(left|right)\b','',s);s=re.sub(r'\s+\d+$','',s.strip());s=re.sub(r'\barteries\b','artery',s);s=re.sub(r'\bveins\b','vein',s);s=re.sub(r'\bbranches\b','branch',s);s=s.replace('esophageal','oesophageal').replace('cecal','caecal');return re.sub(r'\s+',' ',re.sub(r'[^a-z0-9]+',' ',s)).strip()
def ancestry(c):
 result=[];seen=set()
 while c in terms and c not in seen:seen.add(c);r=terms[c];result.append({'id':c,'name_en':r['name_en'],'name_la':r['name_la'],'parent_id':r['parent_id']});c=r['parent_id']
 return result
def compatible(row,c):
 arterial=row['source_group']=='Arterial system';venous=row['source_group']=='Venous system'
 if c.startswith('A12.2.') and not arterial:return False
 if c.startswith('A12.3.') and not venous:return False
 if c.startswith('A15.') and not re.search(r'arter|aort|trunk' if arterial else r'vein|venous',terms[c]['name_en'],re.I):return False
 side=row.get('source_side')
 if side in ['left','right']:
  for a in ancestry(c):
   m=re.search(r'\b(left|right)\b',a['name_en'],re.I)
   if m:return m[1].lower()==side
 return True
lex=collections.defaultdict(list)
for c in permitted:
 labels=[terms[c]['name_en']]+[r['synonym'] for r in synonyms if r['id']==c and r['lang'].lower() in ['en','eng','english']]
 for label in labels:
  lex[canon(label)].append((c,label,'authority English label/synonym'))
  for a in ancestry(c)[1:]:
   if re.search(r'\b(?:artery|arteries|vein|veins|trunk|plexus|sinus)\b',a['name_en']):lex[canon(label+' of '+a['name_en'])].append((c,label+' of '+a['name_en'],'authority label plus explicit vascular ancestor'))
new=[];external=list(previous_rules['cross_root_out_of_denominator']);external_ids={r['source_id'] for r in external};deferred=[]
for row in prior['source_mesh_instances']:
 if row['status']!='unmapped' or row['source_id'] in external_ids:continue
 name=row['source_name'];key=canon(name);hits=[h for h in lex.get(key,[]) if compatible(row,h[0])];codes={h[0] for h in hits};chosen=next(iter(codes)) if len(codes)==1 else None;reason='Exact authority English label/synonym or explicit vascular-ancestor expansion after only side removal, artery/vein/branch number agreement, documented oesophageal/caecal orthography and a trailing source occurrence number. Source laterality and complete source hierarchy remain explicit.';rule_id='literal-morphology-context';labels=hits
 # Reviewed family-specific nomenclature. The authority English/Latin term
 # and actual source-parent path are retained on each individual rule.
 manual=None;parents=' / '.join(p['name'] for p in row['source_hierarchy']).lower()
 if key in ['inferior lingular segmental artery','superior lingular segmental artery'] and 'left pulmonary artery' in parents:manual='A12.2.01.212' if key.startswith('inferior') else 'A12.2.01.213';rule_id='named-lingular-segmental-descriptor';reason='Exact inferior/superior lingular arterial identity under left pulmonary artery; source adds the descriptive word segmental. Authority English and Latin name the same lingular artery, not a broader lobar parent.'
 elif key in ['inferior labial artery','superior labial artery'] and 'facial artery' in parents:manual='A12.2.05.025' if key.startswith('inferior') else 'A12.2.05.026';rule_id='named-labial-artery-branch';reason='Authority English calls this named facial branch; its Latin a. labialis inferior/superior explicitly identifies the source inferior/superior labial artery. Source facial ancestry agrees.'
 elif key=='calcarine artery' and 'posterior cerebral artery' in parents:manual='A12.2.07.102';rule_id='named-calcarine-artery-branch';reason='Source calcarine artery is the specifically named calcarine branch (Latin r. calcarinus), within the posterior cerebral/medial occipital authority ancestry; not an arbitrary posterior cerebral parent.'
 elif key in ['medial posterior choroidal artery','lateral posterior choroidal artery'] and 'posterior cerebral artery' in parents:manual='A12.2.07.091' if key.startswith('medial') else 'A12.2.07.092';rule_id='named-posterior-choroidal-word-order';reason='Source medial/lateral posterior choroidal arterial name and posterior cerebral ancestry correspond to authority posterior medial/lateral choroidal branches and matching Latin rr. choroidei posteriores mediales/laterales. Record as named branch collection; no parent vessel substituted.'
 elif key=='deep femoral vein' and 'deep veins' in parents:manual='A12.3.11.024';rule_id='deep-femoral-profunda-femoris';reason='Source deep femoral vein is the named profunda femoris vein; Latin v. profunda femoris fixes the femoral/deep identity, and source is in the deep veins of thigh/leg, not a generic femoral-vein parent.'
 elif re.fullmatch(r'(posterior|lateral) branch of superior thyroid artery',key) and 'superior thyroid artery' in parents:manual='A12.2.05.008' if key.startswith('posterior') else 'A12.2.05.009';rule_id='qualified-superior-thyroid-glandular-branch';reason='The exact posterior/lateral branch of the source superior thyroid artery corresponds to its sole authority posterior/lateral glandular branch, with matching English/Latin direction and parent. No unspecific arterial branch is accepted.'
 elif re.fullmatch(r'(medial|lateral) dorsal cutaneous branch of posterior intercostal artery',key):manual='A12.2.11.009' if key.startswith('medial') else 'A12.2.11.010';rule_id='intercostal-dorsal-cutaneous-qualified-parent';reason='Source explicitly says medial/lateral dorsal cutaneous branch of a numbered posterior intercostal artery. Authority has medial/lateral cutaneous branches directly under its dorsal branch of posterior intercostal arteries. The source dorsal modifier expresses the authority parent; exact source occurrence IDs remain distinct.'
 if manual:
  assert manual in permitted and compatible(row,manual);chosen=manual;labels=[(manual,terms[manual]['name_en'],'reviewed exact family English/Latin and parent qualification')]
 # Named first/second posterior intercostal arteries have their own TA98
 # identities under the supreme intercostal artery. Do not inflate aortic
 # posterior-intercostal collections with these different named origins.
 occurrence=re.search(r'\s(\d+)$',name);occurrence=int(occurrence[1]) if occurrence else None
 if key=='posterior intercostal artery' and occurrence in [1,2]:
  chosen='A12.2.08.062' if occurrence==1 else 'A12.2.08.063';rule_id='explicit-first-second-posterior-intercostal';reason='The source numbered first/second posterior intercostal artery corresponds to the explicit authority first/second posterior intercostal artery under supreme intercostal artery. Keep this identity separate from the thoracic-aorta posterior collection.';labels=[(chosen,terms[chosen]['name_en'],'explicit source occurrence number and authority ordinal')]
 elif key=='dorsal branch of posterior intercostal artery' and occurrence==2:
  chosen='A12.2.08.064';rule_id='explicit-second-posterior-intercostal-dorsal';reason='Source dorsal branch of the second posterior intercostal artery matches authority dorsal branches directly under second posterior intercostal artery. Number and named parent are explicit.';labels=[(chosen,terms[chosen]['name_en'],'explicit source occurrence number and authority parent')]
 elif key in ['dorsal branch of posterior intercostal artery','medial dorsal cutaneous branch of posterior intercostal artery','lateral dorsal cutaneous branch of posterior intercostal artery','collateral branch of posterior intercostal artery','lateral cutaneous branch of posterior intercostal artery'] and occurrence in [1,2]:
  deferred.append({'source_id':row['source_id'],'source_name':name,'reason':'First/second posterior intercostal origin differs from the aortic collection. No exact indexed child is substituted by a generic posterior intercostal parent.'});continue
 elif key=='lateral cutaneous branch of posterior intercostal artery' and occurrence is not None and 3<=occurrence<=11 and row['source_parent_name'].lower().startswith('lateral cutaneous branches'):
  chosen='A12.2.11.018';rule_id='direct-intercostal-lateral-cutaneous-context';reason='Source is explicitly in the direct lateral cutaneous branch collection of posterior intercostal arteries, not the separately named dorsal cutaneous collection. Match the authority direct lateral cutaneous branch under posterior intercostal arteries; retain source numbers 3–11.';labels=[(chosen,terms[chosen]['name_en'],'exact direct versus dorsal source-parent distinction')]
 elif key=='anterior intercostal artery' and occurrence is not None and 1<=occurrence<=6:
  chosen='A12.2.08.039';rule_id='named-internal-thoracic-anterior-intercostal-collection';reason='Source anterior intercostal arteries numbered 1–6 correspond to authority anterior intercostal branches of internal thoracic artery (Latin rr. intercostales anteriores). Lower source occurrences 7–9 are not forced into this direct-parent collection.';labels=[(chosen,terms[chosen]['name_en'],'explicit numbered anterior intercostal vessel collection and Latin branch terminology')]
 if not chosen:
  # Additional explicit kidney vessel terms are kept outside the fixed
  # historical denominator. Do not relabel capsular veins as arterial tissue.
  m=re.fullmatch(r'Interlobar (arteries|veins) of (left|right) kidney',name)
  if m:
   c='A08.1.03.002' if m[1]=='arteries' else 'A08.1.04.002';t=terms[c];external.append({'rule_id':'explicit-intrarenal-interlobar-organ-context','source_id':row['source_id'],'source_name':name,'source_hierarchy':row['source_hierarchy'],'authority':{'id':c,'name_en':t['name_en'],'name_la':t['name_la'],'parent_id':t['parent_id'],'parent_name':terms[t['parent_id']]['name_en']},'authority_ancestry':ancestry(c),'rationale':'Exact interlobar arteries/veins label qualified by kidney and corresponding same-side renal vessel ancestry; authoritative intrarenal vessel term is outside the historical denominator.','out_of_current_denominator':True,'candidate_identity_only':True,'runtime_binding':False,'anatomical_acceptance':False})
  continue
 if chosen in {r['authority']['id'] for r in previous_rules['in_denominator_rules'] if r['source_id']==row['source_id']}:continue
 t=terms[chosen];number=re.search(r'\s(\d+)$',name)
 item={'rule_id':rule_id,'source_id':row['source_id'],'source_name':name,'source_hierarchy':row['source_hierarchy'],'authority':{'id':chosen,'name_en':t['name_en'],'name_la':t['name_la'],'parent_id':t['parent_id'],'parent_name':terms[t['parent_id']]['name_en']},'authority_ancestry':ancestry(chosen),'literal_normalized_source':key,'matched_authority_expressions':[{'text':label,'method':method} for c,label,method in labels if c==chosen],'rationale':reason,'source_numbered_occurrence':int(number[1]) if number else None,'numbered_collection_identity_only':bool(number),'collective_completeness_anatomically_accepted':False,'out_of_current_denominator':False,'candidate_identity_only':True,'runtime_binding':False,'anatomical_acceptance':False};new.append(item)
all_rules=previous_rules['in_denominator_rules']+new
document={'schema':'reviewed-blood-source-context-rules/v2','inventoryStage':'v4','sourceInventoryV3Sha256':sha(prior_path),'denominatorStable':True,'in_denominator_rules':all_rules,'cross_root_out_of_denominator':external,'deferred_hepatic_numbered_territories':previous_rules['deferred_hepatic_numbered_territories'],'deferred_numbered_origin_mismatches':deferred,'policy':'Only exact morphological/orthographic or explicitly reviewed named anatomical-family rules. Numbered members remain per-occurrence source identity controls; aggregate completeness and anatomical fit are not claimed. First/second posterior intercostal arteries retain their explicit supreme-intercostal authority identities. No generic parent identity is introduced.'}
(output/'reviewed-context-rules.json').write_text(json.dumps(document,indent=2)+'\n')
counts=collections.Counter(r['authority']['id'] for r in new);old_terms={r['id'] for r in prior['terms'] if r['donor_candidate_source_ids']};summary={'newRuleInstances':len(new),'newDistinctIdentityTerms':len(set(counts)-old_terms),'newCodes':[{'term':c,'name':terms[c]['name_en'],'instances':n,'additionalToV3':c not in old_terms} for c,n in sorted(counts.items())],'crossRootInstances':len(external),'crossRootDistinctTerms':len({r['authority']['id'] for r in external}),'elapsedSeconds':time.monotonic()-started};(output/'rule-selection-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
