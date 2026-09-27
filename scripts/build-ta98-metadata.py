#!/usr/bin/env python3
"""Build the viewer's compact TA98 concordance from the published TA98 SQLite extract.

Usage: python3 scripts/build-ta98-metadata.py [path/to/ta98.sqlite]
The default downloads mhalle/ta98-sqlite's extract of the IFAA TA98 pages.
Unmatched or ambiguous atlas concepts deliberately receive no TA/THA code.
"""
import collections
import json
import pathlib
import re
import sqlite3
import sys
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE = 'https://raw.githubusercontent.com/mhalle/ta98-sqlite/master/db/ta98.sqlite'
CATALOGS = [
 'public/models/atlas-male-complete.json', 'public/models/atlas.json',
 'public/models/atlas-hra-female.json', 'public/models/atlas-embryo.json',
 'public/models/atlas-cell.json', '.local-models/reference.json',
 '.local-models/male.json', '.local-models/female.json',
]

def cleaned(name):
 name = name.lower().split(' · ')[0].strip()
 name = re.sub(r'\.[oe]\d*[lr]$', '', name)
 name = re.sub(r'\s+', ' ', name).strip()
 if name == 'skin of body':
  return 'skin'
 return name

def normalized(name):
 name = cleaned(name)
 name = re.sub(r'\s*\((?:left|right)\)$', '', name)
 name = re.sub(r'^(?:left|right)\s+', '', name)
 name = name.strip()
 if name.startswith('(') and name.endswith(')') and name.count('(') == 1:
  name = name[1:-1]
 return name

def variants(name):
 """Conservative atlas spelling variants; never use fuzzy similarity for TA IDs."""
 base = normalized(name)
 candidates = [base]
 aliases = {
  'external abdominal oblique muscle':'external oblique',
  'internal abdominal oblique muscle':'internal oblique',
 }
 if base in aliases:
  candidates.append(aliases[base])
 substitutions = [
  (r'\babducens\b', 'abducent'),
  (r'\bintervertebral disk\b', 'intervertebral disc'),
  (r'\billiacus\b', 'iliacus'),
 ]
 for pattern, replacement in substitutions:
  changed = re.sub(pattern, replacement, base)
  if changed != base:
   candidates.append(changed)
 for candidate in tuple(candidates):
  short = re.sub(r'\s*\((?:[ivx]+|[cclmt][ivx]+|b[ivx]+(?:\+b[ivx]+)?)\)$', '', candidate).strip(' ()')
  if short != candidate:
   candidates.append(short)
 for candidate in tuple(candidates):
  if candidate.endswith(' muscles'):
   candidates.append(candidate[:-1])
  for suffix in (' muscle', ' bone'):
   if candidate.endswith(suffix):
    candidates.append(candidate[:-len(suffix)])
  if candidate.startswith('female '):
   candidates.append(candidate[7:])
  if candidate.startswith('male '):
   candidates.append(candidate[5:])
 return list(dict.fromkeys(candidates))

db_file = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path('/tmp/human-atlas-ta98.sqlite')
if not db_file.exists():
 urllib.request.urlretrieve(SOURCE, db_file)
db = sqlite3.connect(db_file)
db.row_factory = sqlite3.Row
terms = {row['id']: dict(row) for row in db.execute('select id, source_id, name_en, name_la, entity_id_number, fma_id from ta98')}
ancestors = collections.defaultdict(list)
for row in db.execute('select id, ancestor_id from hierarchy order by hierarchy_level desc'):
 ancestors[row['id']].append(row['ancestor_id'])
by_primary = collections.defaultdict(set)
by_synonym = collections.defaultdict(set)
by_base_primary = collections.defaultdict(set)
by_base_synonym = collections.defaultdict(set)
by_fma = collections.defaultdict(set)
for row in terms.values():
 by_primary[cleaned(row['name_en'])].add(row['id'])
 by_base_primary[normalized(row['name_en'])].add(row['id'])
 if row['fma_id']:
  by_fma[row['fma_id']].add(row['id'])
for row in db.execute("select id, synonym from synonyms where lang='en'"):
 by_synonym[cleaned(row['synonym'])].add(row['id'])
 by_base_synonym[normalized(row['synonym'])].add(row['id'])

def record(term):
 return {'ta98': term['source_id'] or term['id'], 'tha': term['entity_id_number'],
         'fma': f"FMA:{term['fma_id']}" if term['fma_id'] else None,
         'latin': term['name_la']}

def resolve(name, concept_id=None, system=None):
 fma_id = re.fullmatch(r'FMA:?([0-9]+)', concept_id or '')
 if fma_id and len(by_fma[fma_id[1]]) == 1:
  return record(terms[next(iter(by_fma[fma_id[1]]))])
 if system == 'attachments':
  # A source origin/insertion marker is not the muscle itself. Its mesh keeps
  # its atlas ID, while the related muscle carries the TA98 identifier.
  return {'ta98':None, 'tha':None, 'fma':None, 'latin':None}
 chapter = {'skeletal':'A02', 'connective':'A03', 'muscular':'A04',
            'digestive':'A05', 'respiratory':'A06', 'urinary':'A08',
            'reproductive':'A09', 'endocrine':'A11', 'cardiac':'A12',
            'arterial':'A12', 'venous':'A12', 'lymphatic':'A13',
            'nervous':'A14', 'sensory':'A15', 'integumentary':'A16',
            'regions':'A01', 'fascia':'A04'}.get(system)
 lookups = [(by_primary,cleaned(name)),(by_synonym,cleaned(name))]
 for candidate in variants(name):
  lookups.extend([(by_base_primary,candidate),(by_base_synonym,candidate)])
 for index, lookup in lookups:
  candidates = index.get(lookup, set())
  if len(candidates) > 1:
   preferred = {id for id in candidates if id.startswith(chapter or '!')}
   if len(preferred) == 1:
    candidates = preferred
  if len(candidates) == 1:
   return record(terms[next(iter(candidates))])
 return {'ta98':None, 'tha':None,
         'fma':f'FMA:{fma_id[1]}' if fma_id else None, 'latin':None}

concepts = {}
concept_paths = {}
source_groups = set()
for catalog in CATALOGS:
 path = ROOT / catalog
 if not path.exists():
  continue
 atlas = json.loads(path.read_text())
 names = {concept['id']: concept['name'] for concept in atlas.get('concepts', [])}
 for part in atlas['parts']:
  source_groups.update(group for group in part.get('groups', []) if not re.match(r'^\d+:\s', group))
  if part.get('suppressed') or part['conceptId'] in concepts:
   continue
  cid = part['conceptId']
  name = names.get(cid) or part['name']
  system = 'regions' if any(re.match(r'^9: regions of human body$', group, re.I) for group in part.get('groups', [])) else part['system']
  match = resolve(name, cid, system)
  source_ontology = re.search(r'Original GLB node \d+; ((?:FMA|UBERON):\d+);', part.get('provenance', {}).get('detail', ''))
  if source_ontology:
   source_id = source_ontology[1]
   if source_id.startswith('FMA:'):
    if len(by_fma[source_id.removeprefix('FMA:')]) == 1:
     match = record(terms[next(iter(by_fma[source_id.removeprefix('FMA:')]))])
    # The HRA node can be side-specific while TA98's FMA cross reference is
    # generic or points to the opposite side. Keep the node's own exact ID.
    match['fma'] = source_id
   else:
    match['ontology'] = source_id
  if match['ta98'] or match['fma'] or match.get('ontology'):
   concepts[cid] = match
   if match['ta98'] in terms:
    concept_paths[cid] = [code for code in ancestors[match['ta98']] if code in terms]

groups = ['Human body', 'Skin', 'Named skin regions', 'Body surface', 'Subcutaneous tissue',
 'Skeletal system', 'Bones', 'Cartilage', 'Joints', 'Ligaments', 'Tendons',
 'Muscular system', 'Muscles', 'Fascia', 'Digestive system', 'Respiratory system',
 'Circulatory system', 'Cardiovascular system', 'Heart', 'Arteries', 'Veins',
 'Urinary system', 'Integumentary system', 'Endocrine system', 'Exocrine system',
 'Lymphatic system', 'Nervous system', 'Central nervous system',
 'Peripheral nervous system', 'Sensory organs', 'Reproductive system',
 'Female genital system', 'Male genital system', 'Mammary gland',
 'Lacrimal glands', 'Salivary glands', 'Skin glands', 'Hair',
 'Head', 'Neck', 'Skull', 'Brain', 'Pelvis', 'Abdomen', 'Chest',
 'Upper arm', 'Forearm', 'Hand', 'Thigh', 'Leg', 'Foot']
groups += [
 'Arm', 'Back', 'Face', 'Upper limb', 'Lower limb',
 'Accessory digestive organs', 'Airways', 'Anal region', 'Antebrachial region',
 'Anterior arm', 'Anterior elbow', 'Anterior thigh', 'Anterior thorax',
 'Auricular region', 'Back & spine', 'Brachial region', 'Carpal region',
 'Concha of auricle', 'Cubital region', 'Digestive tract', 'Digits of foot',
 'Digits of hand', 'Endocrine glands', 'Epicranial regions', 'Epigastric region',
 'Facial regions', 'Femoral region', 'Foot regions', 'Frontal region',
 'Gluteal region', 'Hand regions', 'Infrahyoid region', 'Knee region',
 'Lateral neck', 'Leg regions', 'Longitudinal arch', 'Lungs', 'Lymph nodes',
 'Metatarsal region', 'Nasal region', 'Oral region', 'Orbital region',
 'Palm', 'Pectoral region', 'Perineal regions', 'Posterior ankle',
 'Posterior knee', 'Posterior leg', 'Salivary glands', 'Sole',
 'Sternocleidomastoid region', 'Suprahyoid region', 'Talocrural region',
 'Umbilical region',
]
groups = list(dict.fromkeys(groups + sorted(source_groups)))
group_records = {name:resolve(name) for name in groups}
group_records['All'] = record(terms['A01.0.00.000'])
group_records['Integumentary system'] = record(terms['A16.0.00.001'])
group_records['Skeletal system'] = record(terms['A02.0.00.000'])
group_records['Muscular system'] = record(terms['A04.0.00.000'])
group_records['Digestive system'] = record(terms['A05.0.00.000'])
group_records['Circulatory system'] = record(terms['A12.0.00.000'])
group_records['Lymphatic system'] = record(terms['A13.0.00.000'])
group_records['Endocrine system'] = record(terms['A11.0.00.000'])
group_records['Reproductive system'] = record(terms['A09.0.00.000'])
for group,ta98 in {
 'Muscles':'A04.0.00.000', 'Tendons':'A04.0.00.044',
 'Fascia':'A04.0.00.031', 'Sensory organs':'A15.0.00.000',
 'Lacrimal glands':'A15.2.07.057', 'Hair':'A16.0.00.014',
 'Head':'A01.1.00.001', 'Neck':'A01.1.00.012',
 'Skull':'A02.1.00.001', 'Pelvis':'A01.1.00.017',
 'Chest':'A01.1.00.014', 'Upper arm':'A01.1.00.022',
}.items():
 group_records[group] = record(terms[ta98])
# Wikipedia treats the exocrine glands as a system, while TA98 distributes
# them among several chapters. Keep the Latin descriptor but no invented ID.
group_records['Exocrine system']['latin'] = 'glandulae exocrinae'

output = ROOT / 'app/data/ta98-metadata.json'
output.parent.mkdir(parents=True, exist_ok=True)
parent_codes = {code for path in concept_paths.values() for code in path}
code_records = {code: {'name': terms[code]['name_en'], **record(terms[code])} for code in sorted(parent_codes)}
output.write_text(json.dumps({'source': SOURCE, 'byConcept':concepts,
                              'byConceptPath':concept_paths, 'byCode':code_records,
                              'byGroup':group_records}, separators=(',',':'), ensure_ascii=False) + '\n')
matched = sum(bool(value['ta98']) for value in concepts.values())
print(f'{output}: {matched} atlas concepts matched to a unique TA98 term; {len(concepts)} records include a TA98, FMA, or source ontology identifier')
