#!/usr/bin/env python3
"""Build the Drugs DAG from pinned ChEMBL API snapshots and existing HGNC mappings."""
import argparse, collections, csv, hashlib, json, re, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'work/drugs'
MANIFEST=ROOT/'scripts/drug-sources.json'
API='https://www.ebi.ac.uk/chembl/api/data/'
FIELDS='molecule_chembl_id,pref_name,atc_classifications,molecule_synonyms,molecule_hierarchy,max_phase,veterinary,withdrawn_flag'
REQUESTS={
 'atc': ('atc_class','atc',{}),
 'classified': ('molecule','molecules',{'atc_classifications__isnull':'false','only':FIELDS}),
 'approved': ('molecule','molecules',{'max_phase':'4','only':FIELDS}),
 'mechanisms': ('mechanism','mechanisms',{}),
}
TOP_LEVELS=list('ABCDGHJLMNPRSV')
digest=lambda content:hashlib.sha256(content).hexdigest()
encoded=lambda value:(json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':'))+'\n').encode()

def read_url(url):
 for attempt in range(4):
  try:
   request=urllib.request.Request(url,headers={'User-Agent':'HumanAtlas-data-import/1.0','Accept':'application/json'})
   with urllib.request.urlopen(request,timeout=60) as response:return json.load(response)
  except Exception:
   if attempt==3:raise
   time.sleep(2**attempt)

def fetch_all(endpoint,key,params):
 url=API+endpoint+'.json?'+urllib.parse.urlencode({**params,'limit':1000,'offset':0})
 first=read_url(url);meta=first['page_meta'];total=meta['total_count'];limit=meta['limit']
 offsets=range(limit,total,limit)
 def page(offset):
  data=read_url(API+endpoint+'.json?'+urllib.parse.urlencode({**params,'limit':limit,'offset':offset}))
  assert data['page_meta']['total_count']==total,'Source changed during import'
  return data[key]
 rows=first[key]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for batch in pool.map(page,offsets):rows.extend(batch)
 assert len(rows)==total,f'{endpoint}: incomplete pagination'
 return {'request':url,'total_count':total,'records':rows}

def source_files(download):
 CACHE.mkdir(parents=True,exist_ok=True)
 manifest=json.loads(MANIFEST.read_text()) if MANIFEST.exists() else None
 missing=not manifest or any(not (CACHE/entry['file']).exists() for entry in manifest['sources'].values())
 if missing:
  if not download:raise RuntimeError('Source snapshots are missing. Run npm run build:drugs -- --download.')
  status=read_url(API+'status.json');print('Downloading '+status['chembl_db_version'],flush=True)
  def load(item):
   name,(endpoint,key,params)=item
   data=fetch_all(endpoint,key,params)
   print(f'{name}: {data["total_count"]:,} source rows',flush=True)
   return name,data
  with ThreadPoolExecutor(max_workers=2) as pool:data=dict(pool.map(load,REQUESTS.items()))
  molecules={row['molecule_chembl_id'] for name in ['classified','approved'] for row in data[name]['records']}
  target_ids=sorted({row['target_chembl_id'] for row in data['mechanisms']['records'] if row.get('target_chembl_id') and row['molecule_chembl_id'] in molecules})
  batches=[target_ids[index:index+75] for index in range(0,len(target_ids),75)]
  def targets(ids):return fetch_all('target','targets',{'target_chembl_id__in':','.join(ids)})
  with ThreadPoolExecutor(max_workers=3) as pool:target_batches=list(pool.map(targets,batches))
  target_rows=[row for batch in target_batches for row in batch['records']]
  assert {row['target_chembl_id'] for row in target_rows}==set(target_ids),'Missing drug targets'
  data['targets']={'request':API+'target.json','requests':[batch['request'] for batch in target_batches],'total_count':len(target_rows),'records':target_rows}
  final=read_url(API+'status.json');assert status['chembl_db_version']==final['chembl_db_version'],'ChEMBL release changed during import'
  data['status']=status
  fresh={'release':status['chembl_db_version'],'releaseDate':status['chembl_release_date'],'sources':{}}
  for name,value in data.items():
   content=encoded(value);entry={'file':name+'.json','url':value.get('request',API+'status.json'),'sha256':digest(content)}
   if name=='targets':entry['requests']=value['requests']
   fresh['sources'][name]=entry
   if manifest and entry['sha256']!=manifest['sources'][name]['sha256']:raise RuntimeError('Snapshot differs from the pinned source: '+name)
   (CACHE/entry['file']).write_bytes(content)
  if not manifest:MANIFEST.write_bytes(encoded(fresh))
  manifest=manifest or fresh
 result={}
 for name,entry in manifest['sources'].items():
  content=(CACHE/entry['file']).read_bytes();assert digest(content)==entry['sha256'],f'{name}: source hash mismatch'
  result[name]=json.loads(content)
 return manifest,result

def build(manifest,data):
 genes=json.loads((ROOT/'public/assets/genes.json').read_text())
 cells=json.loads((ROOT/'public/assets/cell-types.json').read_text())
 by_gene={node['id']:node for node in genes['nodes'] if re.fullmatch(r'HGNC:\d+',node['id'])}
 by_cell={node['id']:node for node in cells['nodes']}
 uniprot=collections.defaultdict(set)
 hgnc_snapshot=ROOT/'work/biology/hgnc.tsv'
 expected=json.loads((ROOT/'scripts/biology-sources.json').read_text())['hgnc-genes']['sha256']
 assert digest(hgnc_snapshot.read_bytes())==expected,'HGNC snapshot hash mismatch'
 with (ROOT/'work/biology/hgnc.tsv').open() as file:
  for row in csv.DictReader(file,delimiter='\t'):
   if row['hgnc_id'] not in by_gene:continue
   for accession in row['uniprot_ids'].split('|'):
    if accession:uniprot[accession].add(row['hgnc_id'])
 nodes={};evidence={'sources':manifest,'mechanisms':[],'unmatchedHumanComponents':[],'statistics':{}}
 def add(node):
  assert node['id'] not in nodes,'Duplicate node '+node['id']
  nodes[node['id']]=node
  return node
 def child(parent,id):
  children=parent.setdefault('childIds',[])
  if id not in children:children.append(id)
 atc_rows=data['atc']['records']
 for row in atc_rows:
  for level in range(1,6):
   code=row[f'level{level}'];id='ATC:'+code
   name=row['who_name'] if level==5 else row[f'level{level}_description']
   name=name[:1].upper()+name[1:].lower() if level<5 and name.isupper() else name
   if id not in nodes:add({'id':id,'name':name,'synonyms':[code],'childIds':[]})
   else:assert nodes[id]['name']==name,'Inconsistent ATC name: '+code
   if level>1:child(nodes['ATC:'+row[f'level{level-1}']],id)
 assert sorted(node['id'][4:] for node in nodes.values() if len(node['id'])==5)==TOP_LEVELS,'All 14 main ATC groups are required'
 molecules={row['molecule_chembl_id']:row for row in data['approved']['records']}
 molecules.update({row['molecule_chembl_id']:row for row in data['classified']['records']})
 # Preserve source drug identities, including combinations and multiple ATC assignments.
 unknown_atc=collections.defaultdict(set)
 for molecule in molecules.values():
  for code in molecule['atc_classifications']:
   if 'ATC:'+code not in nodes:unknown_atc[code].add(molecule['molecule_chembl_id'])
 # Some molecule records retain ATC codes absent from the release's category table.
 # Keep their original code under its deepest known parent without inventing a label.
 for code,ids in sorted(unknown_atc.items()):
  names=sorted({molecules[id]['pref_name'] for id in ids if molecules[id]['pref_name']})
  node=add({'id':'ATC:'+code,'name':names[0].lower() if len(names)==1 else code,'synonyms':[code],'childIds':[]})
  parent=next(('ATC:'+code[:length] for length in [5,4,3,1] if 'ATC:'+code[:length] in nodes),None)
  assert parent,'Unknown ATC top group: '+code
  child(nodes[parent],node['id'])
 unclassified=add({'id':'UNCLASSIFIED','name':'Unclassified','childIds':[]})
 for id,molecule in sorted(molecules.items()):
  name=molecule['pref_name'] or next((item['molecule_synonym'] for item in molecule['molecule_synonyms'] if item.get('molecule_synonym')),id)
  aliases=sorted({item['molecule_synonym'] for item in molecule['molecule_synonyms'] if item.get('molecule_synonym') and len(item['molecule_synonym'])<=200})
  node=add({'id':id,'name':name[:500],'synonyms':aliases[:300],'childIds':[]})
  codes=molecule['atc_classifications']
  if codes:
   for code in codes:child(nodes['ATC:'+code],id)
  else:child(unclassified,id)
 targets={row['target_chembl_id']:row for row in data['targets']['records']}
 gene_cache={}
 def target_genes(target):
  if target['target_chembl_id'] in gene_cache:return gene_cache[target['target_chembl_id']]
  ids=set()
  if target.get('tax_id')==9606 and target.get('organism')=='Homo sapiens':
   for component in target['target_components']:
    explicit={xref['xref_id'] for xref in component.get('target_component_xrefs',[]) if xref['xref_src_db']=='HGNC' and xref['xref_id'] in by_gene}
    accession=component.get('accession');fallback=uniprot.get(accession,set())
    mapped=explicit or (fallback if len(fallback)==1 else set());ids.update(mapped)
    if not mapped:evidence['unmatchedHumanComponents'].append({'target':target['target_chembl_id'],'accession':accession,'name':component['component_description']})
  gene_cache[target['target_chembl_id']]=sorted(ids)
  return sorted(ids)
 for record in sorted(data['mechanisms']['records'],key=lambda row:row['mec_id']):
  drug=record['molecule_chembl_id']
  if drug not in molecules:continue
  target=targets.get(record.get('target_chembl_id'));gids=target_genes(target) if target and record.get('molecular_mechanism')==1 else []
  name=record['mechanism_of_action'] or (target['pref_name'] if target else 'Mechanism')
  node=add({'id':'MECHANISM:'+str(record['mec_id']),'name':name[:500],'childIds':gids})
  child(nodes[drug],node['id'])
  target_summary={key:target[key] for key in ['target_chembl_id','pref_name','organism','tax_id','target_type']} if target else None
  evidence['mechanisms'].append({**record,'humanGenes':gids,'target':target_summary})
  for gid in gids:
   if gid in nodes:continue
   original=by_gene[gid]
   links=[{'hierarchy':'genes','id':gid,'name':original['name'].split(' · ')[0]}]
   links.extend(link for link in original.get('links',[]) if link['hierarchy']=='cell-types' and link['id'] in by_cell)
   add({'id':gid,'name':original['name'],'matches':original.get('matches',[]),'links':links})
 # Empty categories and unknown mechanisms stay in the tree without invented anatomy.
 for node in nodes.values():
  if 'childIds' in node:node['childIds'].sort(key=lambda id:(nodes[id]['name'].casefold(),id))
  if not node.get('childIds'):node.pop('childIds',None)
  if not node.get('matches'):node.pop('matches',None)
 stats=evidence['statistics']
 stats.update({'atcSourceSubstances':len(atc_rows),'atcLevels':{str(level):sum(re.fullmatch(pattern,node['id'][4:]) is not None for node in nodes.values() if node['id'].startswith('ATC:')) for level,pattern in enumerate([r'[A-Z]',r'[A-Z]\d{2}',r'[A-Z]\d{2}[A-Z]',r'[A-Z]\d{2}[A-Z]{2}',r'[A-Z]\d{2}[A-Z]{2}\d{2}'],1)},'drugs':len(molecules),'unclassifiedDrugs':len(unclassified.get('childIds',[])),'mechanisms':len(evidence['mechanisms']),'humanTargetGenes':len([id for id in nodes if id.startswith('HGNC:')]),'mappedTargetGenes':sum(bool(node.get('matches')) for node in nodes.values() if node['id'].startswith('HGNC:')),'atcCodesAbsentFromCategoryTable':dict(sorted((code,sorted(ids)) for code,ids in unknown_atc.items()))})
 hierarchy={'schema':'human-atlas-hierarchy/v2','id':'drugs','name':'Drugs','version':manifest['release'],'roots':['ATC:'+code for code in TOP_LEVELS]+['UNCLASSIFIED'],'nodes':list(nodes.values())}
 evidence['biologyInputs']={path:digest((ROOT/path).read_bytes()) for path in ['public/assets/genes.json','public/assets/cell-types.json','work/biology/hgnc.tsv']}
 return hierarchy,evidence

def report(hierarchy,evidence):
 stats=evidence['statistics'];nodes={node['id']:node for node in hierarchy['nodes']}
 view=lambda id:'http://localhost:3016/?'+urllib.parse.urlencode({'model':'male-detail','tree':'guest:drugs','bio':id,'context':'0.18'})
 lines=['# Drugs hierarchy','','Drugs appears before Physical Exercise and Chakras in the Additional hierarchies menu. It loads only when selected.','','## Imported snapshot','',f"Source: {evidence['sources']['release']}, released {evidence['sources']['releaseDate']}. All ATC entries supplied by that ChEMBL release are retained; this is a pinned classification snapshot, not a live WHO index.",'','| Data | Count |','| --- | ---: |']
 for level,count in stats['atcLevels'].items():lines.append(f'| ATC level {level} | {count:,} |')
 for key,label in [('drugs','Drug records'),('unclassifiedDrugs','Drug records without ATC codes'),('mechanisms','Curated mechanism records'),('humanTargetGenes','Human target genes'),('mappedTargetGenes','Target genes with anatomy mappings')]:lines.append(f'| {label} | {stats[key]:,} |')
 lines+=['','## Categories','','| ATC group | Open in viewer |','| --- | --- |']
 for id in hierarchy['roots']:lines.append(f"| {id.removeprefix('ATC:')} | [{nodes[id]['name']}]({view(id)}) |")
 lines+=['','## Mapping rules','', '- The five ATC levels lead to the associated ChEMBL drug records, their curated mechanisms, human target genes, and available anatomical items or sets. Multiple ATC assignments and shared targets retain all their parent branches.','- Drug records include substances with ATC assignments and records classified as phase 4 in ChEMBL. Historical or withdrawn records remain present; this is not a list of currently authorized products. Salts and combination records retain their source identities.','- A gene is linked only for a molecular mechanism with an explicitly human target. Components use source HGNC cross references, or an unambiguous UniProt-to-HGNC match from the pinned HGNC snapshot. Nonhuman pathogen targets are not assigned to human homologues.','- Anatomy comes from the existing gene-location evidence. These are mapped target-gene locations; they do not establish drug action, measured target expression, or clinical effects throughout an entire organ mesh. Cell Types links refer to existing marker associations.','- Unmapped categories, drugs and mechanisms remain searchable without anatomy checkboxes. Missing mappings are not filled with whole systems or guessed organs.','- Share URLs, slides, connection snapshots and embeds use the same built-in hierarchy ID, `guest:drugs`. Search accepts drug names and synonyms, ChEMBL IDs, ATC codes, target names and HGNC IDs.','', '## Sources','','All source credits are in [ATTRIBUTION.md](ATTRIBUTION.md#drugs). Exact request URLs and source hashes are in [scripts/drug-sources.json](scripts/drug-sources.json). Mechanism evidence, original references, mapping gaps and input hashes are in [public/assets/drug-provenance.json](public/assets/drug-provenance.json).','','## Rebuild and validate','','```sh','npm run build:drugs -- --download  # fetch missing pinned source snapshots','npm run build:drugs               # rebuild from the hash-verified cache','npm run build:drugs -- --check','npm run test:drugs','npm run check','npm run build','```','','The ignored `work/drugs` directory stores source snapshots. UniProt matching also uses the pinned HGNC source in `work/biology`, downloaded by `npm run build:biology -- --download`. Deployable data is in `public/assets/drugs.json`.','']
 return '\n'.join(lines)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--download',action='store_true');parser.add_argument('--check',action='store_true');args=parser.parse_args()
 manifest,data=source_files(args.download);hierarchy,evidence=build(manifest,data)
 for filename,value in [('drugs',hierarchy),('drug-provenance',evidence)]:
  path=ROOT/'public/assets'/f'{filename}.json';content=encoded(value)
  if args.check:assert path.read_bytes()==content,f'{filename} needs rebuilding'
  else:path.write_bytes(content)
 path=ROOT/'DRUGS-HIERARCHY.md';content=report(hierarchy,evidence)
 if args.check:assert path.read_text()==content,'Drugs report needs rebuilding'
 else:path.write_text(content)
 print(json.dumps(evidence['statistics'],indent=2))

if __name__=='__main__':main()
