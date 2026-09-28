#!/usr/bin/env python3
"""Build source-backed HGNC/CL DAGs. Downloads are cached; snapshots are hashed."""
import argparse, collections, concurrent.futures, csv, gzip, hashlib, json, re, subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
CACHE=ROOT/'work/biology'
SOURCES=ROOT/'scripts/biology-sources.json'
OBO='http://purl.obolibrary.org/obo/'

def normalized(name):
    return re.sub(r'\s*\((?:left|right)\)\s*$','',name.strip(),flags=re.I).lower().strip()
def ontology_id(value):
    match=re.search(r'(HGNC|CL|UBERON|FMA)[:_](\d+)',value or '',re.I)
    return match[1].upper()+':'+match[2] if match else ''
def read_csv(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig',newline=''),delimiter='\t' if path.suffix=='.tsv' else ','))
def ordered(values):return sorted(set(values))
def compact(data):return json.dumps(data,ensure_ascii=False,separators=(',',':'))+'\n'

def source_files(download):
    manifest=json.loads(SOURCES.read_text());CACHE.mkdir(parents=True,exist_ok=True)
    def fetch(item):
        name,entry=item;path=CACHE/entry['file']
        if not path.exists():
            if not download:raise RuntimeError(f'Missing {path}; run with --download.')
            subprocess.run(['curl','-fsSL','--retry','2','--max-time','300',entry['url'],'-o',str(path)],check=True)
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=entry['sha256']:raise RuntimeError(f'{name}: source hash changed. Preserve the pinned snapshot or review/update biology-sources.json.')
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:list(pool.map(fetch,manifest.items()))
    return manifest

def build(manifest):
    metadata=json.loads((ROOT/'app/data/ta98-metadata.json').read_text())
    catalogues=[]
    # Current anatomy geometry only; VR copies retain the same identities.
    for directory in ['public/models','.local-models']:
        for path in sorted((ROOT/directory).glob('*.json')):
            if path.name.endswith('.vr.json'):continue
            atlas=json.loads(path.read_text())
            if not isinstance(atlas.get('parts'),list) or not isinstance(atlas.get('concepts'),list) or atlas.get('scope')=='cell':continue
            available={p['id'] for p in atlas['parts'] if not p.get('suppressed')};parts_by_id={p['id']:p for p in atlas['parts']}
            names=collections.defaultdict(set);ids=collections.defaultdict(set)
            for concept in atlas['concepts']:
                if not available.intersection(concept['elements']):continue
                cid=concept['id'];members=[parts_by_id[pid] for pid in concept['elements'] if pid in available]
                # An organ concept can bundle a separately typed ligament material.
                # Cell/gene location leaves should select its organ geometry.
                organ=[p for p in members if p.get('material','').rsplit(':',1)[-1]=='Organ']
                targets={p['id'] for p in organ} if organ and len(organ)<len(members) else {cid}
                names[normalized(concept['name'])].update(targets)
                term=metadata['byConcept'].get(cid,{})
                # Parent-level terminology is not an exact anatomical identity.
                if metadata.get('byConceptMatch',{}).get(cid,{}).get('kind')!='parent':
                    for key in ['fma','ontology']:
                        if term.get(key):ids[ontology_id(term[key])].update(targets)
                if ontology_id(cid):ids[ontology_id(cid)].update(targets)
            catalogues.append({'path':str(path.relative_to(ROOT)),'atlas':atlas,'names':names,'ids':ids})
    human=json.loads((CACHE/'cl.json').read_text())['graphs'][0]
    full=json.loads((CACHE/'cl-full.json').read_text())['graphs'][0]
    full_nodes={n['id']:n for n in full['nodes']}
    anatomical_parents=collections.defaultdict(set)
    for edge in full['edges']:
        if edge['pred'] in ['is_a',OBO+'BFO_0000050'] and ontology_id(edge['sub']).startswith('UBERON:') and ontology_id(edge['obj']).startswith('UBERON:'):
            anatomical_parents[edge['sub']].add(edge['obj'])
    location_cache={}
    def anatomy_match(aid,names=(),index=0):
        catalogue=catalogues[index];found=set(catalogue['ids'].get(ontology_id(aid),()))
        term=full_nodes.get(aid) or full_nodes.get(OBO+ontology_id(aid).replace(':','_'),{})
        for xref in term.get('meta',{}).get('xrefs',[]):found.update(catalogue['ids'].get(ontology_id(xref['val']),()))
        for name in list(names)+[term.get('lbl','')]+[synonym['val'] for synonym in term.get('meta',{}).get('synonyms',[]) if synonym.get('pred','').endswith('hasExactSynonym')]:
            if name:found.update(catalogue['names'].get(normalized(name),()))
        return found
    def location_match(aid):
        if aid in location_cache:return location_cache[aid]
        found=set();coarse=False
        for index in range(len(catalogues)):
            frontier={aid};seen=set()
            for depth in range(12):
                matched=set().union(*(anatomy_match(item,index=index) for item in frontier)) if frontier else set()
                if matched:found.update(matched);coarse|=depth>0;break
                seen.update(frontier);frontier={parent for item in frontier for parent in anatomical_parents[item] if parent not in seen}
                if not frontier:break
        location_cache[aid]=(found,coarse);return found,coarse
    cl={ontology_id(n['id']):n for n in human['nodes'] if n['id'].startswith(OBO+'CL_') and n['type']=='CLASS' and not n.get('meta',{}).get('deprecated')}
    cells={cid:{'id':cid,'name':node['lbl'],'synonyms':ordered(x['val'] for x in node.get('meta',{}).get('synonyms',[]) if 0<len(x['val'])<=200),'source':node['id'],'description':node.get('meta',{}).get('definition',{}).get('val','')[:5000],'childIds':[]} for cid,node in cl.items()}
    location_parts=collections.defaultdict(set);location_evidence=collections.defaultdict(set);coarse_cells=set()
    subtype_edges=set();lineage_edges=set();cluster_edges=set()
    def cell_edge(child,predicate,parent):
        child,parent=ontology_id(child),ontology_id(parent)
        if child not in cl or parent not in cl or child==parent:return
        if predicate=='is_a':subtype_edges.add((parent,child))
        elif predicate in [OBO+'RO_0002202',OBO+'RO_0002207']:lineage_edges.add((parent,child))
        elif predicate in [OBO+'RO_0002203',OBO+'RO_0002210']:lineage_edges.add((child,parent))
        elif predicate==OBO+'RO_0015003':cluster_edges.add((parent,child))
    for edge in human['edges']:cell_edge(edge['sub'],edge['pred'],edge['obj'])
    location_predicates={OBO+'BFO_0000050',OBO+'BFO_0000066',OBO+'RO_0001025',OBO+'RO_0002100',OBO+'RO_0002130'}
    def add_location(cid,predicate,aid):
        if cid not in cells or predicate not in location_predicates or not ontology_id(aid).startswith('UBERON:'):return
        parts,coarse=location_match(aid)
        location_parts[cid].update(parts);location_evidence[cid].add(f'{full_nodes.get(aid,{}).get("lbl",ontology_id(aid))} ({predicate.rsplit("/",1)[-1]})')
        if coarse:coarse_cells.add(cid)
    for edge in full['edges']:add_location(ontology_id(edge['sub']),edge['pred'],edge['obj'])
    for axiom in full.get('logicalDefinitionAxioms',[]):
        for restriction in axiom.get('restrictions',[]):add_location(ontology_id(axiom['definedClassId']),restriction['propertyId'],restriction['fillerId'])
    genes={row['hgnc_id']:row for row in read_csv(CACHE/'hgnc.tsv') if row['status']=='Approved'}
    symbols={row['symbol']:cid for cid,row in genes.items()}
    aliases=collections.defaultdict(set)
    for cid,row in genes.items():
        for symbol in (row['alias_symbol']+'|'+row['prev_symbol']).split('|'):
            if symbol:aliases[symbol].add(cid)
    gene_parts=collections.defaultdict(set);gene_cells=collections.defaultdict(set);cell_genes=collections.defaultdict(set)
    gene_tables=collections.defaultdict(set);cell_tables=collections.defaultdict(set)
    unresolved_as={};unclassified=set();rows_used=0;conflicts=[];negative=[];table_stats={}
    for table,entry in manifest.items():
        if not table.startswith('hgnc-') and table not in ['cl-human','cl-full','ensembl']:
            rows=list(csv.reader((CACHE/entry['file']).open(encoding='utf-8-sig',newline='')))
            header_index=next((i for i,row in enumerate(rows) if 'AS/1' in row and any('CT/1'==x for x in row)),None)
            if header_index is None:continue
            header=[h.strip() for h in rows[header_index]]
            count=0
            for values in rows[header_index+1:]:
                row=dict(zip(header,values));path=[]
                for key in header:
                    if re.fullmatch(r'AS/\d+',key) and (row.get(key) or row.get(key+'/ID') or row.get(key+'/LABEL')):
                        path.append((row.get(key+'/ID',''),[row.get(key+'/LABEL',''),row.get(key,'')]))
                cids=[]
                for key in header:
                    if not re.fullmatch(r'CT/\d+',key):continue
                    cid=ontology_id(row.get(key+'/ID',''));label=row.get(key+'/LABEL','') or row.get(key,'')
                    if not cid.startswith('CL:') or not label:continue
                    if cid not in cells:
                        cells[cid]={'id':cid,'name':label+' · ASCT+B','source':entry['url'],'description':'Human ASCT+B cell type absent from the pinned human Cell Ontology view; source table classification only.','childIds':[]};unclassified.add(cid)
                    if cid in unclassified and cids:subtype_edges.add((cids[-1],cid))
                    cids.append(cid)
                # The deepest CT in a table row is the characterized cell type.
                if not cids:continue
                cid=cids[-1];matched=set()
                for index in range(len(catalogues)):
                    for aid,names in reversed(path):
                        found=anatomy_match(aid,names,index)
                        if found:matched.update(found);break
                if path and not matched:
                    aid,names=path[-1];unresolved_as[aid or names[0] or names[1]]=(names[0] or names[1],entry['url'])
                location_parts[cid].update(matched);cell_tables[cid].add(table)
                if path:location_evidence[cid].add(' > '.join((names[0] or names[1] or aid) for aid,names in path))
                markers=set()
                for key in header:
                    if not re.fullmatch(r'BGene/\d+',key):continue
                    text=row.get(key,'').strip();label=row.get(key+'/LABEL','').strip();gid=ontology_id(row.get(key+'/ID',''))
                    if re.search(r'negative|not express|absence|(?:^|\s)[A-Za-z0-9]+\s*[-−](?:\s|$)',text+' '+row.get(key+'/NOTES',''),re.I):negative.append((table,text));continue
                    known_symbol=symbols.get(text) or symbols.get(label)
                    if gid in genes:
                        if known_symbol and known_symbol!=gid:conflicts.append((table,cid,text,label,gid,known_symbol));continue
                    else:
                        gid=known_symbol or (next(iter(aliases[text])) if len(aliases[text])==1 else '')
                    if gid not in genes:continue
                    markers.add(gid)
                for gid in markers:
                    # Gene anatomy comes only from that marker's actual table row,
                    # never the union of all locations for its cell type.
                    gene_parts[gid].update(matched);gene_cells[gid].add(cid);cell_genes[cid].add(gid);gene_tables[gid].add(table)
                rows_used+=1;count+=1
            table_stats[table]=count
    for parent,child in subtype_edges:
        if parent in cells and child in cells:cells[parent]['childIds'].append(child)
    for cid,node in cells.items():
        if location_parts[cid]:node['matches']=ordered(location_parts[cid])
        if cell_genes[cid]:node['links']=[{'hierarchy':'genes','id':gid,'name':genes[gid]['symbol']} for gid in sorted(cell_genes[cid],key=lambda i:genes[i]['symbol'])]
        evidence=ordered(location_evidence[cid]);tables=ordered(cell_tables[cid])
        if evidence:node['description']+='\nAnatomical locations: '+'; '.join(evidence)[:2200]
        if cid in coarse_cells:node['description']+='\nSome locations use the nearest modeled anatomical ancestor; whole meshes represent tissue locations.'
        if tables:node['description']+='\nASCT+B: '+', '.join(tables)
        node['description']=node['description'][:5000]
        node['childIds']=sorted(set(node['childIds']),key=lambda i:cells[i]['name'].lower())
        if not node['childIds']:node.pop('childIds')
        if not node['description']:node.pop('description')
    type_children={child for parent,child in subtype_edges}
    type_roots=sorted(set(cells)-type_children,key=lambda i:cells[i]['name'].lower())
    views=[{'id':'types','name':'Types','roots':type_roots}]
    all_cells=list(cells.values())
    for prefix,title,edges in [('lineage','Lineage',lineage_edges),('cluster','Clusters',cluster_edges)]:
        included={cid for pair in edges for cid in pair};parents={child for parent,child in edges};children=collections.defaultdict(set)
        for parent,child in edges:children[parent].add(child)
        for cid in sorted(included):
            original=cells[cid];node={**original,'id':prefix+':'+cid,'childIds':[prefix+':'+child for child in sorted(children[cid],key=lambda i:cells[i]['name'].lower())]}
            node['links']=[{'hierarchy':'cell-types','id':cid,'name':original['name'][:200]},*original.get('links',[])]
            node['description']=f'{title}: source '+('develops_from / directly develops_from' if prefix=='lineage' else 'subcluster_of')+'.\n'+original.get('description','')
            node['description']=node['description'][:5000];all_cells.append(node)
        if included:views.append({'id':prefix,'name':title,'roots':[prefix+':'+cid for cid in sorted(included-parents,key=lambda i:cells[i]['name'].lower())]})
    family_rows=read_csv(CACHE/'family.csv');families={r['id']:{'id':'HGNC-GROUP:'+r['id'],'name':r['name'],'source':'https://www.genenames.org/data/genegroup/#!/group/'+r['id'],'childIds':[],**({'description':r['desc_comment'][:5000]} if r['desc_comment'] not in ['','NULL'] else {})} for r in family_rows}
    group_edges=[];child_groups=set();members=set()
    for row in read_csv(CACHE/'hierarchy.csv'):
        parent,child=row['parent_fam_id'],row['child_fam_id']
        if parent in families and child in families:families[parent]['childIds'].append(families[child]['id']);child_groups.add(child);group_edges.append((parent,child))
    for row in read_csv(CACHE/'gene_has_family.csv'):
        gid='HGNC:'+row['hgnc_id'];fid=row['family_id']
        if gid in genes and fid in families:families[fid]['childIds'].append(gid);members.add(gid)
    gene_nodes={}
    for gid,row in genes.items():
        node={'id':gid,'name':row['symbol']+' · '+row['name'],'source':'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/'+gid,'description':row['locus_group']+'; '+row['locus_type']+'; '+row['location']}
        node['synonyms']=ordered(x for x in (row['alias_symbol']+'|'+row['prev_symbol']+'|'+row['alias_name']).split('|') if 0<len(x)<=200)[:300]
        if gene_parts[gid]:node['matches']=ordered(gene_parts[gid])
        if gene_cells[gid]:node['links']=[{'hierarchy':'cell-types','id':cid,'name':cells[cid]['name'][:200]} for cid in sorted(gene_cells[gid],key=lambda i:cells[i]['name'].lower())]
        if gene_tables[gid]:node['description']+='\nGene biomarkers in ASCT+B: '+', '.join(ordered(gene_tables[gid]))+'. Anatomy is from the corresponding marker rows.'
        if row['alias_symbol']:node['description']+='\nAliases: '+row['alias_symbol'].replace('|',', ')
        gene_nodes[gid]=node
    transcript_count=protein_count=0;isoforms=[];extensions=collections.defaultdict(list)
    if 'ensembl' in manifest:
        ensembl_genes={row['ensembl_gene_id']:gid for gid,row in genes.items() if row['ensembl_gene_id'] and gene_parts[gid]}
        transcripts={};proteins={}
        entry=manifest['ensembl']
        if entry.get('format')=='gtf':
            with gzip.open(CACHE/entry['file'],'rt') as stream:
                for line in stream:
                    if line.startswith('#'):continue
                    columns=line.rstrip('\n').split('\t')
                    if len(columns)!=9 or columns[2] not in ['transcript','CDS']:continue
                    gene_id=re.search(r'gene_id "([^\"]+)"',columns[8]);gid=ensembl_genes.get(gene_id[1]) if gene_id else None
                    if not gid:continue
                    attrs=dict(re.findall(r'(\w+) "([^\"]*)"',columns[8]))
                    tid=attrs.get('transcript_id');version=attrs.get('transcript_version');tid=tid+('.'+version if version else '')
                    if columns[2]=='transcript':transcripts[tid]=(gid,attrs.get('transcript_name',tid),attrs.get('transcript_biotype','transcript'))
                    elif attrs.get('protein_id'):proteins[tid]=attrs['protein_id']+('.'+attrs['protein_version'] if attrs.get('protein_version') else '')
        else:
            for row in csv.reader((CACHE/entry['file']).open(),delimiter='\t'):
                if len(row)!=5:raise RuntimeError('Invalid Ensembl transcript table')
                eid,tid,pid,biotype,label=row;gid=ensembl_genes.get(eid)
                if not gid:continue
                transcripts[tid]=(gid,label or tid,biotype)
                if pid:proteins[tid]=pid
        for tid,(gid,label,biotype) in sorted(transcripts.items()):
            node={'id':tid,'name':label+' · '+tid,'source':'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t='+tid,'description':biotype+'. Gene-level anatomy inherited from '+genes[gid]['symbol']+'; no isoform-specific spatial evidence is asserted.','matchFrom':gid}
            if tid in proteins:
                pid=proteins[tid];node['childIds']=[pid];extensions[gid].append({'id':pid,'name':'Protein isoform · '+pid,'source':'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?t='+tid,'description':'Ensembl translation of '+tid+'. Anatomy is inherited at gene level.','matchFrom':gid});protein_count+=1
            extensions[gid].append(node);gene_nodes[gid].setdefault('childIds',[]).append(tid);transcript_count+=1
    ensembl_missing=[gid for gid in genes if gene_parts[gid] and gid not in extensions]
    refseq_count=0
    for gid in ensembl_missing:
        for accession in genes[gid]['refseq_accession'].split('|'):
            if not accession:continue
            tid='RefSeq:'+accession
            extensions[gid].append({'id':tid,'name':'RefSeq transcript · '+accession,'source':'https://www.ncbi.nlm.nih.gov/nuccore/'+accession,'description':'Transcript accession curated in the HGNC snapshot. This gene ID is absent from Ensembl 116. Anatomy is inherited at gene level.','matchFrom':gid})
            gene_nodes[gid].setdefault('childIds',[]).append(tid);refseq_count+=1
    for gid,extension_nodes in extensions.items():
        original=gene_nodes[gid];chunk={'schema':'human-atlas-hierarchy/v2','id':'gene-isoforms-'+gid.split(':')[1],'name':genes[gid]['symbol']+' transcripts and proteins','roots':[gid],'nodes':[{**original},*extension_nodes]}
        original['extension']='/assets/biology/genes/'+gid.split(':')[1]+'.json';original.pop('childIds',None)
        isoforms.append((original['extension'],chunk))
    names={node['id']:node['name'] for node in [*families.values(),*gene_nodes.values()]}
    for node in families.values():node['childIds']=sorted(set(node['childIds']),key=lambda i:names[i].lower())
    grouped={'id':'gene-families','name':'Gene families','childIds':[families[fid]['id'] for fid in sorted(set(families)-child_groups,key=lambda i:families[i]['name'].lower())]}
    ungrouped=[]
    for locus in sorted(set(row['locus_group'] for gid,row in genes.items() if gid not in members)):
        ungrouped.append({'id':'ungrouped:'+locus,'name':locus,'childIds':sorted((gid for gid,row in genes.items() if gid not in members and row['locus_group']==locus),key=lambda i:genes[i]['symbol'])})
    gene_roots=[grouped,{'id':'ungrouped','name':'Genes without an HGNC family','childIds':[n['id'] for n in ungrouped]}]
    gene_hierarchy={'schema':'human-atlas-hierarchy/v2','id':'genes','name':'Genes','source':'https://www.genenames.org/data/genegroup/','version':'HGNC 2026-09-25'+(' · Ensembl 116' if 'ensembl' in manifest else ''),'roots':[n['id'] for n in gene_roots],'nodes':[*gene_roots,*ungrouped,*families.values(),*gene_nodes.values()]}
    cell_hierarchy={'schema':'human-atlas-hierarchy/v2','id':'cell-types','name':'Cell Types','source':'https://obophenotype.github.io/cell-ontology/','version':'CL 2026-06-08 · HRA ASCT+B','roots':type_roots,'views':views,'nodes':all_cells}
    stats={'genes':len(genes),'families':len(families),'geneFamilyEdges':len(group_edges),'humanCellTypes':len(cl),'asctbAdditionalCellTypes':len(unclassified),'subtypeEdges':len(subtype_edges),'lineageEdges':len(lineage_edges),'clusterEdges':len(cluster_edges),'geneCellLinks':sum(map(len,gene_cells.values())),'genesWithAnatomy':sum(bool(parts) for parts in gene_parts.values()),'cellTypesWithAnatomy':sum(bool(parts) for parts in location_parts.values()),'genesWithIsoforms':len(extensions),'anatomyGenesWithoutEnsemblTranscripts':ensembl_missing,'refseqFallbackTranscripts':refseq_count,'transcripts':transcript_count+refseq_count,'proteinIsoforms':protein_count,'asctbRows':rows_used,'tables':table_stats,'unmatchedAnatomyLocations':len(unresolved_as),'skippedConflictingGeneIds':len(conflicts),'skippedNegativeMarkers':len(negative),'catalogues':[c['path'] for c in catalogues]}
    evidence={'sources':manifest,'statistics':stats,'conflictingGeneIds':conflicts,'negativeMarkers':negative,'unmatchedAnatomyLocations':[{'id':aid,'name':name,'source':url} for aid,(name,url) in sorted(unresolved_as.items())]}
    return gene_hierarchy,cell_hierarchy,evidence,isoforms

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--download',action='store_true');parser.add_argument('--check',action='store_true');args=parser.parse_args()
    genes,cells,evidence,extensions=build(source_files(args.download))
    for name,value in [('genes',genes),('cell-types',cells),('biology-provenance',evidence)]:
        path=ROOT/'public/assets'/f'{name}.json';text=compact(value)
        if args.check:
            if not path.exists() or path.read_text()!=text:raise RuntimeError(f'{name}.json is stale; rebuild.')
        else:path.write_text(text)
    for url,chunk in extensions:
        path=ROOT/'public'/url.lstrip('/');text=compact(chunk)
        if args.check:
            if not path.exists() or path.read_text()!=text:raise RuntimeError(f'{url} is stale; rebuild.')
        else:path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text)
    print(json.dumps(evidence['statistics'],indent=2))
if __name__=='__main__':main()
