import {readFileSync,writeFileSync} from 'node:fs';
const evidence=JSON.parse(readFileSync('public/assets/biology-provenance.json','utf8')),stats=evidence.statistics;
const genes=JSON.parse(readFileSync('public/assets/genes.json','utf8')),byGene=new Map(genes.nodes.map(node=>[node.id,node]));
const label=text=>text.replace(/[|\[\]]/g,' '),viewer=(hierarchy,id='',view='')=>{const url=new URL('http://localhost:3016/');url.searchParams.set('model','male-detail');url.searchParams.set('tree',`guest:${hierarchy}`);if(id)url.searchParams.set('bio',id);if(view)url.searchParams.set('bioView',view);url.searchParams.set('context','.18');return url.href;};
const geneLink=id=>`[${label(byGene.get(id)?.name??id)}](${viewer('genes',id)})`;
const unmatched=evidence.unmatchedAnatomyLocations.map(item=>`| [${label(item.name)}](https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?iri=${encodeURIComponent('http://purl.obolibrary.org/obo/'+item.id.replace(':','_'))}) | ${item.id} | [ASCT+B table](${item.source}) |`);
const conflicts=evidence.conflictingGeneIds.map(([table,cid,symbol,name,gid,expected])=>`| [${table}](${evidence.sources[table].url}) | [${cid}](${viewer('cell-types',cid)}) | ${label(symbol||name)} | ${geneLink(gid)} | ${geneLink(expected)} |`);
const text=`# Genes and Cell Types guest hierarchies

The Additional hierarchies menu (\`v\`) order is **Genes → Cell Types → Physiology → Dermatomes and Myotomes → Drugs → Physical Exercise → Chakras**. Open [Genes](${viewer('genes')}) or [Cell Types](${viewer('cell-types')}). Both provide search, selectable anatomy leaves, visibility controls, and links to the other hierarchy.

## Imported depth and counts

| Data | Count |
| --- | ---: |
| Approved HGNC genes | ${stats.genes.toLocaleString()} |
| HGNC groups/families, retaining every parent relationship | ${stats.families.toLocaleString()} |
| Group-to-group edges | ${stats.geneFamilyEdges.toLocaleString()} |
| Human Cell Ontology classes | ${stats.humanCellTypes.toLocaleString()} |
| Additional human ASCT+B types absent from that CL view | ${stats.asctbAdditionalCellTypes.toLocaleString()} |
| Type edges | ${stats.subtypeEdges.toLocaleString()} |
| Developmental lineage edges | ${stats.lineageEdges.toLocaleString()} |
| Curated gene–cell marker associations | ${stats.geneCellLinks.toLocaleString()} |
| Genes with modeled anatomy in at least one catalogue | ${stats.genesWithAnatomy.toLocaleString()} |
| Cell types with directly mapped anatomy in at least one catalogue | ${stats.cellTypesWithAnatomy.toLocaleString()} |
| Transcripts for mapped genes | ${stats.transcripts.toLocaleString()} |
| Ensembl protein translations/isoforms | ${stats.proteinIsoforms.toLocaleString()} |

Genes retain the complete HGNC family DAG through its individual gene nodes (up to 10 levels in this snapshot). Ungrouped approved genes remain under their HGNC locus group. Mapped genes have transcript and translated protein branches, loaded per gene; these add two more levels before anatomical leaves. Cell Types retain the human CL subtype DAG (up to 16 levels). Additional ASCT+B entries are visibly marked **ASCT+B** and retain their table parentage without presenting it as a CL assertion.

[Types](${viewer('cell-types','','types')}) uses **is_a**. [Lineage](${viewer('cell-types','','lineage')}) uses **develops_from / directly develops_from**, directed from progenitor to derivative. This is cell development, not a species phylogeny. Experiment-specific cluster individuals imported into CL from other taxa are excluded from the human type hierarchy.

Examples: [ALB / albumin](${viewer('genes','HGNC:399')}), [hepatocyte](${viewer('cell-types','CL:0000182')}), [periportal hepatocyte](${viewer('cell-types','CL:0019026')}), [HOX families](${viewer('genes','homeobox')}).

## Anatomy mapping rules

1. Gene anatomy comes from the anatomical path in that gene's actual **ASCT+B biomarker row**. It does not inherit all other sites associated with the same cell type. These are marker/location associations, not claims that a gene causes the structure to develop.
2. Cell anatomy comes from ASCT+B rows and CL **part_of**, **located_in**, soma-location, and synaptic-location axioms.
3. Resolve exact FMA/UBERON identifiers, exact ontology synonyms, and exact atlas names in each available catalogue. Parent-level TA98 assignments are excluded from exact ID matching. In an ASCT+B path, prefer the deepest modeled structure; for a CL microscopic location, use its nearest modeled UBERON ancestor where available.
4. Anatomical leaves represent the existing whole meshes. They do not add cellular geometry, resolve tissue distributions inside an organ, or claim measured expression throughout every vertex. When a source organ concept bundles a separately typed ligament material, the organ mapping selects its organ geometry.
5. Unmapped source nodes remain visible without anatomy checkboxes. Negative markers and conflicting gene symbol/ID pairs are excluded from positive associations. A shared mesh retains its original ID; Explode assigns it once, even with multiple biological parents.
6. Transcript/protein anatomy is **inherited at gene level**, explicitly labeled in each branch. There is no isoform-specific spatial assertion. Ensembl 116 contains all imported mapped-gene transcripts except ${stats.anatomyGenesWithoutEnsemblTranscripts.map(geneLink).join(' and ')}; their ${stats.refseqFallbackTranscripts??0} HGNC-listed RefSeq transcript accessions are retained as a labeled fallback.

## Source snapshots and attribution

Credits and exact source snapshots are collected in [ATTRIBUTION.md](ATTRIBUTION.md#additional-hierarchies). Input hashes are in [scripts/biology-sources.json](scripts/biology-sources.json); mapping evidence is in [public/assets/biology-provenance.json](public/assets/biology-provenance.json).

## Locations without a modeled match

The following ${stats.unmatchedAnatomyLocations} table locations could not be resolved to available modeled anatomy. Their cell/gene source records remain browseable. Links open the ontology entry and source table.

| Location | ID | Evidence |
| --- | --- | --- |
${unmatched.join('\n')}

## Source gene ID conflicts

${stats.skippedConflictingGeneIds} biomarker entries have a recognized symbol that disagrees with the supplied HGNC ID. They are skipped rather than guessing which column is correct. Both gene candidates and the cell type link to the viewer for review. ${stats.skippedNegativeMarkers} explicitly negative markers are also excluded.

<details>
<summary>Review the conflicting source entries</summary>

| Table | Cell type | Supplied symbol | ID points to | Recognized symbol points to |
| --- | --- | --- | --- | --- |
${conflicts.join('\n')}

</details>

## Build and verification

`;
const commands='```sh\nnpm run build:biology              # use cached, hash-verified source snapshots\nnpm run build:biology -- --download # download missing snapshots, then verify hashes\nnpm run build:biology -- --check\nnpm run test:biology\nnode scripts/report-biological-hierarchies.mjs\nnpm run check\nnpm run build\n```\n\nThe source cache lives in ignored `work/biology`; deployable JSON and per-gene chunks live under `public/assets`. The build copies them into `dist/assets`. Genes and Cell Types load only when selected; transcript/protein chunks load when requested. Shared graph resolution and anatomy arrays are cached, and long child lists are paged. Query, relation view, and requested transcript branches persist in share URLs, slides and connection snapshots.\n\nThe validation checks all eight available anatomy catalogues, every source/cross-hierarchy reference, cycles and dangling references, per-gene chunks, hierarchy navigation, URL restoration, and finite Explode positions with every mesh assigned exactly once. Browser checks cover menu order, marker links, selectable leaves, transcript loading, and lineage controls. No Quest headset was available for a hardware test.\n';
writeFileSync('BIOLOGICAL-HIERARCHIES.md',text+commands);
