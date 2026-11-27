"""Read-only native brain geometry comparison; not anatomy/shape acceptance."""
import collections,hashlib,json,re,sqlite3,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.ta98_authority import corrected_terms
began=time.monotonic();work=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';output=work/'brain-source-comparison-v1/native-brain.json'
if output.exists():raise SystemExit('Refusing to overwrite native brain evidence')
source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';authority=ROOT/'.local-models/ta98-authority.sqlite'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source_sha=sha(source);assert source_sha=='46262d8ab07682fcae3ff0ee634f652abad766c1e010f6e2ec311755aea9abf5'
atlas=json.loads(source.read_bytes());concepts={c['id']:c for c in atlas['concepts']};hidden=set(atlas.get('defaultHidden',[]))
db=sqlite3.connect('file:'+str(authority)+'?mode=ro',uri=True);db.row_factory=sqlite3.Row
terms=corrected_terms({r['id']:dict(r) for r in db.execute('select * from ta98')});children=collections.defaultdict(list)
for code,t in terms.items():children[t['parent_id']].append(code)
scope={'A14.1.03.001'};pending=list(scope)
while pending:
 for c in children[pending.pop()]:
  if c not in scope:scope.add(c);pending.append(c)
def normal(s):
 s=re.sub(r'\.(?:l|r)$','',s.lower());s=re.sub(r'\b(?:left|right)\b','',s)
 return re.sub(r'\s+',' ',re.sub(r'[^a-z0-9]+',' ',s)).strip()
brain_labels={normal(terms[c]['name_en']) for c in scope}
for r in db.execute('select * from synonyms'):
 if r['id'] in scope and r['lang'].lower() in {'en','eng','english'}:brain_labels.add(normal(r['synonym']))
brain_groups={'brain','cerebrum','cerebellum','brainstem','telencephalon','diencephalon','midbrain','pons','medulla oblongata'}
def explicit_brain(p):
 return any(g.lower() in brain_groups for g in p.get('groups',[])) or ':Brain' in p['id'] or p.get('ta98Term') in scope or concepts.get(p['conceptId'],{}).get('ta98Term') in scope
paired_names={normal(p.get('sourceId',p['name'])) for p in atlas['parts'] if explicit_brain(p)}
def candidate(p):
 return explicit_brain(p) or (p.get('system')=='nervous' and (normal(p.get('sourceId',p['name'])) in paired_names or normal(p['name']) in brain_labels))
def category(p):
 name=p['name'].lower();groups={g.lower() for g in p.get('groups',[])};label=normal(p.get('sourceId',p['name']))
 if p.get('system') in ['arterial','venous']:return 'blood-vessel-excluded'
 if re.search(r'\b(?:falx|tentorium|dura|arachnoid|pia|meninges)\b',name):return 'covering-excluded'
 if 'nerve' in name and not re.search(r'\b(?:nucleus|nuclei)\b',name):return 'cranial-nerve-excluded'
 if re.search(r'\b(?:ventricle|aqueduct|central canal)\b',name):return 'cavity-boundary-excluded'
 if 'choroid plexus' in name:return 'choroid-plexus-separate'
 if 'septum pellucidum' in name:return 'septal-partition-separate'
 if re.search(r'\b(?:sulcus|sulci|fissure|fossa)\b',name) and not re.search(r'\b(?:gyrus|gyri|lobule|cortex)\b',name):return 'cortical-or-brainstem-boundary-patch-excluded-from-neural-tissue'
 if p.get('system')=='endocrine' or re.search(r'\b(?:pituitary|pineal|hypophysis)\b',name):return 'endocrine-adjacent-separate'
 if 'cerebellum' in groups or any(x in label for x in ['cerebell','vermis','quadrangular lobule','semilunar lobule','biventral lobule','gracile lobule','central lobule','flocculus','culmen','declive']):return 'cerebellar-neural-tissue'
 if 'brainstem' in groups or any(x in label for x in ['midbrain','pons','medulla oblongata','colliculus','nucleus of','cochlear nucleus','vestibular nuclei','salivatory nucleus','red nucleus','substantia nigra','gracile tubercle','base of peduncle','peduncle of midbrain']):return 'brainstem-neural-tissue'
 if 'neo-cortex' in groups or any(x in label for x in ['gyrus','gyri','cuneus','parietal lobule','temporal pole','temporal plane','insula']):return 'cortical-neural-surface-region'
 return 'deep-cerebral-diencephalic-white-matter-neural-tissue'
reader=AtlasReader(ROOT);rows=[]
for p in atlas['parts']:
 if not candidate(p):continue
 active=bool(p['indexCount'] and not(p.get('suppressed') or p.get('defaultHidden') or p['id'] in hidden))
 row=dict(id=p['id'],name=p['name'],sourceId=p.get('sourceId'),category=category(p),active=active,system=p.get('system'),groups=p.get('groups',[]),ta98Term=p.get('ta98Term') or concepts.get(p['conceptId'],{}).get('ta98Term'),triangles=p['indexCount']//3,vertices=p['vertexCount'],sourceLabel=p.get('provenance',{}).get('label'),reviewStatus=p.get('reviewStatus'))
 if active:
  mesh=reader.read_part(atlas,p,frame='active-male-metres',catalogue_path=source,catalogue_sha=source_sha)
  row.update(geometrySha256=mesh.geometry_sha,orientedSurfaceSha256=mesh.surface_sha,
   storageView=dict(chunkURL=atlas['chunks'][p['chunk']]['url'],positions=p['positions'],indices=p['indices'],vertexCount=p['vertexCount'],indexCount=p['indexCount']),
   triangleTopologySha256=hashlib.sha256(mesh.faces.astype('<u4').tobytes()).hexdigest())
 rows.append(row)
neural_categories={'brainstem-neural-tissue','cerebellar-neural-tissue','cortical-neural-surface-region','deep-cerebral-diencephalic-white-matter-neural-tissue'}
neural=[r for r in rows if r['active'] and r['category'] in neural_categories]
def metrics(rs):
 unique={r['orientedSurfaceSha256']:r for r in rs};sources={r.get('sourceId') or r['id'] for r in rs}
 return dict(placedPartFragments=len(rs),genuineNamedSourceRegionsWithLaterality=len(sources),unsidedSourceRegionNames=len({normal(s) for s in sources}),
  uniqueExactWorldSurfaceMeshes=len(unique),uniqueIndexedWorldGeometryViews=len({r['geometrySha256'] for r in rs}),
  uniqueSurfaceTriangleTotal=sum(r['triangles'] for r in unique.values()),placedFragmentTriangleTotal=sum(r['triangles'] for r in rs),
  exactDuplicateSurfaceInstances=len(rs)-len(unique),geometryIsDistinctFromConceptLabels=True)
categories={k:metrics([r for r in rows if r['active'] and r['category']==k]) for k in sorted({r['category'] for r in rows})}
result=dict(schema='native-active-male-brain-source-audit/v1',status='complete-read-only-native-geometry-census',inputs={str(source):source_sha,str(authority):sha(authority)},
 authorityBrainRoot='A14.1.03.001',authorityBrainDescendants=len(scope),nativeBrainNeuralTissue=metrics(neural),categories=categories,
 selectedCandidateParts=len(rows),activeCandidateParts=sum(r['active'] for r in rows),
 selectionPolicy='Existing explicit Brain/Cerebrum/Cerebellum/Brainstem source groups, their exact laterality-paired source identities, or exact authority brain-descendant labels/metadata. Every selected active mesh is replayed by AtlasReader. Material fragments sharing a sourceId are one named region; exact world surfaces are deduplicated independently.',
 scopeQualification='Counts distinguish physical source geometry from concepts and from sulci, cavities, choroid plexus, septal partition, coverings and cranial nerves. Tissue means named neural source geometry, not verified closed parenchymal solids. Native mirrored left/right packed views are distinct anatomical placements; original pre-bake asset sharing is not inferable from packed catalogue.',
 surfaceQuality='No render, cross-section or self-intersection inspection was performed; mesh density and named subdivision are resolution evidence only. BP3D review-pending fields remain explicit.',
 rows=rows,sourceOrGeometryOrDefaultChanged=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
assert sha(source)==source_sha
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(output=str(output),neural=result['nativeBrainNeuralTissue'],categories=categories,elapsedSeconds=result['elapsedSeconds'])))
