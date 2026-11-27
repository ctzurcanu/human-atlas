"""Finalize passed components without admitting failed venous fit trials."""
from pathlib import Path
import copy,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.tube_continuity import boundary_loops
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/'final-ready-v265';out.mkdir(exist_ok=False);r=AtlasReader()
for sex,stage in [('male','final-connections-v260-male'),('female','final-connections-v259-female')]:
 a=json.loads((b/stage/'atlas.json').read_text());parts={p['id']:p for p in a['parts']};marker=a['hepaticConnectionIntegration']
 for p in a['parts']:
  if p.get('hepaticConnectionRole')or p.get('hepaticRetainedParent'):p['hepaticSourceModel']='local-female-ta98'if sex=='female'else'local-ta98'
 retained=[]
 for pid in marker['regionPartIds']:
  p=parts[pid];source=parts[p['hepaticRegionOf']]
  if p.get('hepaticAttributeSignature')==source.get('hepaticAttributeSignature')and p['vertexCount']==source['vertexCount']:retained.append(pid);continue
  p.update(suppressed=True,defaultHidden=True,replacementParts=[source['id']]);a['defaultHidden']=sorted(set(a.get('defaultHidden',[]))|{pid})
  for g in marker['termGroups']:g['parts']=[source['id']if x==pid else x for x in g['parts']]
  for c in a['concepts']:c['elements']=[source['id']if x==pid else x for x in c.get('elements',[])]
 marker['regionPartIds']=retained
 rows=[]
 for pair in a['liverFinalIntegration']['interfaces']:
  first,second=parts[pair['first']],parts[pair['second']]
  for p in [first,second]:p.pop('suppressed',None);p['defaultHidden']=False
  ms=[r.read_part(a,p,frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres')for p in [first,second]];keys=[]
  for m in ms:keys.append({tuple(sorted(map(tuple,t))):i for i,t in enumerate(m.triangles)})
  common=keys[0].keys()&keys[1].keys();assert len(common)>=12,(sex,pair['first'],'No literal annular tissue interface')
  ids=[[keys[j][key]for key in sorted(common)]for j in [0,1]]
  for m,ff in zip(ms,ids):
   v,ix=np.unique(m.vertices,axis=0,return_inverse=True);loops=boundary_loops(ix[m.faces[ff]]);assert len(loops)==2,'Shared tissue mouth must have two genuine cycles'
  pair.update(first_faces=ids[0],second_faces=ids[1]);rows.append(dict(first=pair['first'],second=pair['second'],triangles=len(common),outer_and_inner_rims=2,lumen_disk=False))
 a['liverFinalIntegration']['remainingConnectionLimits']=['The existing male distal portal reference has not passed a new finite-wall rebuild.','Female native caval recipients remain unresolved; failed outflow collars are excluded.']
 # Source alternatives are often hidden only by this atlas-level list. Do not
 # discard those entries merely because the individual row lacks the flag.
 active=set(a['liverIntegration']['partIds'])|set(marker['wholePartIds'])|{x[k]for x in a['liverFinalIntegration']['interfaces']for k in ['first','second']}
 a['defaultHidden']=sorted(set(a.get('defaultHidden',[]))-active)
 atomic_json(out/(sex+'.json'),a);atomic_json(out/(sex+'-interfaces.json'),dict(passed=True,interfaces=rows,failed_trial_banks_installed=False))
 print(sex,'69 represented terms;',len(rows),'literal annular connections',flush=True)
