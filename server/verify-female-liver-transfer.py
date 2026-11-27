"""Replay the uniform transfer, literal attributes and occupied channel walls."""
from pathlib import Path
import itertools,json,sys
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'female-transfer-v206'
a=json.loads((out/'atlas.json').read_text());plan=json.loads((base/'female-transfer-v202/plan.json').read_text());fit=json.loads((base/'female-transfer-v202/fit-cleared-v204.json').read_text());m=json.loads(Path(plan['sourceCatalogue']).read_text());before=json.loads(Path(plan['targetCatalogue']).read_text());reader=AtlasReader();source={p['id']:p for p in m['parts']};parts={p['id']:p for p in a['parts']};meshes={};checks=[]
for pid in a['liverTransferIntegration']['partIds']:
 p=parts[pid];old=reader.read_part(m,source[p['liverTransferSourcePartId']],frame='source-male');new=reader.read_part(a,p,frame='female');meshes[pid]=new
 assert np.array_equal(old.faces,new.faces),'Source tissue ownership changed'
 assert set(old.attributes)==set(new.attributes)and all(np.array_equal(old.attributes[k],new.attributes[k])for k in old.attributes),'Source corner fields changed'
 expected=old.vertices*fit['scale']+fit['translation_m'];error=np.abs(new.vertices-expected);assert error.max()<1e-8
 checks.append(dict(part=pid,maximumStorageErrorMetres=float(error.max()),indicesAndAttributesLiteral=True))
oldparts={p['id']:p for p in before['parts']};retired=set(a['liverTransferIntegration']['retiredPartIds'])
assert all(p==parts[pid]for pid,p in oldparts.items()if pid not in retired),'Unrelated anatomy changed'
shared=[];alternatives=[];functional=set(a['liverCompletionIntegration']['functionalPartIds'])
def canonical(mesh,groups):
 f=np.concatenate([mesh.faces[g['start']//3:(g['start']+g['count'])//3]for g in groups]);tri=mesh.vertices[f]
 return {tuple(sorted(map(tuple,t)))for t in tri}
for pid in a['liverTransferIntegration']['partIds']:
 p=parts[pid]
 for neighbor in sorted({g['neighbor']for g in p.get('tissueDrawGroups',[])if g.get('neighbor')in meshes}):
  if pid>neighbor:continue
  q=parts[neighbor];g=[g for g in p['tissueDrawGroups']if g.get('neighbor')==neighbor];h=[g for g in q['tissueDrawGroups']if g.get('neighbor')==pid]
  if pid in functional and neighbor=='TA98PART:ZA:Liver:capsule':
   # Functional cells replace both gross lobes and their capsule in the
   # viewer. This inherited neighbour is an occlusion reference, not a
   # simultaneously displayed reciprocal tissue chart.
   alternatives.append(dict(part=pid,reference=neighbor,kind='alternative-presentation-occlusion'));continue
  assert h,'Missing reciprocal tissue ownership'
  first=canonical(meshes[pid],g);second=canonical(meshes[neighbor],h);assert first==second,'Shared tissue coordinates separated'
  shared.append(dict(first=pid,second=neighbor,faces=len(first),literalSharedCoordinates=True))
def solid(mesh):
 v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);s=md.Manifold(md.Mesh64(np.ascontiguousarray(v),np.ascontiguousarray(inv[mesh.faces],dtype=np.uint64)));assert s.status()==md.Error.NoError;return s
roles=['artery','portal','hepatic_vein','biliary'];solids={r:solid(meshes['TA98HEPATIC:female:'+r])for r in roles};pairs=[]
for first,second in itertools.combinations(roles,2):
 overlap=solids[first]^solids[second];assert overlap.status()==md.Error.NoError
 row=dict(first=first,second=second,volumeMetresCubed=float(overlap.volume()),triangles=int(overlap.num_tri()),passed=overlap.is_empty());pairs.append(row);assert row['passed'],'Distinct channel tissue overlaps'
atomic_json(out/'independent-transfer-checks.json',dict(passed=True,scale=fit['scale'],parts=checks,sharedInterfaces=shared,alternativeOcclusionReferences=alternatives,channelPairs=pairs,unrelatedAnatomyLiteral=True,sourceCatalogueSha256=plan['sourceCatalogueSha256'],candidateCatalogueSha256=digest_file(out/'atlas.json'),anatomicalAcceptance=False))
print(json.dumps(dict(passed=True,parts=len(checks),sharedInterfaces=len(shared),channelPairs=len(pairs))))
