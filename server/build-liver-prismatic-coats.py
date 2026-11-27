"""Partition the existing finite capsule envelope into covered tissue coats."""
from pathlib import Path
from collections import defaultdict
import json,sys,time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader
from geometry_tools.tissue import pack,normal_words
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/'final-coats-v216';out.mkdir(exist_ok=False);start=time.monotonic()
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cat.read_text());p=next(p for p in a['parts']if p['id']=='TA98PART:ZA:Liver:capsule');m=AtlasReader().read_part(a,p,frame='model:local-ta98:assembled-metres',catalogue_path=cat)
count=p['tissueDrawGroups'][0]['count']//3;ov,oi=np.unique(m.triangles[:count].reshape(-1,3),axis=0,return_inverse=True);f=oi.reshape(-1,3)
inner=np.concatenate([m.triangles[g['start']//3:(g['start']+g['count'])//3]for g in p['tissueDrawGroups']if g.get('neighbor')]);iv,ii=np.unique(inner.reshape(-1,3),axis=0,return_inverse=True)
_,idx=min_weight_full_bipartite_matching(cKDTree(ov).sparse_distance_matrix(cKDTree(iv),.00015,output_type='coo_matrix').tocsr())
assert set(map(tuple,np.sort(idx[f],axis=1)))==set(map(tuple,np.sort(ii.reshape(-1,3),axis=1)))
iv=iv[idx];n=len(ov);fields=np.array([ov,.8*ov+.2*iv,.55*ov+.45*iv,iv]);v=np.r_[fields.reshape(-1,3),m.vertices]
bare=set(json.loads((b/'final-coats-v215/coats.json').read_text())['bareAreaFaces']);covered=np.array([i not in bare for i in range(count)])
surfaces={};cells={name:[]for name in ['serosa','subserosa','capsule']};external_neighbors={}
def add(sid,tri,first,second=None):
 surfaces[sid]=np.asarray(tri,int);cells[first].append(dict(surface=sid,side=1,neighbor=second))
 if second:cells[second].append(dict(surface=sid,side=-1,neighbor=first))
# Capsule's first `count` triangles retain the original source face order.
begin=0
for end in range(1,count+1):
 if end==count or covered[end]!=covered[begin]:
  sid='capsule:outer:'+str(begin)
  if covered[begin]:add(sid,f[begin:end]+2*n,'capsule','subserosa')
  else:add(sid,f[begin:end],'capsule')
  begin=end
add('serosa:external',f[covered],'serosa')
add('serosa:subserosa',f[covered]+n,'subserosa','serosa')
for j,g in enumerate(p['tissueDrawGroups'][1:]):
 if not g.get('neighbor'):continue # Rebuild only the source port thickness rims below.
 sid='capsule:inner:'+str(j);add(sid,m.faces[g['start']//3:(g['start']+g['count'])//3]+4*n,'capsule');external_neighbors[sid]=g['neighbor']
edges=defaultdict(list)
for i,t in enumerate(f):
 for x,y in zip(t,np.roll(t,-1)):edges[tuple(sorted((int(x),int(y))))].append((i,int(x),int(y)))
rim_groups=defaultdict(list)
for uses in edges.values():
 for k in range(3):
  labels=[(['serosa','subserosa','capsule'][k]if covered[u[0]]else'capsule')for u in uses]
  if len(uses)==2 and labels[0]==labels[1]:continue
  i,x,y=uses[0];first=labels[0];second=labels[1]if len(uses)==2 else None
  # Annular material edge, never a mouth-spanning disk.
  rim_groups[(k,first,second)].extend([[y+k*n,x+k*n,x+(k+1)*n],[y+k*n,x+(k+1)*n,y+(k+1)*n]])
for (k,first,second),tri in rim_groups.items():add('edge:'+str(k)+':'+first+':'+str(second),tri,first,second)
words=normal_words((v-[0,1.25,0]).astype('<f4').astype(float),surfaces)
effective=m.attributes['normals'].astype(np.int64)*m.attributes['normalSigns'][:,None]
for j,g in enumerate(p['tissueDrawGroups'][1:]):
 sid='capsule:inner:'+str(j)
 if sid in surfaces:words[sid]=effective[m.faces[g['start']//3:(g['start']+g['count'])//3]].astype('<i2')
ids={name:'TA98PART:ZA:Liver:'+name for name in cells}
result=pack(v,surfaces,cells,ids,m.frame,out/'packed',origin=[0,1.25,0],corner_words=words,max_pairs=500000,max_tests=8000000,checkpoint_directory=out/'checks')
for g,use in zip(result['parts']['capsule']['tissueDrawGroups'],cells['capsule']):
 if external_neighbors.get(use['surface']):g['neighbor']=external_neighbors[use['surface']]
atomic_json(out/'coats.json',dict(passed=result['passed'],parts=result['parts'],bank=result['bank'],bank_sha256=result['bank_sha256'],referenceSurfaceTriangleCount=count,bareAreaFaces=sorted(bare),envelope_literal=True,core_interfaces_literal=True,layer_fractions=[.2,.25,.55],thickness_status='Illustrative tissue allocation inside the existing finite envelope; not measured human histology',elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(passed=result['passed'],failures=result['checks']['failures'],elapsed_seconds=time.monotonic()-start)),flush=True)
if not result['passed']:raise ValueError('Prismatic coats failed; not installed')
