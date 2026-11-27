"""Finite serosa/subserosa on the covered source liver, preserving bare area."""
from pathlib import Path
from collections import defaultdict
import argparse,importlib.util,json,sys,time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader
from geometry_tools.tissue import pack,normal_words
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim
ap=argparse.ArgumentParser();ap.add_argument('--regularize-mask',action='store_true');args=ap.parse_args()
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/('final-coats-v215' if args.regularize_mask else 'final-coats-v213');out.mkdir(exist_ok=False);start=time.monotonic()
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cat.read_text());p=next(p for p in a['parts']if p['id']=='TA98PART:ZA:Liver:capsule');m=AtlasReader().read_part(a,p,frame='model:local-ta98:assembled-metres',catalogue_path=cat)
count=p['tissueDrawGroups'][0]['count']//3;tri=m.triangles[:count];ov,oi=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);f=oi.reshape(-1,3)
inner=np.concatenate([m.triangles[g['start']//3:(g['start']+g['count'])//3]for g in p['tissueDrawGroups']if g.get('neighbor')]);iv,ii=np.unique(inner.reshape(-1,3),axis=0,return_inverse=True)
matrix=cKDTree(ov).sparse_distance_matrix(cKDTree(iv),.00015,output_type='coo_matrix');_,idx=min_weight_full_bipartite_matching(matrix.tocsr())
assert set(map(tuple,np.sort(idx[f],axis=1)))==set(map(tuple,np.sort(ii.reshape(-1,3),axis=1)))
direction=ov-iv[idx];direction/=np.linalg.norm(direction,axis=1)[:,None]
fields=np.array([ov+.00005*direction,ov+.00004*direction,ov]);n=len(ov);v=np.r_[fields.reshape(-1,3),m.vertices]
bare=set(next(g for g in a['liverFeatureIntegration']['termGroups']if g['term']=='A05.8.01.008')['surface']['faces']);covered=np.array([i not in bare for i in range(count)]);assert bare and covered.any()
mask_check=None
if args.regularize_mask:
 spec=importlib.util.spec_from_file_location('region_labels',ROOT/'server/hepatic-region-labels.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 # Close only the three actual source ports in this label-topology helper.
 # Their triangles are never emitted as anatomy.
 from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank
 from geometry_tools.algorithms.convex_surface_channel import earcut
 faces=f.copy();extra=[]
 for loop in boundary_loops(f):
  points=ov[loop];center=points.mean(0);_,_,axes=np.linalg.svd(points-center);xy=(points-center)@axes[:2].T;coords={int(i):q for i,q in zip(loop,xy)}
  extra.extend(earcut(loop[::-1].tolist(),coords))
 closed=np.r_[faces,np.asarray(extra,int)]
 labels,mask_check=module.regularize_labels(ov,closed,np.r_[covered.astype(int),np.zeros(len(extra),int)])
 covered=labels[:count].astype(bool);bare=set(np.flatnonzero(~covered).tolist())
 atomic_json(out/'mask-check.json',mask_check)
surfaces={'serosa:external':f[covered], 'serosa:subserosa':f[covered]+n};cells={'serosa':[dict(surface='serosa:external',side=1,neighbor=None),dict(surface='serosa:subserosa',side=-1,neighbor='subserosa')], 'subserosa':[dict(surface='serosa:subserosa',side=1,neighbor='serosa')], 'capsule':[]};external_neighbors={}
begin=0
for end in range(1,count+1):
 if end==count or covered[end]!=covered[begin]:
  sid='capsule:outer:'+str(begin);surfaces[sid]=f[begin:end]+2*n
  cells['capsule'].append(dict(surface=sid,side=1,neighbor='subserosa'if covered[begin]else None))
  if covered[begin]:cells['subserosa'].append(dict(surface=sid,side=-1,neighbor='capsule'))
  begin=end
for j,g in enumerate(p['tissueDrawGroups'][1:]):
 sid='capsule:inner:'+str(j);surfaces[sid]=m.faces[g['start']//3:(g['start']+g['count'])//3]+3*n;cells['capsule'].append(dict(surface=sid,side=1,neighbor=None));external_neighbors[sid]=g.get('neighbor')
for j,mouth in enumerate(boundary_loops(f[covered])):
 for k,name in enumerate(['serosa','subserosa']):
  sid='edge:'+name+':'+str(j);surfaces[sid]=annular_rim(mouth+k*n,mouth+(k+1)*n);cells[name].append(dict(surface=sid,side=1,neighbor=None))
ids={'serosa':'TA98PART:ZA:Liver:serosa','subserosa':'TA98PART:ZA:Liver:subserosa','capsule':p['id']}
result=pack(v,surfaces,cells,ids,m.frame,out/'packed',origin=[0,1.25,0],max_pairs=500000,max_tests=8000000,checkpoint_directory=out/'checks')
for g,use in zip(result['parts']['capsule']['tissueDrawGroups'],cells['capsule']):
 if external_neighbors.get(use['surface']):g['neighbor']=external_neighbors[use['surface']]
atomic_json(out/'coats.json',dict(passed=result['passed'],parts=result['parts'],bank=result['bank'],bank_sha256=result['bank_sha256'],capsuleOuterFacesPreserveOrder=True,referenceSurfaceTriangleCount=count,bareAreaTriangles=len(bare),bareAreaFaces=sorted(bare),coveredTriangles=int(covered.sum()),serosa_thickness_m=.00001,subserosa_thickness_m=.00004,thicknessStatus='Illustrative normal human-tissue display parameters, not measured histology',sourceCapsule=m.geometry_sha,external_neighbors=external_neighbors,mask_check=mask_check,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(passed=result['passed'],failures=result['checks']['failures'],bytes=result['total_bytes'],elapsed_seconds=time.monotonic()-start)),flush=True)
if not result['passed']:raise ValueError('Coat candidate failed; not installed')
