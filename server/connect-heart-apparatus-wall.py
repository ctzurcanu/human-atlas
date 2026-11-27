"""Resolve current source basal/annular overlap without moving cardiac anatomy.

The whole fourteen-cell wall/root host and fifteen-cell AV apparatus retain
their reciprocal interfaces. This construction does not infer that every new
contact is the correct anatomical attachment; its ownership is reviewed next.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.finite_tissue_boolean import same_oriented_triangles
from geometry_tools.algorithms.conforming_assembly_ownership import subtract_assembly
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args()
start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());marker=a['finiteCardiacIntegration'];pp={p['id']:p for p in a['parts']};origin=np.array([0.,1.375,.03125]);frame='model:local-ta98:assembled-metres';reader=AtlasReader()
hosts=marker['wallPartIds']+marker['rootWallPartIds'];cutters=[pid for pid in marker['valvePartIds']+marker['chordalPartIds']+marker['papillaryPartIds']if pp[pid].get('cardiacApparatusOwnership')];assert len(hosts)==14 and len(cutters)==15
meshes={pid:reader.read_part(a,pp[pid],frame=frame,catalogue_path=path)for pid in hosts+cutters};offsets={};all_points=[];n=0
for pid,m in meshes.items():offsets[pid]=n;all_points.extend(m.vertices-origin);n+=len(m.vertices)
p,inv=np.unique(np.asarray(all_points),axis=0,return_inverse=True);surfaces={};fields={};textures={};cells={};shared={}
for pid in sorted(meshes):
 m=meshes[pid];f=inv[offsets[pid]+m.faces];cells[pid]=[]
 for gi,g in enumerate(pp[pid]['tissueDrawGroups']):
  sl=slice(g['start']//3,(g['start']+g['count'])//3);actual=f[sl];peer=g.get('neighbor');assert not peer or peer in meshes
  pair=frozenset([pid,peer])if peer else None
  if pair and pair in shared:
   sid=shared[pair];assert same_oriented_triangles(actual,surfaces[sid][:,[0,2,1]]);cells[pid].append(dict(surface=sid,side=-1,neighbor=peer));continue
  sid=f'{pid}/source-chart-{gi}';surfaces[sid]=actual;words=m.attributes['normals'][m.faces[sl]].copy();signs=m.attributes.get('normalSigns')
  if signs is not None:words*=signs[m.faces[sl],None]
  fields[sid]=words.astype('<i2');textures[sid]=m.attributes['uvs'][m.faces[sl]].copy()if'uvs'in m.attributes else np.zeros((len(actual),3,2),'<f4');cells[pid].append(dict(surface=sid,side=1,**({'neighbor':peer}if peer else{})))
  if pair:shared[pair]=sid
atomic_npz(out/'source-arrangement.npz',positions=p,**{f'faces_{i}':f for i,f in enumerate(surfaces.values())})
atomic_json(out/'source.json',dict(catalogue_sha256=digest_file(path),origin=origin.tolist(),hosts=hosts,cutters=cutters,cells=cells,surfaces=[dict(id=sid,faces=f'faces_{i}')for i,sid in enumerate(surfaces)],sources={pid:m.source for pid,m in meshes.items()}))
print('Decoded',len(meshes),'cells;',sum(map(len,surfaces.values())),'unique chart triangles',flush=True)
result=subtract_assembly(p,surfaces,cells,hosts,cutters);q=result['positions'];arrays=dict(positions=q);definitions=[]
print('Boolean ownership and shared-edge refinement complete',round(time.monotonic()-start,2),'s',flush=True)
for i,(sid,faces)in enumerate(result['surfaces'].items()):
 if time.monotonic()-start>110:raise TimeoutError('Attachment field stage deadline')
 chart=result['charts'][sid];refs=np.asarray(chart['source_surfaces']);parents=chart['source_faces'];sides=chart['normal_signs'];bary=chart['barycentrics'];words=np.empty((len(faces),3,3),'<i2');uv=np.empty((len(faces),3,2),'<f4')
 for family in np.unique(refs):
  ids=np.flatnonzero(refs==family);par=parents[ids];w=bary[ids];source_words=fields[family][par];values=np.einsum('tik,tkj->tij',w,source_words.astype(float));values/=np.linalg.norm(values,axis=2)[:,:,None];ww=np.rint(values*32767).astype('<i2');uu=np.einsum('tik,tkj->tij',w,textures[family][par]).astype('<f4');tri=p[surfaces[family][par]];exact=np.all(q[faces[ids]][:,:,None,:]==tri[:,None,:,:],axis=3)
  for corner in range(3):
   same=exact[:,corner].any(axis=1);j=exact[:,corner].argmax(axis=1);ti=np.flatnonzero(same);ww[ti,corner]=source_words[ti,j[ti]];uu[ti,corner]=textures[family][par[ti],j[ti]]
  words[ids]=ww*sides[ids,None,None];uv[ids]=uu
 arrays[f'surface_{i}_faces']=faces;arrays[f'surface_{i}_normals']=words;arrays[f'surface_{i}_uvs']=uv;definitions.append(dict(id=sid,faces=f'surface_{i}_faces',normals=f'surface_{i}_normals',uvs=f'surface_{i}_uvs'))
atomic_npz(out/'arrangement.npz',**arrays);contacts=[]
for pid in hosts:
 for u in result['cells'][pid]:
  if u.get('neighbor')in cutters:contacts.append(dict(host=pid,apparatus=u['neighbor'],surface=u['surface'],triangles=len(result['surfaces'][u['surface']])))
atomic_json(out/'ownership.json',dict(constructed=True,numerical_acceptance_pending=True,catalogue_sha256=digest_file(path),source_stage=args.stage,origin=origin.tolist(),hosts=hosts,cutters=cutters,surfaces=definitions,cells=result['cells'],contacts=contacts,evidence=result['evidence'],source_parts={pid:pp[pid]for pid in meshes},elapsed_seconds=time.monotonic()-start,installed=False))
print(json.dumps(dict(constructed=True,cells=len(result['cells']),contacts=contacts,elapsed_seconds=time.monotonic()-start)),flush=True)
