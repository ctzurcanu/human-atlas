"""Retain a complete shared cardiac arrangement while removing source overlap."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.finite_tissue_boolean import same_oriented_triangles
from geometry_tools.algorithms.conforming_assembly_ownership import subtract_assembly
from geometry_tools.tissue import validate_cells
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--hosts',type=int,nargs='+',required=True);ap.add_argument('--cutters',type=int,nargs='+',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.input;out=h/args.stage;out.mkdir(exist_ok=False);start=time.monotonic();surfaces={};fields={};textures={};cells={}
if(src/'ownership.json').exists():
 meta=json.loads((src/'ownership.json').read_text());z=np.load(src/'arrangement.npz');p=z['positions'];cells=meta['cells'];sources=meta['sources']
 for r in meta['surfaces']:surfaces[r['id']]=z[r['faces']];fields[r['id']]=z[r['normals']];textures[r['id']]=z[r['uvs']]
else:
 meta=json.loads((src/'emission.json').read_text());assert json.loads((src/'checked-geometry.json').read_text())['all_emitted_cells_checked'];assert all(c['passed']for c in json.loads((src/'checked-geometry.json').read_text())['cells'].values());assert digest_file(src/'cells.bin')==meta['bank_sha256'];sources=meta['sources'];reader=AtlasReader();mini=dict(chunks=[dict(url='/local-models/'+str((src/'cells.bin').relative_to(ROOT/'.local-models')))]);meshes={pid:reader.read_part(mini,dict(r,chunk=0,geometryOrigin=[0.,0.,0.]),frame='heart-native-local-metres')for pid,r in meta['parts'].items()};offsets={};all_points=[];n=0
 for pid,m in meshes.items():offsets[pid]=n;all_points.extend(m.vertices);n+=len(m.vertices)
 p,inv=np.unique(np.array(all_points),axis=0,return_inverse=True);shared={}
 for pid in sorted(meshes):
  m=meshes[pid];record=meta['parts'][pid];f=inv[offsets[pid]+m.faces];cells[pid]=[]
  for gi,g in enumerate(record['tissueDrawGroups']):
   sl=slice(g['start']//3,(g['start']+g['count'])//3);actual=f[sl];peer=g.get('neighbor');pair=frozenset([pid,peer])if peer else None
   if pair and pair in shared:
    sid=shared[pair];assert same_oriented_triangles(actual,surfaces[sid][:,[0,2,1]]);cells[pid].append(dict(surface=sid,side=-1,neighbor=peer));continue
   sid=f'{pid}/chart-{gi}';surfaces[sid]=actual;words=m.attributes['normals'][m.faces[sl]].copy();signs=m.attributes.get('normalSigns')
   if signs is not None:words*=signs[m.faces[sl],None]
   fields[sid]=words.astype('<i2');textures[sid]=m.attributes['uvs'][m.faces[sl]].copy();cells[pid].append(dict(surface=sid,side=1,**({'neighbor':peer}if peer else{})))
   if pair:shared[pair]=sid
def group(numbers):return[r['id']for s in sources if s['source']in numbers for r in s['parts'].values()]
result=subtract_assembly(p,surfaces,cells,group(args.hosts),group(args.cutters));q=result['positions'];arrays=dict(positions=q);definitions=[]
for i,(sid,faces)in enumerate(result['surfaces'].items()):
 chart=result['charts'][sid];words=np.empty((len(faces),3,3),dtype='<i2');uv=np.empty((len(faces),3,2),dtype='<f4')
 for ti,(family,parent,side,bary)in enumerate(zip(chart['source_surfaces'],chart['source_faces'],chart['normal_signs'],chart['barycentrics'])):
  values=bary@fields[family][parent].astype(float);values/=np.linalg.norm(values,axis=1)[:,None];words[ti]=np.rint(values*32767).astype('<i2');uv[ti]=(bary@textures[family][parent]).astype('<f4');tri=p[surfaces[family][parent]]
  for k in range(3):
   exact=np.flatnonzero(np.all(tri==q[faces[ti,k]],axis=1))
   if len(exact):words[ti,k]=fields[family][parent,int(exact[0])];uv[ti,k]=textures[family][parent,int(exact[0])]
  words[ti]*=side
 arrays[f'surface_{i}_faces']=faces;arrays[f'surface_{i}_normals']=words;arrays[f'surface_{i}_uvs']=uv;definitions.append(dict(id=sid,faces=f'surface_{i}_faces',normals=f'surface_{i}_normals',uvs=f'surface_{i}_uvs'))
decoded={name:Mesh(q,np.concatenate([result['surfaces'][u['surface']]if u['side']==1 else result['surfaces'][u['surface']][:,[0,2,1]]for u in uses]),'heart-native-local-metres')for name,uses in result['cells'].items()};checks=validate_cells(decoded,max_pairs=1500000,max_tests=16000000,checkpoint_directory=out/'checkpoints');atomic_json(out/'checked-f64-geometry.json',checks);atomic_npz(out/'arrangement.npz',**arrays);report=dict(passed=checks['passed'],source_stage=args.input,source_metadata_sha256=digest_file(src/('ownership.json'if(src/'ownership.json').exists()else'emission.json')),hosts=args.hosts,cutters=args.cutters,sources=sources,surfaces=definitions,cells=result['cells'],evidence=result['evidence'],elapsed_seconds=time.monotonic()-start,installed=False);atomic_json(out/'ownership.json',report);print(json.dumps(dict(passed=checks['passed'],failures=checks['failures'],elapsed_seconds=report['elapsed_seconds'])),flush=True)
