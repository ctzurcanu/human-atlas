"""Complementary finite caudate processes and fibrous left-tip tissue.

All cuts propagate through the shared gross/coating assembly; no extra organ
is superimposed on its parent. Anatomical regional extents remain provisional.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.tissue import pack,normal_words,decode_packed
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank,split_closed_cell
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces
ap=argparse.ArgumentParser();ap.add_argument('--canonical-cuts',action='store_true');ap.add_argument('--preserve-hila',action='store_true');ap.add_argument('--normal-guard',action='store_true');args=ap.parse_args()
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/('final-processes-v223'if args.normal_guard else'final-processes-v221'if args.preserve_hila else'final-processes-v220'if args.canonical_cuts else'final-processes-v218');out.mkdir(exist_ok=False);start=time.monotonic()
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cat.read_text());reader=AtlasReader();frame='model:local-ta98:assembled-metres';old={p['id']:p for p in a['parts']}
coat=json.loads((b/'final-coats-v216/coats.json').read_text());assert coat['passed']
cp={'chunks':[dict(url='/local-models/'+str(Path(coat['bank']).relative_to(ROOT/'.local-models')))]}
meshes={name:reader.read_part(a,old['TA98PART:ZA:Liver:'+name],frame=frame,catalogue_path=cat)for name in ['right','left','quadrate','caudate']}
for name,row in coat['parts'].items():meshes[name]=reader.read_part(cp,{**row,'chunk':0},frame=frame)
names=list(meshes);allv=np.concatenate([meshes[n].vertices for n in names]);v,inv=np.unique(allv,axis=0,return_inverse=True);offset=0;faces={};ancestry={}
for name in names:
 m=meshes[name];faces[name]=inv[offset:offset+len(m.vertices)][m.faces];ancestry[name]=np.arange(len(m.faces));offset+=len(m.vertices)
bank=SharedCutBank(v);cuts=[]
def fragment_faces(f,normal,offset,plane):
 length=np.linalg.norm(normal);normal=normal/length;offset=offset/length;points=np.asarray(bank.positions);distance=points@normal-offset;result=[];parents=[]
 for fi,tri in enumerate(f):
  if np.all(abs(distance[tri])<1e-12):result.append(tri.tolist());parents.append(fi);continue
  for side in [-1,1]:
   poly=[]
   for x,y in zip(tri,np.roll(tri,-1)):
    dx,dy=distance[x],distance[y]
    if side*dx>=-1e-12:poly.append(int(x))
    if dx*dy<0 and min(abs(dx),abs(dy))>1e-12:poly.append(bank.intersection(x,y,normal,offset,plane))
   poly=list(dict.fromkeys(poly))
   if poly:
    k=poly.index(min(poly));poly=poly[k:]+poly[:k]
   for j in range(1,len(poly)-1):result.append([poly[0],poly[j],poly[j+1]]);parents.append(fi)
 return np.asarray(result,int),np.asarray(parents,int)
plans=[('left','fibrous-appendix',[1,0,0],.094,1),('caudate','caudate-process',[.5,1,0],1.252,-1),('caudate','papillary-process',[-.8,1,0],1.266,-1)]
for target,child,normal,offset,child_side in plans:
 normal=np.asarray(normal,float);evidence={}
 for name in list(faces):
  f=faces[name];points=np.asarray(bank.positions);d=points[np.unique(f)]@normal-offset
  if d.min()>=0 or d.max()<=0:continue
  prior=ancestry[name]
  if args.preserve_hila and name!=target:
   faces[name],parent=fragment_faces(f,normal,offset,child);ancestry[name]=prior[parent];continue
  try:split=split_closed_cell(bank,f,normal,offset,child,allow_holes=True)
  except Exception as exc:
   atomic_json(out/'failed-cut.json',dict(cell=name,plane=child,error=str(exc)));raise
  def mapped(parent):
   ids=np.full(len(parent),-1,int);take=parent>=0;ids[take]=prior[parent[take]];return ids
  if name==target:
   lower,upper=split.negative_faces,split.positive_faces;la,ua=mapped(split.negative_ancestry),mapped(split.positive_ancestry)
   faces[child],ancestry[child]=(upper,ua)if child_side==1 else(lower,la)
   faces[name],ancestry[name]=(lower,la)if child_side==1 else(upper,ua)
   evidence=split.evidence
  else:
   nf=np.r_[split.negative_faces,split.positive_faces];pa=np.r_[split.negative_ancestry,split.positive_ancestry];take=pa>=0;order=np.argsort(pa[take],kind='stable');faces[name]=nf[take][order];ancestry[name]=prior[pa[take][order]]
 assert child in faces and len(faces[child])>12
 if args.canonical_cuts:
  # The stored representation is Float32. Canonicalize the actual emitted
  # coordinates between stages, never proximity-weld anatomical controls.
  origin=np.array([0,1.25,0]);encoded=(np.asarray(bank.positions)-origin).astype('<f4').astype(float)+origin
  quantized,remap=np.unique(encoded,axis=0,return_inverse=True)
  for name in faces:faces[name]=remap[faces[name]]
  bank=SharedCutBank(quantized)
 cuts.append(dict(parent=target,child=child,normal=normal.tolist(),offset=offset,**evidence));print(child,len(faces[child]),flush=True)
v=np.asarray(bank.positions);owners=defaultdict(list);outer_count=coat['referenceSurfaceTriangleCount']
for name,f in faces.items():
 for i,t in enumerate(f):owners[tuple(sorted(map(int,t)))].append((name,t,int(ancestry[name][i])))
grouped=defaultdict(list);host_parents=defaultdict(list)
for key,uses in owners.items():
 assert 1<=len(uses)<=2 and len({u[0]for u in uses})==len(uses)
 uses.sort(key=lambda x:x[0]);first=uses[0];second=uses[1]if len(uses)==2 else None
 host=next((u for u in uses if u[0]=='capsule'),None);is_outer=host is not None and 0<=host[2]<outer_count
 group=(first[0],second[0]if second else'',bool(is_outer));grouped[group].append(first[1]);host_parents[group].append(host[2]if host else-1)
surfaces={};cells={name:[]for name in faces};group_for={}
for j,(group,tris)in enumerate(grouped.items()):
 first,second,outer=group;sid=('outer:'if outer else'chart:')+str(j);surfaces[sid]=np.asarray(tris);group_for[sid]=group
 cells[first].append(dict(surface=sid,side=1,neighbor=second or None))
 if second:cells[second].append(dict(surface=sid,side=-1,neighbor=first))
cells['capsule'].sort(key=lambda u:0 if u['surface'].startswith('outer:')else 1)
surfaces,_,conform=conform_quantized_surfaces(v,surfaces,[0,1.25,0])
ids={name:'TA98PART:ZA:Liver:'+name for name in cells}
if args.preserve_hila:
 origin=np.array([0,1.25,0]);encoded=(v-origin).astype('<f4').astype(float);words=normal_words(encoded,surfaces)
 if args.normal_guard:
  for sid,f in surfaces.items():
   tri=encoded[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);bad=np.einsum('tki,ti->tk',words[sid].astype(float),fn)<=0
   if np.any(bad):
    flat=np.rint(fn/np.linalg.norm(fn,axis=1)[:,None]*32767).astype('<i2');words[sid][bad]=np.broadcast_to(flat[:,None],words[sid].shape)[bad]
 packed=pack_finite_cells(v,surfaces,words,cells,ids,origin);decoded=decode_packed(packed,frame,origin);checks={}
 for name,m in decoded.items():
  checks[name]=audit(m,require_closed=name not in ['right','quadrate'],max_pairs=1000000,max_tests=12000000);print(name,checks[name]['passed'],checks[name]['failures'],flush=True)
 directory=out/'packed';directory.mkdir();bank_path=directory/'cells.bin';atomic_bytes(bank_path,packed.data)
 for name,m in decoded.items():m.save(directory/(name+'.npz'))
 failures=[name for name,q in checks.items()if not q['passed']]
 result=dict(passed=not failures,parts=packed.cells,bank=str(bank_path),bank_sha256=digest_file(bank_path),checks=dict(failures=failures,cells=checks),cells=cells,source_hilar_boundaries_preserved=['right','quadrate'])
 atomic_json(directory/'tissue.json',result)
else:result=pack(v,surfaces,cells,ids,frame,out/'packed',origin=[0,1.25,0],max_pairs=1000000,max_tests=12000000,checkpoint_directory=out/'checks')
# Each capsule face can be rebound to its original source facet for labels.
host_sources=[]
for use in cells['capsule']:
 group=group_for[use['surface']];parents=host_parents[group]
 assert len(parents)==len(surfaces[use['surface']]),'Conformity changed host facet cardinality; rebind before install'
 host_sources.extend(parents)
atomic_json(out/'processes.json',dict(passed=result['passed'],parts=result['parts'],bank=result['bank'],bank_sha256=result['bank_sha256'],cuts=cuts,host_source_faces=host_sources,referenceSurfaceTriangleCount=sum(0<=i<outer_count for i in host_sources),source_reference_faces=outer_count,positions_literal_before_new_cut_controls=True,regional_boundaries_provisional=True,conform=conform,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(passed=result['passed'],failures=result['checks']['failures'],parts=len(cells),elapsed_seconds=time.monotonic()-start)),flush=True)
if not result['passed']:raise ValueError('Finite process partition failed; not installed')
