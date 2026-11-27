"""Bounded native equivalent of the canonical labelled-tetrahedra algorithm."""
from pathlib import Path
import hashlib
import json,struct,subprocess,time
import numpy as np
from .common import CACHE,atomic_json,digest_file

def native_executable():
 source=Path(__file__).with_name('native')/'labelled_tetrahedra.cpp'
 source_hash=digest_file(source);folder=CACHE/'native-labelled-tetrahedra'/source_hash;folder.mkdir(parents=True,exist_ok=True);binary=folder/'partition';manifest=folder/'compile.json'
 if binary.exists()and manifest.exists():
  m=json.loads(manifest.read_text())
  if m['source_sha256']==source_hash and m['binary_sha256']==digest_file(binary):return binary,m
 start=time.monotonic();subprocess.run(['c++','-std=c++17','-O2',str(source),'-o',str(binary)],check=True,timeout=20)
 m=dict(source_sha256=source_hash,binary_sha256=digest_file(binary),compiler='c++ -std=c++17 -O2',deadline_seconds=20,elapsed_seconds=time.monotonic()-start);atomic_json(manifest,m);return binary,m

def partition_native(positions,tetrahedra,scores,output,*,seconds=60.,maximum_tetrahedra=500000):
 p=np.asarray(positions,dtype='<f8');t=np.asarray(tetrahedra,dtype='<u8');s=np.asarray(scores,dtype='<f8');output=Path(output);output.mkdir(parents=True,exist_ok=False)
 if p.ndim!=2 or p.shape[1]!=3 or t.ndim!=2 or t.shape[1]!=4 or s.ndim!=2 or len(s)!=len(p):raise ValueError('Indexed positions, tetrahedra and score vectors required')
 if len(t)>maximum_tetrahedra or not 0<seconds<=300 or not np.isfinite(p).all()or not np.isfinite(s).all()or t.max()>=len(p):raise ValueError('Finite fields and explicit partition budgets required')
 binary,build=native_executable();input_file=output/'input.bin';result_file=output/'partition.bin'
 with input_file.open('wb')as stream:
  stream.write(b'ATLTET01'+struct.pack('<QQQ',len(p),len(t),s.shape[1]));stream.write(p.tobytes());stream.write(t.tobytes());stream.write(s.tobytes())
 start=time.monotonic();result=subprocess.run([str(binary),str(input_file),str(result_file),str(seconds)],stdout=subprocess.PIPE,check=True,timeout=seconds+2);native=json.loads(result.stdout)
 with result_file.open('rb')as stream:
  if stream.read(8)!=b'ATLPAR01':raise ValueError('Invalid native partition output')
  n,labels=struct.unpack('<QQ',stream.read(16));v=np.fromfile(stream,dtype='<f8',count=n*3).reshape(n,3);cells={}
  for label in range(labels):
   count=struct.unpack('<Q',stream.read(8))[0];cells[label]=np.fromfile(stream,dtype='<u8',count=count*3).reshape(count,3).astype('int64')
  if stream.read(1):raise ValueError('Trailing native partition data')
 if not np.array_equal(v[:len(p)],p):raise ValueError('Native partition changed source controls')
 evidence=dict(**native,native_backend=True,source_controls_literal=True,source_code_sha256=build['source_sha256'],binary_sha256=build['binary_sha256'],input_sha256=digest_file(input_file),output_sha256=digest_file(result_file),wall_deadline_seconds=seconds,process_elapsed_seconds=time.monotonic()-start,identical_affine_score_priority=list(range(s.shape[1])),anatomical_acceptance=False)
 atomic_json(output/'evidence.json',evidence);return dict(positions=v,cells=cells,evidence=evidence)

def partition_constraints_native(positions,tetrahedra,fields,cells,zero_side_priority,output,*,seconds=90.,maximum_tetrahedra=500000,tetrahedron_domains=None,cell_domains=None,surface_refinement_planes=()):
 """Intersect each tetrahedron with explicitly owned linear inequalities.

 ``cells`` contains (field index, positive/negative side) references. Zero
 side priority prevents two owners of a field identically zero over a whole
 calculation tetrahedron. The caller must establish coverage/nonoverlap of
 its anatomical classification and audit the actual emitted cells.
 """
 p=np.asarray(positions,dtype='<f8');t=np.asarray(tetrahedra,dtype='<u8');s=np.asarray(fields,dtype='<f8');zero=np.asarray(zero_side_priority,dtype='<i8');output=Path(output);output.mkdir(parents=True,exist_ok=False)
 if p.ndim!=2 or p.shape[1]!=3 or t.ndim!=2 or t.shape[1]!=4 or s.ndim!=2 or len(s)!=len(p)or zero.shape!=(s.shape[1],):raise ValueError('Indexed positions, tetrahedra, fields and zero-side priorities required')
 if len(t)>maximum_tetrahedra or not 0<seconds<=300 or not np.isfinite(p).all()or not np.isfinite(s).all()or t.max()>=len(p)or s.shape[1]>64 or not np.isin(zero,[-1,1]).all():raise ValueError('Finite constrained fields and explicit budgets required')
 domains=tetrahedron_domains is not None or cell_domains is not None
 if domains:
  td=np.asarray(tetrahedron_domains);cd=np.asarray(cell_domains)
  if td.shape!=(len(t),)or cd.shape!=(len(cells),)or not np.issubdtype(td.dtype,np.integer)or not np.issubdtype(cd.dtype,np.integer)or td.min()<0 or cd.min()<0 or td.max()>63 or cd.max()>63:raise ValueError('Explicit integer tetrahedral and cell material domains required')
  if set(map(int,td))!=set(map(int,cd)):raise ValueError('Every material domain requires calculation tetrahedra and owned cells')
  td=td.astype('<u8');cd=cd.astype('<u8')
 # Different anatomical classifications may use literally the same affine
 # comparison (for example a chamber's lining and its underlying muscle).
 # Give such a plane one canonical identity before clipping. No coefficient
 # or regional extent is changed and hash matches are checked literally.
 columns=[];priorities=[];aliases=[];lookup={}
 for i in range(s.shape[1]):
  column=s[:,i].copy();column[column==0]=0.;key=hashlib.sha256(column.tobytes()).hexdigest();index=lookup.get(key)
  if index is None:index=len(columns);lookup[key]=index;columns.append(column);priorities.append(int(zero[i]))
  elif not np.array_equal(column,columns[index])or zero[i]!=priorities[index]:raise ValueError('Incompatible coincident field priority')
  aliases.append(index)
 s=np.column_stack(columns);zero=np.array(priorities,dtype='<i8');cells=[list(dict.fromkeys((aliases[field],side)for field,side in cell))for cell in cells]
 if any(not isinstance(i,(int,np.integer))or not 0<=i<len(aliases)for i in surface_refinement_planes):raise ValueError('Explicit valid surface refinement field IDs required')
 refinements=list(dict.fromkeys(aliases[int(i)]for i in surface_refinement_planes))
 binary,build=native_executable();input_file=output/'input.bin';result_file=output/'partition.bin'
 with input_file.open('wb')as stream:
  mode=b'ATLPDR01'if domains and refinements else b'ATLPLR01'if refinements else b'ATLPDC01'if domains else b'ATLPLC01';stream.write(mode+struct.pack('<QQQQ',len(p),len(t),s.shape[1],len(cells)));stream.write(p.tobytes());stream.write(t.tobytes());stream.write(s.tobytes());stream.write(zero.tobytes())
  if domains:stream.write(td.tobytes());stream.write(cd.tobytes())
  if refinements:stream.write(struct.pack('<Q',len(refinements)));stream.write(np.asarray(refinements,dtype='<u8').tobytes())
  for cell in cells:
   stream.write(struct.pack('<Q',len(cell)))
   for field,side in cell:
    if not 0<=field<s.shape[1]or side not in [-1,1]:raise ValueError('Invalid cell constraint')
    stream.write(struct.pack('<Qq',field,side))
 start=time.monotonic();result=subprocess.run([str(binary),str(input_file),str(result_file),str(seconds)],stdout=subprocess.PIPE,check=True,timeout=seconds+2);native=json.loads(result.stdout)
 with result_file.open('rb')as stream:
  if stream.read(8)!=b'ATLPAR01':raise ValueError('Invalid native constrained output')
  n,labels=struct.unpack('<QQ',stream.read(16));v=np.fromfile(stream,dtype='<f8',count=n*3).reshape(n,3);result_cells={}
  for label in range(labels):
   count=struct.unpack('<Q',stream.read(8))[0];result_cells[label]=np.fromfile(stream,dtype='<u8',count=count*3).reshape(count,3).astype('int64')
  if stream.read(1):raise ValueError('Trailing constrained partition data')
 if not np.array_equal(v[:len(p)],p):raise ValueError('Native partition changed source controls')
 evidence=dict(**native,native_backend=True,mode='material-domain-cell-constraints'if domains else'explicit-piecewise-linear-cell-constraints',source_controls_literal=True,material_domains=sorted(set(map(int,td)))if domains else None,surface_refinement_planes=refinements,coincident_field_aliases=aliases,canonical_field_count=len(columns),source_code_sha256=build['source_sha256'],binary_sha256=build['binary_sha256'],input_sha256=digest_file(input_file),output_sha256=digest_file(result_file),wall_deadline_seconds=seconds,process_elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False)
 atomic_json(output/'evidence.json',evidence);return dict(positions=v,cells=result_cells,evidence=evidence)
