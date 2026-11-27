"""Staged PG covering from original finite-face ancestry, never fragment planes.

Only new local candidate files are written. No generic Boolean predicate,
source catalogue, ordinary atlas, fixed anatomy or installation is modified.
"""
from __future__ import annotations
import argparse,json,sys,time
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance
from palatoglossal_original_provenance import original_controls
W=ROOT/'.local-models/geometry-workbench/tongue-completion';FRAME='model:local-ta98:assembled-metres'

def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def qwords(mesh):
    if 'normals'not in mesh.attributes:
        t=mesh.vertices[mesh.faces];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);n/=np.linalg.norm(n,axis=1)[:,None];return np.repeat(np.rint(n*32767).astype('<i2')[:,None,:],3,axis=1)
    return (mesh.attributes['normals'][mesh.faces].astype(np.int32)*mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices),dtype='i1'))[mesh.faces][:,:,None]).astype('<i2')
def qkey(q):return tuple(str(x)for x in q)
def checked_plan(output):
    p=read(output/'plan.json')
    for path,sha in p['bindings'].items():require(digest_file(path)==sha,'Bound input/code changed: '+path)
    return p
def plan(output):
    require(output.is_relative_to(W)and not output.exists(),'Use one new original-parent stage directory')
    old=read(W/'palatoglossal-mucosal-entry-v2/plan.json');paths=[Path(p)for p in old['bindings']]
    paths += [Path(__file__).resolve(),Path(__file__).with_name('palatoglossal_incidence_fixture.py'),Path(__file__).with_name('palatoglossal_original_provenance.py'),ROOT/'geometry_tools/algorithms/finite_tissue_boolean.py',ROOT/'geometry_tools/algorithms/exact_source_intersections.py']
    source_files=dict(mucosa=old['mucosa'],whole=old['whole'],palate=old['palate'],arch_left=str(W/'palatoglossal-mucosal-entry-v2/left-uncut-expanded-arch.npz'),arch_right=str(W/'palatoglossal-mucosal-entry-v2/right-uncut-expanded-arch.npz'))
    paths += [Path(p)for p in source_files.values()]
    output.mkdir(parents=True);p=dict(schema='PG-original-parent-arrangement-plan/v1',source_files=source_files,addon=old['addon'],aponeurosis=old['aponeurosis'],superior_longitudinal=old['superior_longitudinal'],authority_inventory=str(ROOT/'TA98-PALATOGLOSSAL-ARCH.md'),
        canonical_arch_term='A05.2.01.005',canonical_arch_parent='A05.2.01.003',tongue_registry_terms_added=0,covering_estimate_m=.00035,measured_histology=False,
        method='Retain original source IDs through both Boolean operations; rational finite controls then one arrangement per identical original source face',
        fixed=['PG meshes','Whole tongue','Native palate','Skeleton','Unchanged source mucosal faces and reciprocal apo chart'],
        expressions=dict(arch='(expandedPG - originalM) - (PG union palate union (whole - originalM))',tongue='originalM - (PG union palate union (whole - originalM))'),
        bindings={str(p.resolve()):digest_file(p)for p in paths},installed=False,anatomical_acceptance=False)
    atomic_json(output/'plan.json',p);print(json.dumps(dict(stage='plan',output=str(output))))

def raw_stage(output):
    began=time.monotonic();p=checked_plan(output);require(not(output/'raw.json').exists(),'Do not replace a completed stage')
    meshes={name:Mesh.load(path)for name,path in p['source_files'].items()};addon=read(p['addon']);reader=AtlasReader()
    for row in addon['parts']:meshes['PG_'+row['side']]=reader.read_part(addon,row,frame=FRAME)
    sources={};solids={};source_ids={};source_words={};arrays={};source_meta={}
    for i,(name,m)in enumerate(meshes.items()):
        word=qwords(m);solid,oid=make_solid(m.vertices,m.faces,properties=word);solids[name]=solid;sources[oid]=m;source_ids[oid]=name;source_words[name]=word
        arrays[name+'_positions']=m.vertices;arrays[name+'_faces']=m.faces;arrays[name+'_words']=word;source_meta[name]=dict(positions=name+'_positions',faces=name+'_faces',words=name+'_words',geometry_sha256=m.geometry_sha)
    inherited={};intermediate={};provenance_events=[]
    def retained(name,solid):
        require(time.monotonic()-began<45,'Raw stage approached finite deadline')
        rational,events=original_controls(solid.to_mesh64(),sources,inherited,stage=name)
        intermediate[name]=dict(controls=len(rational),maximum_replay_m=max((r['replay_m']for r in events),default=0.),retained_prior_events=sum(r['method']=='literal-preserved-intermediate-event'for r in events))
        provenance_events.extend(dict(**{k:v for k,v in r.items()if k!='point'},point=qkey(r['point']))for r in events)
        return solid
    core=retained('inner-core',solids['whole']-solids['mucosa'])
    occupied=retained('bilateral-PG-union',solids['PG_left']+solids['PG_right'])
    palatal=retained('PG-palate-union',occupied+solids['palate'])
    cutter=retained('PG-palate-inner-core-union',palatal+core)
    left=retained('left-owner-difference',solids['arch_left']-solids['mucosa'])
    right=retained('right-owner-difference',solids['arch_right']-solids['mucosa'])
    cells=dict(arch_left=left-cutter,arch_right=right-cutter,tongue=solids['mucosa']-cutter)
    face_rows=defaultdict(list);checks={};controls={};family_counts=defaultdict(int)
    for cell,solid in cells.items():
        require(time.monotonic()-began<45,'Raw stage approached finite deadline')
        raw=solid.to_mesh64();v=np.asarray(raw.vert_properties)[:,:3];f=np.asarray(raw.tri_verts);oid,pid=provenance(raw)
        rational,replay=original_controls(raw,sources,inherited,stage='final:'+cell)
        cell_q={tuple(point):q for point,q in rational.items()};triangles=[]
        for face,original,parent in zip(f,oid,pid):
            name=source_ids[int(original)];source=sources[int(original)];sf=int(parent);q=tuple(cell_q[tuple(point)]for point in v[face]);bank=tuple(qkey(point)for point in q)
            t=np.array([[float(x)for x in point]for point in q]);st=source.vertices[source.faces[sf]];sign=1 if np.dot(np.cross(t[1]-t[0],t[2]-t[0]),np.cross(st[1]-st[0],st[2]-st[0]))>0 else-1
            face_rows[name].append(dict(cell=cell,source_face=sf,sign=sign,corners=bank));family_counts[name]+=1
            for key,qpoint in zip(bank,q):controls[key]=qpoint
            triangles.append(t)
        checks[cell]=dict(faces=len(f),original_source_reconstruction=True,maximum_original_ancestry_replay_m=max((r['replay_m']for r in replay),default=0.),manifold_volume_m3=float(solid.volume()))
    arrays['rational_controls_float']=np.array([[float(x)for x in q]for q in controls.values()]);atomic_npz(output/'original-sources.npz',**arrays)
    families={}
    for family,rows in face_rows.items():
        path=output/('raw-'+family+'.json');atomic_json(path,dict(family=family,source_faces=rows));families[family]=dict(path=str(path),sha256=digest_file(path),faces=len(rows),parents=len(set(r['source_face']for r in rows)))
    atomic_json(output/'intermediate-exact-events.json',dict(operations=intermediate,events=provenance_events))
    result=dict(schema='PG-original-parent-raw/v1',plan_sha256=digest_file(output/'plan.json'),original_sources=str(output/'original-sources.npz'),original_sources_sha256=digest_file(output/'original-sources.npz'),sources=source_meta,families=families,cells=checks,intermediate_operations=intermediate,
        every_boolean_face_has_original_finite_parent=True,rounded_fragments_reimported=False,distinct_original_faces_merged=False,shared_arrangement_pending=True,qualifying_candidate=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'raw.json',result);print(json.dumps(result))

if __name__=='__main__':
    q=argparse.ArgumentParser(description=__doc__);q.add_argument('--stage',choices=['plan','raw'],required=True);q.add_argument('--output',type=Path,required=True);a=q.parse_args();out=a.output.resolve()
    (plan if a.stage=='plan'else raw_stage)(out)
