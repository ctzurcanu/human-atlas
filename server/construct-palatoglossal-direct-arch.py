"""One capped PG arch cut with direct actual inner-sheet classifier ancestry.

No whole-minus-M re-Boolean, source remeshing or generic predicate changes.
Keeps original finite source faces, including the direct classifier's M-owned
inner/rim faces, through every Boolean operation. New local stages only.
"""
import argparse,importlib.util,json,sys,time
from collections import defaultdict
from pathlib import Path
import manifold3d as md
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance
from palatoglossal_original_provenance import original_controls
W=ROOT/'.local-models/geometry-workbench/tongue-completion';FRAME='model:local-ta98:assembled-metres'
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def words(m):return (m.attributes['normals'][m.faces].astype(np.int32)*m.attributes.get('normalSigns',np.ones(len(m.vertices),dtype='i1'))[m.faces][:,:,None]).astype('<i2')
def qkey(q):return tuple(str(x)for x in q)
def direct_classifier(negative,mucosa_id,basal_id,source_mucosa_faces):
    positions=negative.triangles.reshape(-1,3);props=words(negative).reshape(-1,3);n=len(positions);vp=np.c_[positions,props];f=np.arange(n,dtype=np.uint64).reshape(-1,3)
    _,first,inverse=np.unique(positions,axis=0,return_index=True,return_inverse=True);target=first[inverse];source=np.flatnonzero(target!=np.arange(n))
    count=len(source_mucosa_faces);face_ids=np.r_[source_mucosa_faces,np.arange(len(f)-count)].astype(np.uint64)
    mesh=md.Mesh64(np.ascontiguousarray(vp),np.ascontiguousarray(f),run_index=np.array([0,count*3,n],np.uint64),run_original_id=np.array([mucosa_id,basal_id],np.uint32),face_id=face_ids,
        merge_from_vert=source.astype(np.uint64),merge_to_vert=target[source].astype(np.uint64))
    solid=md.Manifold(mesh);require(solid.status()==md.Error.NoError,'Literal direct classifier was rejected by Manifold');return solid
def plan(output):
    require(output.is_relative_to(W)and not output.exists(),'Use a new direct arch stage directory')
    prior=W/'palatoglossal-original-parent-v2/plan.json';p=read(prior);negative=W/'palatoglossal-direct-inner-domain-v1/evidence.json';rim=W/'palatoglossal-direct-rim-proof-v1/evidence.json';n=read(negative);r=read(rim)
    require(n['classifier_geometry_passed']and r['passed']and r['off_rim_points_or_segments']==0 and r['basal_arch_contacts']==0,'Direct domain boundary proof failed')
    paths=[prior,negative,rim,Path(n['negative_mesh']),W/'palatoglossal-direct-inner-domain-v1/basal-classifier-only.npz',Path(__file__).resolve(),Path(__file__).with_name('palatoglossal_original_provenance.py'),Path(__file__).with_name('palatoglossal_incidence_fixture.py'),ROOT/'geometry_tools/algorithms/exact_source_intersections.py',ROOT/'geometry_tools/algorithms/finite_tissue_boolean.py']
    bindings={**p['bindings'],**r['source_bindings'],**{str(f):digest_file(f)for f in paths}}
    output.mkdir(parents=True);atomic_json(output/'plan.json',dict(schema='PG-direct-original-arch-plan/v1',sources={k:v for k,v in p['source_files'].items()if k!='whole'},addon=p['addon'],direct_classifier=str(negative),rim_proof=str(rim),
        canonical_arch_term='A05.2.01.005',canonical_arch_parent='A05.2.01.003',tongue_registry_terms_added=0,bindings=bindings,original_replay_bound_m=2e-10,
        fixed=['PG','Mucosal source except explicit PG entry cuts','Native palate','Whole tongue','Skeleton','Literal source APO chart'],installed=False,anatomical_acceptance=False));print(json.dumps(dict(stage='plan',output=str(output))))
def cut(output):
    began=time.monotonic();p=read(output/'plan.json');require(not(output/'raw.json').exists(),'Do not replace a completed direct arch stage')
    for path,sha in p['bindings'].items():require(digest_file(path)==sha,'Bound source/code changed: '+path)
    meshes={name:Mesh.load(path)for name,path in p['sources'].items()};addon=read(p['addon']);reader=AtlasReader()
    for row in addon['parts']:meshes['PG_'+row['side']]=reader.read_part(addon,row,frame=FRAME)
    source={};family={};solids={};ids={};arrays={};source_meta={};events=[];operations={};inherited={};stage='input'
    for name,m in meshes.items():
        solid,oid=make_solid(m.vertices,m.faces,properties=words(m));solids[name]=solid;source[oid]=m;family[oid]=name;ids[name]=oid
        arrays[name+'_positions']=m.vertices;arrays[name+'_faces']=m.faces;arrays[name+'_words']=words(m);source_meta[name]=dict(positions=name+'_positions',faces=name+'_faces',words=name+'_words',geometry_sha256=m.geometry_sha)
    negative_meta=read(p['direct_classifier']);negative=Mesh.load(negative_meta['negative_mesh']);basal=Mesh.load(W/'palatoglossal-direct-inner-domain-v1/basal-classifier-only.npz');basal_id=md.Manifold.reserve_ids(1);source[basal_id]=basal;family[basal_id]='classifier-only-basal'
    actual_m_faces=np.r_[negative_meta['original_inner_source_faces'],negative_meta['original_rim_source_faces']];core=direct_classifier(negative,ids['mucosa'],basal_id,actual_m_faces)
    def retained(name,solid):
        nonlocal stage
        stage=name;require(time.monotonic()-began<45,'Direct arch cut approached finite deadline')
        rational,ev=original_controls(solid.to_mesh64(),source,inherited,stage=name,maximum_replay_m=p['original_replay_bound_m'])
        operations[name]=dict(controls=len(rational),maximum_replay_m=max((r['replay_m']for r in ev),default=0.),inherited_events=sum(r['method']=='literal-preserved-intermediate-event'for r in ev))
        events.extend(dict(**{k:v for k,v in r.items()if k!='point'},point=qkey(r['point']))for r in ev);return solid
    try:
        retained('literal-inner-classifier',core)
        pg_union=retained('bilateral-PG-union',solids['PG_left']+solids['PG_right'])
        pg_pal=retained('PG-palate-union',pg_union+solids['palate'])
        cutter=retained('PG-palate-direct-inner-union',pg_pal+core)
        left=retained('left-owner-partition',solids['arch_left']-solids['mucosa']);right=retained('right-owner-partition',solids['arch_right']-solids['mucosa'])
        cells=dict(arch_left=left-cutter,arch_right=right-cutter,tongue=solids['mucosa']-cutter);by_family=defaultdict(list);cell_checks={}
        for cell,solid in cells.items():
            stage='final:'+cell;raw=solid.to_mesh64();v=np.asarray(raw.vert_properties)[:,:3];f=np.asarray(raw.tri_verts);oid,pid=provenance(raw)
            rational,ev=original_controls(raw,source,inherited,stage=stage,maximum_replay_m=p['original_replay_bound_m']);require(not np.any(oid==basal_id),'A classifier-only basal face contributed to retained tissue')
            triangles=[];corner_words=[]
            for face,original,parent in zip(f,oid,pid):
                original=int(original);parent=int(parent);name=family[original];m=source[original];q=tuple(rational[tuple(point)]for point in v[face]);t=np.array([[float(x)for x in point]for point in q]);st=m.triangles[parent]
                sign=1 if np.dot(np.cross(t[1]-t[0],t[2]-t[0]),np.cross(st[1]-st[0],st[2]-st[0]))>0 else-1
                by_family[name].append(dict(cell=cell,source_face=parent,sign=sign,corners=tuple(qkey(point)for point in q)));triangles.append(t)
            t=np.array(triangles);path=output/('raw-'+cell+'.npz');Mesh(t.reshape(-1,3),np.arange(t.size//3).reshape(-1,3),FRAME).save(path)
            cell_checks[cell]=dict(raw_mesh=str(path),raw_sha256=digest_file(path),faces=len(f),manifold_volume_m3=float(solid.volume()),maximum_original_replay_m=max((r['replay_m']for r in ev),default=0.))
    except Exception as error:
        atomic_json(output/'failed-stage.json',dict(schema='PG-direct-arch-exact-failure/v1',stage=stage,error=str(error),successful_operations=operations,elapsed_seconds=time.monotonic()-began,qualifying_candidate=False,installed=False));raise
    atomic_npz(output/'original-sources.npz',**arrays);families={}
    for name,rows in by_family.items():
        path=output/('raw-'+name+'.json');atomic_json(path,dict(family=name,source_faces=rows));families[name]=dict(path=str(path),sha256=digest_file(path),faces=len(rows),parents=len(set(r['source_face']for r in rows)))
    atomic_json(output/'exact-intermediate-events.json',dict(operations=operations,events=events));result=dict(schema='PG-direct-original-arch-raw/v1',plan_sha256=digest_file(output/'plan.json'),source_arrays=str(output/'original-sources.npz'),source_arrays_sha256=digest_file(output/'original-sources.npz'),sources=source_meta,families=families,cells=cell_checks,operations=operations,
        classifier_only_basal_faces_emitted=0,original_finite_parents_preserved=True,distance_welding=False,shared_arrangement_pending=True,qualifying_candidate=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'raw.json',result);print(json.dumps(result))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['plan','cut'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();(plan if a.stage=='plan'else cut)(a.output.resolve())
