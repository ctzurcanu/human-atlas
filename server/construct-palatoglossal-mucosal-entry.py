"""One finite, separately owned faucial-covering/tongue-entry construction.

The ordinary canonical palatoglossal arch is not a tongue muscle or a new
tongue term. This stage writes original finite-cell construction data only;
it never modifies the source anatomy, registry or catalogue and never installs.
"""
from __future__ import annotations
import argparse,json,sqlite3,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit,topology
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,subtract_occupied
from geometry_tools.algorithms.conforming_fields import replay_corner_fields

W=ROOT/'.local-models/geometry-workbench/tongue-completion';ORIGIN=np.array([0.,1.625,0.])
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def bound(p):return dict(path=str(Path(p).resolve()),sha256=digest_file(p))
def deadline():
    end=time.monotonic()+48
    def check():
        if time.monotonic()>end:raise TimeoutError('Finite arch construction reached its48second worker limit')
    return check
def literal_words(mesh):return mesh.attributes['normals'][mesh.faces].astype(np.int32)*mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices),dtype='i1'))[mesh.faces][:,:,None]
def normals(v,f):
    t=v[f];a=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);n=np.zeros_like(v)
    for i in range(3):np.add.at(n,f[:,i],a)
    n/=np.linalg.norm(n,axis=1)[:,None];return n
def decoded(solid,frame,source_ids):
    raw=solid.to_mesh64();v=np.asarray(raw.vert_properties,float);f=np.asarray(raw.tri_verts,int);orig,parents=provenance(raw);t=v[f,:3]
    g=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);g/=np.linalg.norm(g,axis=1)[:,None]
    word=np.repeat(np.rint(g*32767).astype('<i2')[:,None,:],3,axis=1)
    if v.shape[1]>=6:
        field=v[f,3:6];sign=np.array([source_ids.get(int(o),1)for o in orig]);field*=sign[:,None,None]
        length=np.linalg.norm(field,axis=2);good=(length>0)&((field*g[:,None,:]).sum(2)>0)
        word[good]=np.rint(np.clip(field[good],-32767,32767)).astype('<i2')
    return Mesh(t.reshape(-1,3),np.arange(f.size).reshape(-1,3),frame,dict(normals=word.reshape(-1,3))),dict(provenance_ids=orig.tolist(),parent_faces=parents.tolist())

def plan_stage(output):
    require(output.is_relative_to(W)and not output.exists(),'Use one new local arch stage')
    authority=ROOT/'.local-models/ta98-authority.sqlite';db=sqlite3.connect(authority);db.row_factory=sqlite3.Row
    ids=['A05.2.01.003','A05.2.01.004','A05.2.01.005','A05.2.01.006','A05.2.01.007','A05.2.01.008','A05.1.01.002','A05.1.04.011','A05.1.04.110']
    terms=[dict(db.execute('SELECT id,name_en,name_la,parent_id,parent_name,immaterial,type_of_entity FROM ta98 WHERE id=?',(term,)).fetchone())for term in ids]
    sources=read(W/'source-v1/source-audit.json');palate=next(r for r in sources['contextParts']if r['partId']=='ZA:Soft palate')
    bindings={str(authority):digest_file(authority),str(ROOT/'server/construct-palatoglossal-mucosal-entry.py'):digest_file(__file__)}
    addon=W/'palatoglossus-addon-v1/addon.json';a=read(addon);assert a['passed'];paths=[addon,Path(a['output_bank']['path']),Path(palate['meshSavePath']),W/'basal-whole-protected-v1/whole-tongue.npz']
    shared=W/'shared-coverings-exposed-apo-sl-v3/packed'
    paths += [shared/('TA98TONGUE_male_'+name+'.npz')for name in ['mucosa','aponeurosis','intrinsic_superior-longitudinal']]
    for path in paths:bindings[str(path)]=digest_file(path)
    for term in terms:
        term.update(representation_class='finite-composite-mucosal-fold'if term['id']=='A05.2.01.005'else'finite-mucosal-coat'if term['id']in ['A05.1.01.002','A05.1.04.011']else'finite-muscle-envelope'if term['id']=='A05.1.04.110'else'whole-organ'if term['id']=='A05.2.01.003'else'material-fold-or-region',
            existing_source_meshes=['ZA:Soft palate']if term['id']=='A05.2.01.003'else['TA98TONGUE:male:mucosa']if term['id']=='A05.1.04.011'else a['termGroup']['partIds']if term['id']=='A05.1.04.110'else[],
            planned_representation='Separate finite arch-covering owner'if term['id']=='A05.2.01.005'else'Local PG-entry conformation only; all other surfaces fixed'if term['id']=='A05.1.04.011'else'Unchanged source'if term['id']in ['A05.2.01.003','A05.1.04.110']else'Not constructed by this stage')
    plan=dict(schema='palatoglossal-arch-finite-entry-plan/v1',authority=bound(authority),terms=terms,term_count=len(terms),bindings=bindings,
        addon=str(addon),palate=palate['meshSavePath'],whole=str(W/'basal-whole-protected-v1/whole-tongue.npz'),mucosa=str(paths[-3]),aponeurosis=str(paths[-2]),superior_longitudinal=str(paths[-1]),
        canonical_arch_term='A05.2.01.005',canonical_arch_parent='A05.2.01.003',canonical_muscle_parent='A04.1.06.001',tongue_registry_terms_added=0,
        covering_estimate_m=.00035,covering_histology_measured=False,connective_and_mucosal_composite_not_separate_micro_layers=True,
        primary_references=['https://pubmed.ncbi.nlm.nih.gov/281278/','https://pubmed.ncbi.nlm.nih.gov/23686451/'],
        method='Original expanded reviewed-PG envelope minus PG occupied muscle, native palate and actual inner tongue compartment; subtract/conform existing tongue mucosa using one reciprocal Boolean chart',
        cells=dict(tongue='A05.1.04.011; existing tongue coat with physically open PG entry',arch='A05.2.01.005; independent faucial mucosal/connective covering, never tongue mucosa or muscle'),
        fixed=['Both original PG meshes','Native palate geometry','Whole tongue source','Skeleton','Tongue–aponeurosis reciprocal facet geometry and normal fields'],
        adjoining_interfaces=['Arch↔tongue mucosa: literal shared reciprocal interface, one chart','Arch↔PG and tongue mucosa↔PG: finite covering bounds around muscular core; no disk/plug across muscle',
            'Palatal termination within fixed native palate: physical covering termination only; histological enthesis remains unestablished'],
        missing_sources=['No donor-specific pillar thickness measurement','Triangular fold extent unprovided; no representation inferred','Current native palate is an undivided whole, not separate palatal mucosal micro layers'],
        installed=False,anatomical_acceptance=False)
    output.mkdir(parents=True);atomic_json(output/'plan.json',plan)
    report=ROOT/'TA98-PALATOGLOSSAL-ARCH.md'
    table='\n'.join('| '+r['id']+' | '+r['name_en']+' | '+r['parent_id']+' | '+r['representation_class']+' | '+r['planned_representation']+' |'for r in terms)
    if not report.exists():
        report.write_text('# Palatoglossal arch: authority and finite tissue contract\n\nAuthority: local TA98 database SHA256 '+digest_file(authority)+'. No source geometry is modified or installed by this stage.\n\n| Term | Name | Canonical parent | Representation | Stage scope |\n|---|---|---|---|---|\n'+table+'\n\nThe arch remains a soft-palate child, separate from the 36 tongue terms. The estimated 0.35 mm composite covering is not measured histology. The triangular fold remains unrepresented.\n\nAdult cadaver references: https://pubmed.ncbi.nlm.nih.gov/281278/ and https://pubmed.ncbi.nlm.nih.gov/23686451/. The muscle lies inside the connective/mucosal pillar. The finite arch and tongue coats must have distinct owners and reciprocal shared coordinates/normals; no covering disk may cross the muscular core. Native palate, PG, tongue source and bones stay fixed. The existing tongue–aponeurosis seam must remain literal.\n\nMachine-readable inventory and exact input bindings: '+str(output/'plan.json')+'.\n')
    else:
        with report.open('a')as stream:stream.write('\nBound producer revision: '+str(output/'plan.json')+'. The initialv1 separate-cut ordering failed exact-source replay; this revision partitions owners first and cuts both simultaneously.\n')
    return dict(output=str(output),terms=len(terms),planned=True,installed=False)

def construct_stage(output):
    began=time.monotonic();check=deadline();plan=read(output/'plan.json');require(not(output/'construction.json').exists(),'Do not repeat a completed construction')
    for p,s in plan['bindings'].items():require(digest_file(p)==s,'A bound plan input changed: '+p)
    addon=read(plan['addon']);reader=AtlasReader();pg={p['side']:reader.read_part(addon,p,frame='model:local-ta98:assembled-metres')for p in addon['parts']}
    mucosa=Mesh.load(plan['mucosa']);whole=Mesh.load(plan['whole']);palate=Mesh.load(plan['palate'])
    m,mid=make_solid(mucosa.vertices,mucosa.faces,properties=literal_words(mucosa));w,_=make_solid(whole.vertices,whole.faces);pal,_=make_solid(palate.vertices,palate.faces)
    pgs={};pid={};outer={};outer_inputs={}
    for side,p in pg.items():
        pgs[side],pid[side]=make_solid(p.vertices,p.faces,properties=literal_words(p))
        v,first,alias=np.unique(p.vertices,axis=0,return_index=True,return_inverse=True);f=alias[p.faces];n=normals(v,f);expanded=v+plan['covering_estimate_m']*n
        em=Mesh(expanded,f,p.frame,dict(normals=np.rint(n*32767).astype('<i2')));em.save(output/(side+'-expanded-covering-envelope.npz'))
        ea=audit(em,require_closed=True,max_pairs=100000,max_tests=1000000);atomic_json(output/(side+'-expanded-preflight.json'),ea);require(ea['passed'],'Expanded composite arch envelope is invalid; do not tune/retry silently')
        outer[side],oid=make_solid(expanded,f,properties=em.attributes['normals'][em.faces]);outer_inputs[side]=(em,oid)
    occupied=pgs['left']+pgs['right'];core=w-m
    # Partition source owners before any PG cuts. Otherwise independently
    # rounded contact planes appear in the two derived host boundaries.
    mmesh=mucosa
    host_meshes={};check()
    for side in ['left','right']:
        host=outer_inputs[side][0];host.save(output/(side+'-uncut-expanded-arch.npz'))
        ah=audit(host,require_closed=True,max_pairs=100000,max_tests=1000000);atomic_json(output/(side+'-arch-preflight.json'),ah);require(ah['passed'],'Actual arch shell failed before ownership partition')
        host_meshes[side]=host
    # One canonical difference simultaneously returns arch-owned boundary
    # fragments and the reverse/free boundary of the occupied tongue coat.
    positions=np.vstack([host_meshes[s].vertices for s in ['left','right']]);offset=0;surfaces={};rows={};cells={};parents={};uses={};host_words={}
    for side in ['left','right']:
        host=host_meshes[side];sid='arch:'+side;surfaces[sid]=host.faces+offset;rows[sid]=dict(id=sid,kind='outer',array=sid);cells[side]=surfaces[sid];parents[side]=[(sid,i,1)for i in range(len(host.faces))];uses[side]=[dict(surface=sid,side=1)];host_words[sid]=literal_words(host).astype('<i2');offset+=len(host.vertices)
    check();result=subtract_occupied(positions,surfaces,rows,cells,parents,uses,mmesh.vertices,mmesh.faces,
        occupied_face_owners=['tongue-mucosa']*len(mmesh.faces),partition_free_outer=True,exact_source_replay=True,
        reconcile_shared_subdivisions=True,diagnostic_directory=output/'boolean-proof')
    # Replay the first owner partition before the one simultaneous anatomical
    # cut. Every shared facet has a single source chart at this point.
    chart_rows={r['id']:r for r in result['surfaceRows']};fields={};all_words={**host_words,'tongue:source':literal_words(mmesh).astype('<i2')};all_surfaces={**surfaces,'tongue:source':mmesh.faces+len(positions)}
    for sid,chart in result['charts'].items():
        source='tongue:source'if chart_rows[sid]['parentSurface']=='occupied-outer'else chart_rows[sid]['parentSurface'];sign=chart_rows[sid].get('parentNormalSign',1)
        fields[sid]=dict(source_surfaces=[source]*len(result['surfaces'][sid]),source_faces=chart['parentFaces'],normal_signs=np.full(len(result['surfaces'][sid]),sign),barycentrics=chart['barycentrics'])
    replay=dict(positions=result['positions'],surfaces=result['surfaces'],charts=fields);words,_=replay_corner_fields(np.vstack([positions,mmesh.vertices]),all_surfaces,all_words,{k:np.zeros((*fs.shape,2),'<f4')for k,fs in all_surfaces.items()},replay)
    cells_out=result['cells'];cells_out['tongue-mucosa']=[dict(surface=u['surface'],side=u['side'],neighbor=u.get('neighbor'))for u in result['vascularContacts']+result['vascularFreeSurfaces']]
    cut=occupied+pal+core;raw=cut.to_mesh64();cut_positions=np.asarray(raw.vert_properties,float)[:,:3];cut_faces=np.asarray(raw.tri_verts,int)
    cell_faces={};parent_faces={};source_rows={r['id']:r for r in result['surfaceRows']}
    for name,uses in cells_out.items():
        pieces=[];parent=[]
        for use in uses:
            fs=result['surfaces'][use['surface']];pieces.append(fs if use['side']==1 else fs[:,[0,2,1]])
            parent.extend((use['surface'],i,use['side'])for i in range(len(fs)))
        cell_faces[name]=np.vstack(pieces);parent_faces[name]=parent
    check();final=subtract_occupied(result['positions'],result['surfaces'],source_rows,cell_faces,parent_faces,cells_out,
        cut_positions,cut_faces,occupied_face_owners=['fixed-PG-palate-inner-tongue']*len(cut_faces),
        exact_source_replay=True,reconcile_shared_subdivisions=True,diagnostic_directory=output/'simultaneous-cut-proof')
    final_fields={};final_rows={r['id']:r for r in final['surfaceRows']};tri=cut_positions[cut_faces];gn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);gn/=np.linalg.norm(gn,axis=1)[:,None]
    cutter_words=np.repeat(np.rint(gn*32767).astype('<i2')[:,None,:],3,axis=1)
    for sid,chart in final['charts'].items():
        source='cutter'if final_rows[sid]['parentSurface']=='occupied-outer'else final_rows[sid]['parentSurface'];sign=final_rows[sid].get('parentNormalSign',1)
        final_fields[sid]=dict(source_surfaces=[source]*len(final['surfaces'][sid]),source_faces=chart['parentFaces'],normal_signs=np.full(len(final['surfaces'][sid]),sign),barycentrics=chart['barycentrics'])
    all_final_surfaces={**result['surfaces'],'cutter':cut_faces+len(result['positions'])};all_final_words={**words,'cutter':cutter_words}
    words,_=replay_corner_fields(np.vstack([result['positions'],cut_positions]),all_final_surfaces,all_final_words,
        {k:np.zeros((*fs.shape,2),'<f4')for k,fs in all_final_surfaces.items()},dict(positions=final['positions'],surfaces=final['surfaces'],charts=final_fields))
    final['evidence']['ownerPartition']=result['evidence'];result=final;cells_out=final['cells']
    arrays=dict(positions=result['positions']);meta_surfaces=[]
    for i,(sid,faces)in enumerate(result['surfaces'].items()):arrays['surface'+str(i)]=faces;arrays['words'+str(i)]=words[sid];meta_surfaces.append(dict(id=sid,faces='surface'+str(i),words='words'+str(i)))
    atomic_npz(output/'finite-arch-entry-arrangement.npz',**arrays)
    check();meta=dict(schema='finite-palatoglossal-entry-arrangement/v1',positions='positions',surfaces=meta_surfaces,cells=cells_out,source_plan_sha256=digest_file(output/'plan.json'),
        boolean_evidence=result['evidence'],muscle_geometry_modified=False,native_palate_modified=False,whole_tongue_modified=False,
        arch_term='A05.2.01.005',arch_canonical_parent='A05.2.01.003',tongue_mucosa_term='A05.1.04.011',covering_estimate_m=.00035,
        actual_raw_self_checks_pending=True,packing_pending=True,apo_reciprocal_seam_replay_pending=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'construction.json',meta)
    print(json.dumps(dict(output=str(output),saved_arrangement=True,surfaces=len(meta_surfaces),cells=list(cells_out),packing_pending=True,elapsed_seconds=meta['elapsed_seconds'])))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['plan','construct'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve()
    if a.stage=='plan':print(json.dumps(plan_stage(out)))
    else:construct_stage(out)
