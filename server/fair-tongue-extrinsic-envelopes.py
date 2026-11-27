"""One constrained actual-surface fairing trial for four source-bound muscles.

Each command performs one saved finite stage. Superior insertion triangles and
their corner fields, cortical contact controls and cortical cap faces are fixed.
Only other physical tissue surfaces receive an adaptive refinement and a single
discrete biharmonic solve. No skeleton, source bank or catalogue is modified.
"""
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit,topology
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
FRAME='model:local-ta98:assembled-metres'
ORIGIN=np.array([0.,1.625,0.])
PARTS={'GG-left':'ZA:Genioglossus muscle.l','GG-right':'ZA:Genioglossus muscle.r',
       'HG-left':'ZA:Hyoglossus muscle.l','HG-right':'ZA:Hyoglossus muscle.r'}

def require(condition,message):
    if not condition:raise ValueError(message)
def read(path):return json.loads(Path(path).read_text())
def timer(seconds=48.):
    end=time.monotonic()+seconds
    def check():
        if time.monotonic()>end:raise TimeoutError('Extrinsic fairing stage reached its finite worker deadline')
    return check
def checked_plan(output):
    p=output/'plan.json';plan=read(p)
    require(plan.get('schema')=='tongue-single-constrained-fairing-plan/v1','Wrong fairing plan')
    for path,sha in plan['bindings'].items():require(digest_file(path)==sha,'A planned input changed: '+path)
    return plan
def canonical(tri):return min(tuple(map(tuple,np.roll(tri,i,axis=0)))for i in range(3))
def compact(mesh):
    v,first,alias=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True)
    return v,alias[mesh.faces]

def quality(v,f,normals=None):
    t=v[f];raw=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);length=np.linalg.norm(raw,axis=1)
    require(np.all(length>0),'Quality check found a zero-area face')
    fn=raw/length[:,None];edges=defaultdict(list)
    for i,tri in enumerate(f):
        for x,y in zip(tri,np.roll(tri,-1)):edges[tuple(sorted((int(x),int(y))))].append(i)
    pairs=np.array([r for r in edges.values()if len(r)==2]);angles=np.rad2deg(np.arccos(np.clip(np.einsum('ij,ij->i',fn[pairs[:,0]],fn[pairs[:,1]]),-1,1)))
    el=np.linalg.norm(v[np.array(list(edges))[:,0]]-v[np.array(list(edges))[:,1]],axis=1)
    root=np.arange(len(f))
    def find(i):
        while root[i]!=i:root[i]=root[root[i]];i=root[i]
        return i
    for a,b in pairs[angles<.1]:root[find(a)]=find(b)
    regions=defaultdict(list)
    for i in range(len(f)):regions[find(i)].append(i)
    largest=sorted(regions.values(),key=lambda r:length[r].sum(),reverse=True)[:5]
    result=dict(triangles=len(f),vertices=len(v),surface_area_mm2=float(length.sum()/2*1e6),
        volume_ml=float(np.einsum('ij,ij->i',t[:,0]-v.mean(0),np.cross(t[:,1]-v.mean(0),t[:,2]-v.mean(0))).sum()/6*1e6),
        largest_triangle_mm2=float(length.max()/2*1e6),edge_length_m_quantiles=np.quantile(el,[.5,.95,1]).tolist(),
        geometric_dihedral_degrees_quantiles=np.quantile(angles,[.5,.9,.99,1]).tolist(),
        geometric_edges_over45=int((angles>45).sum()),geometric_edges_over150=int((angles>150).sum()),
        largest_coplanar_patches=[dict(faces=len(r),area_mm2=float(length[r].sum()/2*1e6),
            span_m=np.ptp(t[r].reshape(-1,3),axis=0).tolist())for r in largest])
    if normals is not None:
        nn=normals.astype(float);nn/=np.linalg.norm(nn,axis=2)[:,:,None]
        an=np.rad2deg(np.arccos(np.clip(np.einsum('fci,fi->fc',nn,fn),-1,1)))
        result.update(corner_normal_angle_degrees_quantiles=np.quantile(an,[.5,.9,.99,1]).tolist(),
                      backward_corner_normals=int((an>90).sum()))
    return result

def plan_stage(output):
    require(not output.exists()and output.is_relative_to(WORK),'Use one new local fairing directory')
    selected=WORK/'hyoglossus-origin-v7-nonplanar-rim-topology/selected-extrinsics.json'
    rows=read(selected)['parts'];source_audit=WORK/'source-v1/source-audit.json';context=read(source_audit)
    bones=[r for r in context['contextParts']if r['sourcePart'].get('system')=='skeletal']
    require(len(bones)==8,'All eight fixed skeletal sources are required')
    covers=WORK/'shared-coverings-exposed-apo-sl-v3/packed'
    bindings={str(selected):digest_file(selected),str(source_audit):digest_file(source_audit),str(Path(__file__).resolve()):digest_file(__file__)}
    sources={}
    for name,pid in PARTS.items():
        row=next(r for r in rows if r['originalPartId']==pid);p=Path(row['mesh'])
        require(digest_file(p)==row['sha256'],'Selected producer changed')
        bindings[str(p)]=row['sha256'];sources[name]=dict(part_id=pid,path=str(p),sha256=row['sha256'],
            fixed_superior_face_minimum_y=1.630 if name.startswith('GG') else 1.621,
            maximum_displacement_m=.002 if name.startswith('GG')else.004,
            maximum_absolute_volume_fraction_change=.15 if name.startswith('GG')else.25)
    fixed=[]
    for row in bones:
        p=Path(row['meshSavePath']);bindings[str(p)]=digest_file(p);fixed.append(dict(part_id=row['partId'],path=str(p),sha256=bindings[str(p)]))
    coverings=[]
    for key in ('mucosa','aponeurosis'):
        p=covers/('TA98TONGUE_male_'+key+'.npz');bindings[str(p)]=digest_file(p);coverings.append(dict(key=key,path=str(p),sha256=bindings[str(p)]))
    inventory=ROOT/'reports/ta98-tongue-inventory.json';bindings[str(inventory)]=digest_file(inventory)
    require(read(inventory)['term_count']==36,'Complete classified inventory is required')
    plan=dict(schema='tongue-single-constrained-fairing-plan/v1',sources=sources,bones=fixed,coverings=coverings,
        bindings=bindings,inventory=str(inventory),method='Single discrete biharmonic lower-surface fairing with fixed literal insertion faces and cortical controls',
        fidelity_weight=.02,refine_free_edges_longer_than_m=.001,protected_mucosal_distance_m=.00025,
        fixed_cortical_contact_distance_m=1e-6,maximum_triangles_per_muscle=12000,
        anatomical_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
        adjoining_interfaces=['Literal superior GG/HG insertion charts remain fixed','Original cortical attachment sites stay fixed at reviewed adjacent-unjoined gap',
            'Every bone and both actual finite covering cells are immutable obstacles','Muscle envelopes may interdigitate; they are not exclusive histological cells'],
        one_parameter_set_no_search=True,source_geometry_modified=False,installed=False,anatomical_acceptance=False)
    output.mkdir(parents=True);atomic_json(output/'plan.json',plan);return dict(output=str(output),stage='planned',parts=4)

def fit_stage(output,name):
    began=time.monotonic();check=timer();plan=checked_plan(output);spec=plan['sources'][name];dest=output/(name+'.npz')
    require(not dest.exists(),'The one candidate is immutable; do not rerun an existing fit')
    original=Mesh.load(spec['path']);v,f=compact(original);raw_words=original.attributes['normals'][original.faces].copy()
    signs=original.attributes.get('normalSigns',np.ones(len(original.vertices),dtype='i1'))[original.faces].copy()
    bone_meshes=[Mesh.load(r['path'])for r in plan['bones']];bone_tree=BVH(Mesh(np.concatenate([m.vertices for m in bone_meshes]),
        np.concatenate([m.faces+sum(len(q.vertices)for q in bone_meshes[:i])for i,m in enumerate(bone_meshes)]),FRAME))
    near,_=bone_tree.closest_points(v,backend='open3d');bone_control=near['distance_m']<=plan['fixed_cortical_contact_distance_m']
    mucosa=Mesh.load(plan['coverings'][0]['path']);muco_near,_=BVH(mucosa).closest_points(v,backend='open3d')
    superior=np.all(v[f][:,:,1]>=spec['fixed_superior_face_minimum_y'],axis=1)
    insertion=np.any(muco_near['distance_m'][f]<plan['protected_mucosal_distance_m'],axis=1)
    cortical=np.all(bone_control[f],axis=1);protected=superior|insertion|cortical
    fixed=np.zeros(len(v),bool);fixed[np.unique(f[protected])]=True;fixed|=bone_control
    edge_faces=defaultdict(list)
    for face,tri in enumerate(f):
        for x,y in zip(tri,np.roll(tri,-1)):edge_faces[tuple(sorted((int(x),int(y))))].append(face)
    points=v.tolist();constraints=fixed.tolist();edge_mid={}
    for (x,y),ii in edge_faces.items():
        if len(ii)==2 and not np.any(protected[ii])and np.linalg.norm(v[x]-v[y])>plan['refine_free_edges_longer_than_m']:
            edge_mid[(x,y)]=len(points);points.append(((v[x]+v[y])/2).tolist());constraints.append(False)
    faces=f[protected].tolist();parent=np.flatnonzero(protected).tolist();fixed_count=len(faces)
    for i,tri in enumerate(f):
        if protected[i]:continue
        loop=[]
        for x,y in zip(tri,np.roll(tri,-1)):
            loop.append(int(x));mid=edge_mid.get(tuple(sorted((int(x),int(y)))))
            if mid is not None:loop.append(mid)
        center=len(points);points.append(v[tri].mean(0).tolist());constraints.append(False)
        for x,y in zip(loop,np.roll(loop,-1)):faces.append([center,int(x),int(y)]);parent.append(i)
    v0=np.asarray(points);nf=np.asarray(faces,int);fixed=np.asarray(constraints,bool);free=np.flatnonzero(~fixed)
    require(len(nf)<=plan['maximum_triangles_per_muscle']and len(free)>0,'Refinement resource cap or free surface exhausted')
    e=np.concatenate([nf[:,[0,1]],nf[:,[1,2]],nf[:,[2,0]]]);e=np.unique(np.sort(e,axis=1),axis=0)
    adjacency=sparse.coo_matrix((np.ones(len(e)*2),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(v0),len(v0))).tocsr()
    degree=np.asarray(adjacency.sum(1)).ravel();L=sparse.eye(len(v0))-sparse.diags(1/degree)@adjacency
    # Hard boundary coordinates are eliminated, not assigned a soft penalty.
    B=(L.T@L).tocsr()+plan['fidelity_weight']*sparse.eye(len(v0));locked=np.flatnonzero(fixed)
    rhs=plan['fidelity_weight']*v0[free]-B[free][:,locked]@v0[locked]
    check();solved=spsolve(B[free][:,free].tocsc(),rhs);delta=solved-v0[free]
    max_raw=float(np.linalg.norm(delta,axis=1).max());fraction=min(1.,spec['maximum_displacement_m']/max_raw)
    actual=v0.copy();actual[free]+=fraction*delta
    stored=(actual-ORIGIN).astype('<f4');actual=stored.astype(float)+ORIGIN
    require(np.array_equal(actual[locked],v0[locked]),'Fixed emitted coordinates changed on encoding')
    tr=actual[nf];fn=np.cross(tr[:,1]-tr[:,0],tr[:,2]-tr[:,0]);smooth=np.zeros_like(actual)
    for c in range(3):np.add.at(smooth,nf[:,c],fn)
    nn=np.linalg.norm(smooth,axis=1);require(np.all(nn>0),'Fairing created a zero shading normal')
    smooth/=nn[:,None];words=np.rint(smooth[nf]*32767).astype('<i2');corner_sign=np.ones((len(nf),3),dtype='i1')
    words[:fixed_count]=raw_words[protected];corner_sign[:fixed_count]=signs[protected]
    mesh=Mesh(actual[nf].reshape(-1,3),np.arange(nf.size).reshape(-1,3),FRAME,
        dict(normals=words.reshape(-1,3),normalSigns=corner_sign.reshape(-1)),
        dict(source_part_id=spec['part_id'],source_geometry_sha256=original.geometry_sha,proposal='Single hard-boundary biharmonic physical fairing',anatomical_acceptance=False))
    check();mesh.save(dest);atomic_npz(output/(name+'-stored.npz'),positions=stored[nf].reshape(-1,3),faces=mesh.faces,
        geometryOrigin=ORIGIN,attribute_normals=mesh.attributes['normals'],attribute_normalSigns=mesh.attributes['normalSigns'],frame=np.array(FRAME))
    before=quality(v,f,raw_words*signs[:,:,None]);after=quality(actual,nf,words*corner_sign[:,:,None])
    keys={canonical(t):i for i,t in enumerate(mesh.triangles)};literal=all(canonical(t)in keys for t in original.triangles[protected])
    volume_change=after['volume_ml']/before['volume_ml']-1
    proof=dict(schema='tongue-single-constrained-extrinsic-fit/v1',part=name,source=spec['path'],source_sha256=spec['sha256'],
        candidate=str(dest),candidate_sha256=digest_file(dest),source_plan_sha256=digest_file(output/'plan.json'),before=before,after=after,
        fixed_superior_faces=int(superior.sum()),fixed_mucosal_insertion_faces=int(insertion.sum()),fixed_cortical_faces=int(cortical.sum()),
        all_fixed_faces=int(protected.sum()),fixed_cortical_controls=int(bone_control.sum()),all_fixed_triangles_literal=literal,
        raw_protected_corner_words_and_signs_literal=bool(np.array_equal(words[:fixed_count],raw_words[protected])and np.array_equal(corner_sign[:fixed_count],signs[protected])),
        fixed_contact_controls_literal=bool(np.array_equal(actual[:len(v)][bone_control],v[bone_control])),
        free_controls=len(free),maximum_unconstrained_displacement_m=max_raw,maximum_applied_displacement_m=float(np.linalg.norm(actual-v0,axis=1).max()),
        common_bounded_solve_fraction=fraction,volume_fraction_change=volume_change,
        volume_retention_pass=abs(volume_change)<=spec['maximum_absolute_volume_fraction_change'],
        physically_changed_controls=int((np.linalg.norm(actual-v0,axis=1)>1e-8).sum()),
        source_geometry_modified=False,skeleton_modified=False,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/(name+'-fit.json'),proof);checked_plan(output)
    return dict(part=name,triangles=len(nf),fixed_faces=int(protected.sum()),moved=proof['physically_changed_controls'],
        max_displacement_mm=proof['maximum_applied_displacement_m']*1000,volume_change=volume_change,
        volume_pass=proof['volume_retention_pass'],elapsed_seconds=proof['elapsed_seconds'])

def audit_stage(output,name):
    began=time.monotonic();check=timer();plan=checked_plan(output);fit=read(output/(name+'-fit.json'));dest=output/(name+'-audit.json')
    require(not dest.exists()and digest_file(fit['candidate'])==fit['candidate_sha256'],'Use the one unchanged fit for one saved audit')
    m=Mesh.load(fit['candidate']);checks=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000)
    contacts=[]
    for row in plan['bones']+plan['coverings']:
        check();other=Mesh.load(row['path']);arrays,summary=pair_contacts(m,other,0.,max_pairs=250000,max_tests=2000000)
        key=row.get('part_id',row.get('key'));path=output/(name+'-'+key.replace(':','_')+'-contacts.npz');atomic_npz(path,**arrays)
        contacts.append(dict(counterpart=key,path=str(path),sha256=digest_file(path),summary=summary))
    topo=checks['topology'];passing=checks['passed']and len(topo['components'])==1 and topo['components'][0]['closed_orientable_genus']==0
    passing=passing and fit['volume_retention_pass']and fit['all_fixed_triangles_literal']and fit['raw_protected_corner_words_and_signs_literal']and fit['fixed_contact_controls_literal']
    passing=passing and all(r['summary']['complete_within_pair_budget']and not r['summary']['truncated']and r['summary']['crossing_face_pairs']==0 and r['summary']['coplanar_overlap_face_pairs']==0 for r in contacts)
    result=dict(schema='tongue-single-constrained-extrinsic-audit/v1',part=name,passed=bool(passing),candidate=fit['candidate'],candidate_sha256=fit['candidate_sha256'],fit_sha256=digest_file(output/(name+'-fit.json')),
        actual_geometry_audit=checks,actual_contacts=contacts,quality=fit['after'],anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    check();checked_plan(output);atomic_json(dest,result)
    return dict(part=name,passed=bool(passing),self=checks['passed'],contacts={r['counterpart']:r['summary']['crossing_face_pairs']for r in contacts},elapsed_seconds=result['elapsed_seconds'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['plan','fit','audit'],required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--part',choices=PARTS)
    a=p.parse_args();o=a.output.resolve();require(o.is_relative_to(WORK),'Local tongue workbench only')
    if a.stage=='plan':require(a.part is None,'Plan is shared by four muscles');result=plan_stage(o)
    else:require(a.part is not None,'Explicit part required');result=(fit_stage if a.stage=='fit'else audit_stage)(o,a.part)
    print(json.dumps(result))
