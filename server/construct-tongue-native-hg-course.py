"""One original native-HG course correction; no posterior origin fan.

The native curved belly receives the existing whole-tongue basal cage, an
explicit upper lateral contraction and two short C2 origin-site fields. Only
the fixed hyoid occupied volume is subtracted. This is a numerical proposal,
not a measured attachment or accepted anatomical reconstruction.
"""
from __future__ import annotations
import argparse,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit,topology
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.longitudinal_organ_cage import LongitudinalOrganCage
from geometry_tools.algorithms.finite_tissue_boolean import make_solid

W=ROOT/'.local-models/geometry-workbench/tongue-completion';ORIGIN=np.array([0.,1.625,0.])
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def module(name,file):
    s=importlib.util.spec_from_file_location(name,ROOT/'server'/file);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def fall(y,a,b):
    u=np.clip((y-a)/(b-a),0,1);return 1-u*u*(3-2*u)
def kernel(points,centers,radius):
    r=np.linalg.norm(points[:,None,:]-centers[None,:,:],axis=2)/radius;q=np.maximum(1-r,0)
    return q**4*(4*r+1)
def ordinary_normals(v,f):
    t=v[f];c=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);n=np.zeros_like(v)
    for j in range(3):np.add.at(n,f[:,j],c)
    n/=np.maximum(np.linalg.norm(n,axis=1)[:,None],1e-30);return n

def main(output,side):
    began=time.monotonic();end=began+48
    def check():
        if time.monotonic()>end:raise TimeoutError('Native HG correction exceeded its finite worker budget')
    require(not output.exists()and output.is_relative_to(W),'Use a new local source-based candidate directory')
    output.mkdir(parents=True)
    src=W/'source-v1/meshes'/('native-03.npz'if side=='left'else'native-04.npz');native=Mesh.load(src)
    cps=W/'basal-whole-protected-v1/fit-input.json';cage=LongitudinalOrganCage(**read(cps)['cage'])
    row=read(W/'source-v1/source-audit.json');bones=[r for r in row['contextParts']if r['sourcePart'].get('system')=='skeletal'];bm={r['partId']:Mesh.load(r['meshSavePath'])for r in bones}
    hyoid=[m for k,m in bm.items()if k.startswith('TA98HYOID:')];bv=np.concatenate([m.vertices for m in hyoid]);bf=np.concatenate([m.faces+sum(len(q.vertices)for q in hyoid[:i])for i,m in enumerate(hyoid)]);bone=Mesh(bv,bf,native.frame)
    body=bm['TA98HYOID:A02.1.16.002:whole'];horn=bm['TA98HYOID:A02.1.16.004:'+('l'if side=='left'else'r')]
    p=native.vertices;first=cage.map(p[:,[2,0,1]])[:,[1,2,0]]
    # Source-site witnesses, not an invented30mm posterior footprint.
    native_anchor=[int(np.argmin(p[:,1])),int(np.argmin(p[:,2]))];centers=first[native_anchor];targets=[];witnesses=[]
    for ordinal,host in enumerate([body,horn]):
        t=host.triangles;n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);n/=np.linalg.norm(n,axis=1)[:,None];c=t.mean(1);ids=np.flatnonzero(n[:,1]>.5)
        face=int(ids[np.argmin(np.linalg.norm(c[ids]-p[native_anchor[ordinal]],axis=1))]);target=c[face]+.00015*n[face]
        targets.append(target);witnesses.append(dict(native_control=native_anchor[ordinal],native_point=p[native_anchor[ordinal]].tolist(),host_part_id=host.source['record']['id'],host_face=face,target_point=target.tolist(),cortical_face_center=c[face].tolist(),initial_tissue_side_clearance_m=.00015))
    targets=np.array(targets);radius=.009;weights=np.linalg.solve(kernel(centers,centers,radius),targets-centers)
    def deformation(points):
        q=points.copy();q+=kernel(points,centers,radius)@weights
        # A smooth contraction only where actual source feasibility reported
        # mucosal crossings, reaching20% inside the upper lateral strip.
        blend=1-fall(points[:,1],1.618,1.625);q[:,0]*=1-.20*blend
        return q
    proposed=deformation(first)
    jac=np.empty((len(p),3,3));step=1e-6
    for axis in range(3):
        d=np.eye(3)[axis]*step;jac[:,:,axis]=(deformation(first+d)-deformation(first-d))/(2*step)
    determinant=np.linalg.det(jac);maxmove=np.linalg.norm(proposed-first,axis=1).max()
    plan=dict(schema='native-hg-curved-source-short-origin-plan/v1',side=side,source=str(src),source_sha256=digest_file(src),source_geometry_sha256=native.geometry_sha,
        original_source_faces=len(native.faces),source_whole_cage=str(cps),source_whole_cage_sha256=digest_file(cps),origin_field_radius_m=radius,
        origin_field='Two fixed C2 Wendland displacement sites with actual superior cortical-facet witnesses',witnesses=witnesses,
        upper_insertion='Smooth20% lateral contraction at Y1.618–1.625m; old artificial insertion charts deliberately replaced',
        origin_sites='Adjacent anterior hyoid-body and greater-horn cortex only; earlier30mm spanning footprint rejected by measured29.6mm source mismatch',
        discrete_jacobian_determinant_range=[float(determinant.min()),float(determinant.max())],maximum_local_origin_and_insertion_motion_m=float(maxmove),
        microscopic_entheses_established=False,source_geometry_modified=False,skeleton_modified=False,anatomical_acceptance=False,installed=False,
        bindings={str(src):digest_file(src),str(cps):digest_file(cps),str(W/'source-v1/source-audit.json'):digest_file(W/'source-v1/source-audit.json'),str(Path(__file__).resolve()):digest_file(__file__)})
    atomic_json(output/'plan.json',plan)
    require(determinant.min()>.05 and maxmove<=.007,'Source-site field reverses or exceeds explicit7mm local displacement bound')
    # Carry source shading through both continuous maps before the small bone
    # trim; new cut surfaces use actual geometric normals.
    cagej=cage.jacobians(p[:,[2,0,1]]);oldn=native.attributes['normals'].astype(float)[:,[2,0,1]]
    n=np.linalg.solve(cagej.transpose(0,2,1),oldn[:,:,None])[:,:,0][:,[1,2,0]];n=np.linalg.solve(jac.transpose(0,2,1),n[:,:,None])[:,:,0];n/=np.linalg.norm(n,axis=1)[:,None]
    attrs=dict(normals=np.rint(np.clip(n,-1,1)*32767).astype('<i2'));deformed=Mesh(proposed,native.faces,native.frame,attrs)
    deformed.save(output/'native-course-before-short-bone-trim.npz');check()
    host,oid=make_solid(proposed,native.faces,properties=attrs['normals'][native.faces]);occupied,_=make_solid(bv,bf);result=host-occupied
    decode=module('native_hg_decode','fit-tongue-mucosal-insertions.py');before,ancestry,_=decode.decoded_candidate(result,deformed,oid,native.frame)
    raw=result.to_mesh64();tri=np.asarray(raw.vert_properties,float)[np.asarray(raw.tri_verts,int),:3];before=Mesh(tri.reshape(-1,3),before.faces,native.frame,before.attributes,before.source)
    enc=module('native_hg_bone_encoding','repair-tongue-hyoglossus-origin.py');candidate,stored,adjacency=enc.encode_bone_adjacency(before,bone,BVH(bone));path=output/('HG-'+side+'.npz');candidate.save(path)
    atomic_npz(output/('HG-'+side+'-stored.npz'),positions=stored,faces=candidate.faces,geometryOrigin=ORIGIN,frame=np.array(candidate.frame),attribute_normals=candidate.attributes['normals'])
    check();checks=audit(candidate,require_closed=True,max_pairs=250000,max_tests=2000000);contacts=[]
    coverdir=W/'shared-coverings-exposed-apo-sl-v3/packed';others={**bm,**{k:Mesh.load(coverdir/('TA98TONGUE_male_'+k+'.npz'))for k in ['mucosa','aponeurosis']}}
    for key,other in others.items():
        check();arrays,summary=pair_contacts(candidate,other,0.,max_pairs=250000,max_tests=2000000);dest=output/(key.replace(':','_')+'-contacts.npz');atomic_npz(dest,**arrays);contacts.append(dict(counterpart=key,summary=summary,witnesses=str(dest),witness_sha256=digest_file(dest)))
    qs=module('native_hg_shape_measure','fair-tongue-extrinsic-envelopes.py');v,f=qs.compact(candidate);words=candidate.attributes['normals'][candidate.faces]
    beforeq=qs.quality(native.vertices,native.faces,native.attributes['normals'][native.faces]);afterq=qs.quality(v,f,words)
    require(all(digest_file(p)==s for p,s in plan['bindings'].items()),'An input changed during construction')
    passed=checks['passed']and all(r['summary']['complete_within_pair_budget']and not r['summary']['truncated']and r['summary']['crossing_face_pairs']==0 and r['summary']['coplanar_overlap_face_pairs']==0 for r in contacts)
    report=dict(schema='native-hg-curved-source-short-origin-candidate/v1',side=side,passed=bool(passed),candidate=str(path),candidate_sha256=digest_file(path),
        plan_sha256=digest_file(output/'plan.json'),source_triangle_ancestry=ancestry,adjacency=adjacency,actual_geometry_audit=checks,contacts=contacts,
        native_shape=beforeq,actual_stored_shape=afterq,source_geometry_modified=False,skeleton_modified=False,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'candidate.json',report)
    spec=dict(frame=native.frame,resolution=800,wall_seconds=45,meshes=[dict(mesh=str(path),name='Native-based HG '+side,color=[.74,.32,.34,1])]+[dict(mesh=r['meshSavePath'],name='Fixed '+r['partId'],color=[.65,.62,.53,.3])for r in bones if r['partId'].startswith('TA98HYOID:')],views=[dict(name='front',direction=[0,0,1]),dict(name='right',direction=[1,0,0]),dict(name='sagittal',direction=[1,0,0],section=dict(normal=[1,0,0],point=[.010 if side=='left'else-.010,1.620,.045],clear_outer=True))])
    atomic_json(output/'review-spec.json',spec)
    print(json.dumps(dict(output=str(output),passed=bool(passed),self=checks['passed'],contacts={r['counterpart']:r['summary']['crossing_face_pairs']for r in contacts},triangles=len(candidate.faces),elapsed_seconds=report['elapsed_seconds'])))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--side',choices=['left','right'],required=True);a=p.parse_args();main(a.output.resolve(),a.side)
