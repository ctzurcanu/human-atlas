"""Finite collars built around the actual joined urinary passage.

The fluid collar is independently closed and audited. Its conservative finite
offset supplies the external wall; source body/ureter envelopes are then united
before subtracting the one connected fluid domain. No solid lumen plug is used.
"""
from pathlib import Path
import json

import manifold3d as md
import numpy as np

from .algorithms.annular_pipe_envelopes import boundary_loops
from .algorithms.finite_tissue_boolean import make_solid, provenance
from .audit import audit
from .bladder_carriers import recover
from .bladder_workbench import crop
from .common import atomic_json, atomic_npz, digest_file
from .db import Database
from .model import Mesh
from .renal_master import hull_solid
from .review import prepare


def arrays(solid):
    m=solid.to_mesh64()
    return np.asarray(m.vert_properties)[:,:3].copy(),np.asarray(m.tri_verts,int).copy()


def finite_offset(solid,minimum_m):
    unit=md.Manifold.sphere(1.,12);v,f=arrays(unit);t=v[f]
    normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    inradius=float(np.min(abs(np.einsum('ij,ij->i',normal,t[:,0]))))
    if inradius<=0:raise ValueError('The conservative offset kernel is invalid')
    kernel=unit.scale([minimum_m/inradius]*3)
    return solid.minkowski_sum(kernel),dict(minimum_offset_m=minimum_m,kernel_faces=len(f),kernel_inradius_before_scale=inradius,
        kernel_circumradius_m=float(np.linalg.norm(v,axis=1).max()*minimum_m/inradius),
        scope='Finite construction clearance, not a measured histological coat thickness')


def build(join_dir,source_carriers,bladder_wall,output,*,source_cache):
    join_dir,output=Path(join_dir).resolve(),Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    joined=Mesh.load(join_dir/'urinary-wall.npz');bladder=Mesh.load(bladder_wall)
    pose_report=json.loads((join_dir/'junctions.json').read_text())
    if pose_report['bladder_source_sha256']!=bladder.geometry_sha:raise ValueError('The bladder fit changed')
    sources=Path(source_carriers).resolve()
    with np.load(sources,allow_pickle=False)as z:
        if not np.array_equal(z['wall_faces'],bladder.faces):raise ValueError('The bladder source face binding changed')
        bv=bladder.vertices;bo=z['outer_faces'].copy();bf=z['fluid_faces'].copy();bd=z['computational_disks'].copy()
    outer_parts=[];fluid_parts=[];checks={};rows=[]
    def source(name,v,f):
        mesh=Mesh(v,f,joined.frame);check=audit(mesh,require_closed=True,expected_components=1,expected_genus=[0],max_pairs=1000000)
        checks[name]=check
        if not check['passed']:
            mesh.save(output/(name+'-failure.npz'));atomic_json(output/(name+'-failure.json'),check)
            raise ValueError('A collar classifier fails actual geometry checks: '+name)
        solid,_=make_solid(v,f);return solid
    outer_parts.append(source('bladder-outer',bv,bo));fluid_parts.append(source('bladder-fluid',bv,bf))
    for side in ['left','right']:
        # The transport is source-bound in the preceding construction. Reuse
        # only independently passed physical tubes, never its rejected outer loft.
        row=next(q for q in pose_report['ureters']if q['side']==side)
        if not row['transported_wall_check']['passed'] or not row['fixed_tail_literal']:
            raise ValueError('A rejected or root-moving ureter pose cannot be used')
        moved=Mesh.load(join_dir/(side+'-transported-wall.npz'))
        cache=Path(source_cache).resolve()
        lineage=json.loads((cache/(side+'-sheet-lineage.json')).read_text())
        if not lineage.get('explicit_root_and_mucosal_lineage')or lineage['source_geometry_sha256']!=row['source_geometry_sha256']:
            raise ValueError('The complete, source-bound native lining is required')
        with np.load(cache/(side+'-carriers/carriers.npz'),allow_pickle=False)as z:
            uv=moved.vertices;uo=z['outer_faces'].copy();uf=z['fluid_faces'].copy();ud=z['computational_disks'].copy();local_inner=z['distal_ureter_ring'].copy()
        with np.load(join_dir/(side+'-transport.npz'),allow_pickle=False)as z:
            inner_bridge=z['inner_bridge'].copy();body_inner=z['bladder_inner_ring'].copy();global_inner=z['ureter_inner_ring'].copy();tube_outer=z['ureter_outer_ring'].copy()
        offset=int(global_inner.min()-local_inner.min())
        if set(global_inner-offset)!=set(local_inner)or not np.array_equal(joined.vertices[global_inner],uv[global_inner-offset]):
            raise ValueError('The collar/ureter coordinates disagree')
        bcap=bd[np.all(np.isin(bd,body_inner),axis=1)]
        ucap=ud[np.all(np.isin(ud,local_inner),axis=1)]+offset
        if not len(bcap)or not len(ucap):raise ValueError('Explicit paired fluid mouth caps are missing')
        # Overlap the computational end into the receiving fluid, rather than
        # closing an oblique nonplanar loft on its own mouth. The full cap is
        # placed beyond the loft's bounding plane; it never becomes anatomy.
        old_cap=bcap[:,[0,2,1]];cap_points=joined.vertices[np.unique(old_cap)]
        cq=joined.vertices[old_cap];direction=np.cross(cq[:,1]-cq[:,0],cq[:,2]-cq[:,0]).sum(0);direction/=np.linalg.norm(direction)
        loft_points=joined.vertices[np.unique(inner_bridge)]
        extension=max(0.,float((loft_points@direction).max()-(cap_points@direction).min()))+.00005
        cap_loop=boundary_loops(old_cap)[0]
        new_ids=np.arange(len(joined.vertices),len(joined.vertices)+len(cap_loop));mapping=dict(zip(map(int,cap_loop),map(int,new_ids)))
        collar_positions=np.r_[joined.vertices,joined.vertices[cap_loop]+extension*direction]
        sides=[]
        for a,b in zip(cap_loop,np.roll(cap_loop,-1)):
            sides.extend([[int(a),int(b),mapping[int(b)]],[int(a),mapping[int(b)],mapping[int(a)]]])
        new_cap=np.array([[mapping[int(i)]for i in t]for t in old_cap],int)
        collar_faces=np.r_[inner_bridge[:,[0,2,1]],np.asarray(sides,int),new_cap,ucap[:,[0,2,1]]]
        # The retained right ostial polygon contains a small re-entrant source
        # notch. A ruled strip folds across that notch. Regularize this short
        # opening locally as the minimum convex passage containing both actual
        # profiles; the complete source walls elsewhere remain untouched.
        inward_ureter=np.asarray(row['outward_normal'])
        support=np.r_[joined.vertices[np.unique(inner_bridge)],
                      collar_positions[new_ids],
                      joined.vertices[global_inner]+.0003*inward_ureter]
        collar_hull=hull_solid(support);cv,cf=arrays(collar_hull)
        collar=source(side+'-fluid-collar',cv,cf)
        tube_o=source(side+'-ureter-outer',uv,uo);tube_f=source(side+'-ureter-fluid',uv,uf)
        outer_parts.append(tube_o);fluid_parts.extend([tube_f,collar])
        center=joined.vertices[global_inner].mean(0)
        outer_radius=float(np.linalg.norm(joined.vertices[tube_outer]-center,axis=1).max())
        inner_radius=float(np.linalg.norm(joined.vertices[global_inner]-center,axis=1).min())
        thickness=max(.00035,1.15*(outer_radius-inner_radius))
        sleeve,offset_evidence=finite_offset(collar,thickness);sv,sf=arrays(sleeve)
        outer_parts.append(source(side+'-finite-collar',sv,sf))
        Mesh(sv,sf,joined.frame).save(output/(side+'-collar-outer.npz'))
        rows.append(dict(side=side,collar=offset_evidence,computational_cap_extension_m=extension,
                         fluid_collar_construction='minimum convex support of both native profiles with computational overlaps',
                         local_reentrant_ostial_notches_regularized=True,source_ring_positions_used_literally=True,
                         source_renal_root_fixed=True,source_pose_sha256=digest_file(join_dir/(side+'-transport.npz'))))
    outer=outer_parts[0]
    for p in outer_parts[1:]:outer=outer+p
    fluid=fluid_parts[0]
    for p in fluid_parts[1:]:fluid=fluid+p
    ov,of=arrays(outer);fv,ff=arrays(fluid)
    # Rebind the two complete classifiers once. The wall then retains exact
    # face lineage to fluid versus external boundary for the coat constructor.
    outer,_=make_solid(ov,of);fluid,fid=make_solid(fv,ff)
    wall=outer-fluid;v,f=arrays(wall);wm=wall.to_mesh64();oids,parent=provenance(wm)
    inner_ids=np.flatnonzero(oids==fid)
    whole=Mesh(v,f,joined.frame);whole.save(output/'urinary-wall.npz')
    Mesh(ov,of,joined.frame).save(output/'outer.npz');Mesh(fv,ff,joined.frame).save(output/'fluid.npz')
    atomic_npz(output/'inner-face-ids.npz',inner_face_ids=inner_ids)
    final=audit(whole,require_closed=True,expected_components=1,expected_genus=[2],max_pairs=1000000)
    fluid_check=audit(Mesh(fv,ff,joined.frame),require_closed=True,expected_components=1,expected_genus=[0],max_pairs=1000000)
    outer_check=audit(Mesh(ov,of,joined.frame),require_closed=True,expected_components=1,expected_genus=[0],max_pairs=1000000)
    escape=float((fluid-outer).volume())
    scope=max(float(np.ptp(v,axis=0).max()),np.finfo(float).tiny)
    numerical=2048*np.finfo(float).eps*scope**3*(len(ff)+len(of))
    result=dict(schema='urinary-finite-offset-collars/v1',source_join_sha256=digest_file(join_dir/'junctions.json'),
        sources=checks,collars=rows,whole_check=final,fluid_check=fluid_check,outer_check=outer_check,
        fluid_outside_outer_m3=escape,numerical_volume_bound_m3=numerical,
        passed=bool(final['passed']and fluid_check['passed']and outer_check['passed']and escape<=numerical),
        rejected_outer_strips_used=False,source_renal_roots_fixed=True,skeleton_modified=False,
        installed=False,anatomical_acceptance=False)
    atomic_json(output/'junctions.json',result)
    if result['passed']:
        result['carriers']=recover(whole,output/'verified-carriers',inner_face_ids=inner_ids,
            mouth_labels=['urethral-neck','right-renal-root','left-renal-root'],cap_forbid_source_diagonals=True)
        cropped=crop(whole,bladder.bounds+[[-.01,-.01,-.01],[.01,.04,.01]])
        cropped.save(output/'junction-review.npz')
        center=bladder.bounds.mean(0)
        result['review']=prepare(Database(),dict(frame=whole.frame,
            meshes=[dict(mesh=str(output/'junction-review.npz'),name='Continuous urinary wall; uncapped upper crop',color=[.78,.43,.4,1])],
            views=[dict(name='posterior-joins',direction=[0,.2,-1]),dict(name='front-joins',direction=[0,.2,1]),
                   dict(name='bladder-and-joins-section',direction=[0,.1,1],section=dict(normal=[0,0,1],point=center.tolist(),clear_outer=True))],resolution=1000))
        atomic_json(output/'review-request.json',result)
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('join_dir');p.add_argument('--source-carriers',required=True)
    p.add_argument('--bladder-wall',required=True);p.add_argument('--output',required=True);p.add_argument('--source-cache',required=True)
    a=p.parse_args();r=build(a.join_dir,a.source_carriers,a.bladder_wall,a.output,source_cache=a.source_cache)
    print(json.dumps(dict(passed=r['passed'],failures=r['whole_check']['failures'],review=r.get('review')),indent=2))
