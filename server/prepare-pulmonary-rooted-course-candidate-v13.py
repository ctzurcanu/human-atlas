"""One source-bound pulmonary course initialization; never an install gate.

The declared field moves source section centres and rotates each literal ring
properly. It does not create an inner wall, open a collector, or establish a
communicating vascular junction. Existing finite lingular cells move as one
rigid block, including their actual main-arterial cut profile.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_station_topology import recover_station_rings, StationTopologyError

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OUT=BASE/'rooted-course-initialization-v13'
FRAME='model:male:assembled-metres'
PV=[f'ZA:{side} {level} pulmonary vein'for side in ['Left','Right']for level in ['inferior','superior']]
FIXED=PV+['ZA:Pulmonary trunk','ZA:Bifurcation of pulmonary trunk']
LINGULAR=['ZA:Superior lingular artery of left lung','ZA:Inferior lingular artery of left lung']


def read_sources():
    family=json.loads((BASE/'family-contract-v2.json').read_text())
    path=Path(family['source_catalogue']);assert digest_file(path)==family['source_catalogue_sha256']
    atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']}
    previous=json.loads((BASE/'takeoff-and-root-controls-v10/evidence.json').read_text())
    ids=previous['whole_runtime_qualified_blood_source_ids'];assert len(ids)==54
    reader=AtlasReader();meshes={pid:reader.read_part(atlas,parts[pid],frame=FRAME)for pid in ids}
    return family,previous,atlas,parts,reader,meshes


def contract():
    began=time.monotonic();family,previous,atlas,parts,reader,meshes=read_sources()
    OUT.mkdir()
    inp=np.load(BASE/'takeoff-and-root-controls-v10/literal-correspondence-and-root-profiles.npz',allow_pickle=False)
    curves=json.loads((BASE/'native-curves-v3/evidence.json').read_text())
    center=np.load(BASE/'native-curves-v3/literal-section-centres.npz',allow_pickle=False)
    roots=[];arrays={}
    # Protect the source-medial open boundary of each named main artery. The
    # literal profiles are held, but their inherited gap to the bifurcation is
    # not declared an established anatomical communication.
    for pid,boundary in [('ZA:Left pulmonary artery',0),('ZA:Right pulmonary artery',1)]:
        r=next(r for r in previous['native_arterial_root_profiles']if r['part_id']==pid and r['boundary']==boundary)
        pp=inp[r['profile_array']];c=pp.mean(0);key=f'root{len(roots)}';arrays[key]=pp
        roots.append(dict(part_id=pid,profile_array=key,centre_m=c.tolist(),
                          radius_m=float(np.linalg.norm(pp-c,axis=1).max()+.00025),transition_m=.012,
                          ownership='Literal source-medial main-arterial rim held for later matched recipient-port reconstruction',
                          anatomical_communication_established=False))
    # Current closed PV casts have no open LA annuli. Hold them literally.
    # Retain nearby branch-end profiles as conservative preservation anchors,
    # explicitly NOT as topological proof of their existing intersections.
    indices={pid:BVH(meshes[pid])for pid in PV}
    for row in curves['curves']:
        if row['system']!='venous':continue
        sections=center[row['key']+'_literal_sections'];cc=sections.mean(1)
        for end in [0,-1]:
            pp=sections[end];c=cc[end];radius=float(np.linalg.norm(pp-c,axis=1).max())
            nearest=min([(index.closest_point(c),pid)for pid,index in indices.items()],key=lambda q:q[0]['distance_m'])
            witness,pid=nearest
            if witness['distance_m']>radius+.001:continue
            if any(np.linalg.norm(c-np.array(r['centre_m']))<1e-6 for r in roots):continue
            key=f'root{len(roots)}';arrays[key]=pp
            roots.append(dict(part_id=row['part_id'],component=row['component'],endpoint=end,profile_array=key,
                              centre_m=c.tolist(),radius_m=radius+.00025,transition_m=.008,
                              nearby_fixed_cast_id=pid,closest_cast_face=int(witness['face_id']),
                              closest_cast_point_m=witness['point'].tolist(),closest_cast_barycentric=witness['barycentric'].tolist(),
                              center_to_cast_surface_m=witness['distance_m'],
                              ownership='Conservative literal branch-end preservation near a fixed closed cast; not an anatomical interface',
                              anatomical_communication_established=False))
    lingular=LINGULAR+[pid for pid in meshes if pid.startswith('TA98CIRCULATION:male:')]
    assert len(lingular)==9
    shared=next(r for r in previous['native_arterial_root_profiles']if r['part_id']=='ZA:Left pulmonary artery'and r['boundary']==1)
    assert shared['complete_literal_opposed_profile']and shared['opposed_edges']==33
    arrays['left_PA_lingular_shared_cut']=inp[shared['profile_array']]
    np.savez_compressed(OUT/'literal-protected-and-shared-profiles.npz',**arrays)
    evidence=dict(schema='pulmonary-whole-source-rooted-course-initialization-contract/v13',
      source_catalogue=family['source_catalogue'],source_catalogue_sha256=family['source_catalogue_sha256'],
      current_runtime_evidence_sha256=digest_file(BASE/'actual-runtime-context-v9/runtime-context.json'),
      physical_context_evidence_sha256=digest_file(BASE/'full-context-participation-v11/evidence.json'),
      physical_context_bank_hashes_sha256=digest_file(BASE/'full-context-participation-v11/bank-content-hashes-v12.json'),
      course_control_evidence_sha256=digest_file(BASE/'diverged-course-controls-v12/evidence.json'),
      protected_profile_arrays_sha256=digest_file(OUT/'literal-protected-and-shared-profiles.npz'),
      source_geometry=[dict(id=pid,geometry_sha256=m.geometry_sha,vertices=len(m.vertices),faces=len(m.faces))for pid,m in meshes.items()],
      fixed_whole_family_source_ids=FIXED,root_profile_protection=roots,
      rigid_lingular_closure=lingular,shared_lingular_main_PA_cut=shared,
      planned_field=dict(kind='New regularized Gaussian stationary velocity initialization, RK4 flow; not historical HRA replay',
        gaussian_sigma_m=.045,regularization=.35,steps=16,
        course_chart_weight=1.,common_parent_chart_weight=3.,
        exact_control_interpolation=False,proper_rotation_of_each_literal_serial_ring=True,
        exact_literal_ring_area_and_pairwise_distances_required=True,
        common_source_section_centres_share_one_flow_and_rotation=True,
        finite_lingular_cells_move_together_rigidly=True,
        left_PA_nonserial_surface='Smooth local rigid blend at its actual shared 33-edge cut and proper transport of its distal 12-point rim; actual quality gate remains pending'),
      fixed_outside_family='All other current 1524-part qualified context, including all heart, skeleton, coronary, airway and lung geometry; no catalogue writes',
      constraints=['Closed main pulmonary-vein source casts stay literal; no LA annulus or blood-space topology inferred.',
                   'Existing arterial root gaps and all outer-cast-only sources remain unaccepted until whole communicating finite-wall construction.',
                   'Shared finite lingular tissue coordinates and source face ownership stay together; no per-cell shear.',
                   'No exact anatomical registration or patency is claimed by this initialization candidate.',
                   'Positive sampled flow Jacobians do not establish discrete triangle quality or absence of contacts.'],
      geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,
      elapsed_seconds=time.monotonic()-began)
    (OUT/'contract.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),stage='contract',roots=len(roots),lingular_pieces=len(lingular),
      contract_sha256=digest_file(OUT/'contract.json'),elapsed_seconds=evidence['elapsed_seconds'])))


def smoothstep(t):
    t=np.clip(t,0.,1.);return t*t*t*(10.+t*(-15.+6.*t))


class Field:
    def __init__(self,c):
        data=json.loads((BASE/'diverged-course-controls-v12/evidence.json').read_text())
        rows=data['course_controls']+data['common_parent_controls']
        self.source=np.array([r['source_point_m']for r in rows]);self.target=np.array([r['target_point_m']for r in rows])
        self.sigma=c['planned_field']['gaussian_sigma_m'];weight=np.array([1.]*len(data['course_controls'])+[3.]*len(data['common_parent_controls']))
        delta=self.source[:,None]-self.source[None,:];kernel=np.exp(-np.sum(delta*delta,axis=2)/self.sigma**2)
        self.coeff=np.linalg.solve(kernel+np.diag(.35/weight),self.target-self.source)
        self.roots=np.array([r['centre_m']for r in c['root_profile_protection']])
        self.radius=np.array([r['radius_m']for r in c['root_profile_protection']])
        self.width=np.array([r['transition_m']for r in c['root_profile_protection']])

    def velocity(self,x):
        delta=x[:,None]-self.source[None,:];kernel=np.exp(-np.sum(delta*delta,axis=2)/self.sigma**2)
        value=kernel@self.coeff
        derivative=np.einsum('ni,nib,ia->nab',kernel,-2.*delta/self.sigma**2,self.coeff,optimize=True)
        radial=x[:,None]-self.roots[None,:];rr=np.linalg.norm(radial,axis=2)
        t=np.clip((rr-self.radius)/self.width,0.,1.);h=smoothstep(t);g=np.prod(h,axis=1)
        hp=30.*t*t*(1.-t)*(1.-t)/self.width
        loggrad=np.einsum('ni,nib->nb',hp/np.maximum(h,1e-100),radial/np.maximum(rr[:,:,None],1e-100),optimize=True)
        grad=g[:,None]*loggrad
        return g[:,None]*value,g[:,None,None]*derivative+value[:,:,None]*grad[:,None,:]

    def flow(self,points):
        outputs=[];jacobians=[]
        for start in range(0,len(points),512):
            x=np.array(points[start:start+512],float);j=np.tile(np.eye(3),(len(x),1,1));dt=1./16.
            for _ in range(16):
                def rate(xx,jj):
                    v,d=self.velocity(xx);return v,np.einsum('nab,nbc->nac',d,jj,optimize=True)
                a,aj=rate(x,j);b,bj=rate(x+.5*dt*a,j+.5*dt*aj)
                cc,cj=rate(x+.5*dt*b,j+.5*dt*bj);d,dj=rate(x+dt*cc,j+dt*cj)
                x+=dt/6.*(a+2*b+2*cc+d);j+=dt/6.*(aj+2*bj+2*cj+dj)
            outputs.append(x);jacobians.append(j)
        return np.concatenate(outputs),np.concatenate(jacobians)


def rotations(jac):
    u,s,vh=np.linalg.svd(jac);r=u@vh
    assert np.all(np.linalg.det(r)>0.),'Flow polar factor is not a proper rotation'
    return r,s


def vertex_normals(v,f):
    n=np.zeros_like(v);cross=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])
    for k in range(3):np.add.at(n,f[:,k],cross)
    return n/np.maximum(np.linalg.norm(n,axis=1,keepdims=True),1e-100)


def minimal_rotations(a,b):
    axis=np.cross(a,b);s=np.linalg.norm(axis,axis=1);c=np.sum(a*b,axis=1)
    k=np.zeros((len(a),3,3));k[:,0,1]=-axis[:,2];k[:,0,2]=axis[:,1];k[:,1,0]=axis[:,2]
    k[:,1,2]=-axis[:,0];k[:,2,0]=-axis[:,1];k[:,2,1]=axis[:,0]
    assert np.all(c>-.999999),'Antipodal source/target vertex normal needs explicit review'
    return np.eye(3)+k+(k@k)/np.maximum(1.+c,1e-12)[:,None,None]


def candidate():
    began=time.monotonic();c=json.loads((OUT/'contract.json').read_text())
    assert c['geometry_changed']is False
    family,previous,atlas,parts,reader,meshes=read_sources()
    assert family['source_catalogue_sha256']==c['source_catalogue_sha256']
    assert [(pid,m.geometry_sha)for pid,m in meshes.items()]==[(r['id'],r['geometry_sha256'])for r in c['source_geometry']]
    field=Field(c);profiles=np.load(OUT/'literal-protected-and-shared-profiles.npz',allow_pickle=False)
    anchor=profiles['left_PA_lingular_shared_cut'].mean(0)
    anchor_new,anchor_j=field.flow(anchor[None]);anchor_r,anchor_s=rotations(anchor_j);anchor_new,anchor_r=anchor_new[0],anchor_r[0]
    rigid=lambda x:(x-anchor)@anchor_r.T+anchor_new
    out=OUT/'candidate';out.mkdir();rows=[];ring_data={};all_j=[];metrics=[]
    for pi,(pid,m)in enumerate(meshes.items()):
        v,inv=np.unique(m.vertices,axis=0,return_inverse=True);faces=inv[m.faces]
        mode=None;local_rings=[];max_area=0.;max_radius=0.;max_profile=0.;sample_j=[]
        if pid in c['fixed_whole_family_source_ids']:
            vv=v.copy();mode='literal fixed proximal source envelope; no open LA/PA patency claim'
        elif pid in c['rigid_lingular_closure']:
            vv=rigid(v);mode='one shared proper rigid transform of nine-piece lingular closure'
        else:
            try:
                tubes=recover_station_rings(v,faces)
            except StationTopologyError:
                assert pid=='ZA:Left pulmonary artery',pid
                tubes=None
            if tubes is not None:
                centers=np.concatenate([v[np.array(t['rings'])].mean(1)for t in tubes])
                yy,jj=field.flow(centers);rr,ss=rotations(jj);sample_j.extend(jj)
                vv=np.full_like(v,np.nan);offset=0;mode='proper rigid transport of every literal source section'
                for ti,tube in enumerate(tubes):
                    ids=np.array(tube['rings']);source=v[ids];cs=centers[offset:offset+len(ids)]
                    local=source-cs[:,None];target=yy[offset:offset+len(ids),None]+np.einsum('nij,nkj->nki',rr[offset:offset+len(ids)],local)
                    vv[ids]=target
                    # Properly rotated literal polygons retain their section
                    # area, all pairwise distances and every local radius.
                    src_area=np.linalg.norm(np.cross(local,np.roll(local,-1,axis=1)).sum(1),axis=1)*.5
                    new_local=target-target.mean(1)[:,None]
                    dst_area=np.linalg.norm(np.cross(new_local,np.roll(new_local,-1,axis=1)).sum(1),axis=1)*.5
                    area=float(np.max(np.abs(dst_area/src_area-1.)))
                    radius=float(np.max(np.abs(np.linalg.norm(new_local,axis=2)-np.linalg.norm(local,axis=2))))
                    pair=float(np.max(np.abs(np.linalg.norm(target[:,:,None]-target[:,None,:],axis=3)-np.linalg.norm(source[:,:,None]-source[:,None,:],axis=3))))
                    max_area=max(max_area,area);max_radius=max(max_radius,radius);max_profile=max(max_profile,pair)
                    key=f'part{pi}_tube{ti}';ring_data[key+'_source_sections']=source;ring_data[key+'_candidate_sections']=target
                    ring_data[key+'_flow_jacobians']=jj[offset:offset+len(ids)];ring_data[key+'_proper_rotations']=rr[offset:offset+len(ids)]
                    local_rings.append(dict(key=key,stations=len(ids),ring_points=ids.shape[1],source_area_min_m2=float(src_area.min()),maximum_relative_area_error=area,maximum_pairwise_profile_error_m=pair))
                    offset+=len(ids)
                assert np.isfinite(vv).all()
            else:
                yy,jj=field.flow(v);rr,ss=rotations(jj);sample_j.extend(jj);vv=yy
                # Couple the existing literal left-PA/lingular 33-edge cut to
                # the one rigid finite-cell block, never cell-by-cell shear.
                radius=float(np.linalg.norm(profiles['left_PA_lingular_shared_cut']-anchor,axis=1).max()+.0005)
                weight=1.-smoothstep((np.linalg.norm(v-anchor,axis=1)-radius)/.012)
                vv=(1.-weight[:,None])*vv+weight[:,None]*rigid(v)
                # Preserve the actual distal source-rim profile by a proper
                # transform; this is not a newly constructed junction mouth.
                distal=next(r for r in previous['native_arterial_root_profiles']if r['part_id']==pid and r['boundary']==2)
                source_profiles=np.load(BASE/'takeoff-and-root-controls-v10/literal-correspondence-and-root-profiles.npz',allow_pickle=False)
                rim=source_profiles[distal['profile_array']];cc=rim.mean(0);ff,fj=field.flow(cc[None]);fr,_=rotations(fj)
                radius=float(np.linalg.norm(rim-cc,axis=1).max()+.0005)
                w=1.-smoothstep((np.linalg.norm(v-cc,axis=1)-radius)/.010)
                dest=(v-cc)@fr[0].T+ff[0];vv=(1.-w[:,None])*vv+w[:,None]*dest
                mode='nonserial left-PA velocity field plus declared shared-cut and distal-profile rigid blending; quality unaccepted'
        assert np.isfinite(vv).all()
        transformed=vv[inv];attributes={k:a.copy()for k,a in m.attributes.items()}
        if 'normals'in attributes and pid not in FIXED:
            old=vertex_normals(m.vertices,m.faces);new=vertex_normals(transformed,m.faces)
            r=minimal_rotations(old,new);raw=np.asarray(attributes['normals'],float)
            moved=np.einsum('nij,nj->ni',r,raw)
            attributes['normals']=np.rint(np.clip(moved,-32767,32767)).astype(np.int16)
        source=dict(part_id=pid,source_geometry_sha256=m.geometry_sha,source_attribute_sha256=m.attributes_sha,
                    original_source_faces_unchanged=True,mode=mode,initialization_only=True,
                    outer_cast_is_not_a_finite_wall=pid not in c['rigid_lingular_closure'],
                    contract_sha256=digest_file(OUT/'contract.json'),anatomical_acceptance=False)
        mesh=Mesh(transformed,m.faces,m.frame,attributes,source);path=out/f'part-{pi:02d}.npz';mesh.save(path)
        old_cross=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0])
        new_cross=np.cross(mesh.triangles[:,1]-mesh.triangles[:,0],mesh.triangles[:,2]-mesh.triangles[:,0])
        ratio=np.linalg.norm(new_cross,axis=1)/np.maximum(np.linalg.norm(old_cross,axis=1),1e-100)
        all_j.extend(sample_j)
        metric=dict(id=pid,mode=mode,section_tubes=len(local_rings),literal_section_count=sum(r['stations']for r in local_rings),
          maximum_relative_literal_section_area_error=max_area,maximum_radius_error_m=max_radius,
          maximum_pairwise_literal_profile_error_m=max_profile,
          displacement_max_m=float(np.linalg.norm(transformed-m.vertices,axis=1).max()),
          minimum_triangle_area_ratio=float(ratio.min()),maximum_triangle_area_ratio=float(ratio.max()),
          zero_candidate_area_faces=int((np.linalg.norm(new_cross,axis=1)<=1e-15).sum()),
          normal_sign_words_unchanged=all(np.array_equal(m.attributes[k],attributes[k])for k in ['normalSigns']if k in attributes),
          source_UV_and_color_words_unchanged=all(np.array_equal(m.attributes[k],attributes[k])for k in ['uvs','colors']if k in attributes),
          source_faces_identical=np.array_equal(m.faces,mesh.faces),source_geometry_sha256=m.geometry_sha,
          candidate_geometry_sha256=mesh.geometry_sha,mesh=str(path.resolve()),rings=local_rings)
        metrics.append(metric);rows.append(dict(id=pid,mesh=str(path.resolve()),geometry_sha256=mesh.geometry_sha))
    # Control residual is a coarse initialization measure, not proof that a
    # named artery supplies a recipient or a vein is intersegmental.
    initial=np.linalg.norm(field.source-field.target,axis=1);moved,jj=field.flow(field.source);residual=np.linalg.norm(moved-field.target,axis=1)
    all_j=np.asarray(all_j);det=np.linalg.det(all_j);singular=np.linalg.svd(all_j,compute_uv=False)
    np.savez_compressed(OUT/'literal-ring-transport-and-jacobians.npz',**ring_data,control_source=field.source,control_target=field.target,
                        control_flow=moved,control_flow_jacobians=jj,lingular_anchor_source=anchor,lingular_anchor_candidate=anchor_new,
                        lingular_proper_rotation=anchor_r,field_coefficients=field.coeff)
    assert digest_file(Path(c['source_catalogue']))==c['source_catalogue_sha256']
    result=dict(schema='pulmonary-source-rooted-coarse-course-candidate/v13',
      source_catalogue_sha256=c['source_catalogue_sha256'],contract_sha256=digest_file(OUT/'contract.json'),
      candidate_meshes=rows,part_metrics=metrics,
      literal_ring_transport_arrays_sha256=digest_file(OUT/'literal-ring-transport-and-jacobians.npz'),
      course_control_initial_rms_m=float(np.sqrt(np.mean(initial**2))),course_control_candidate_rms_m=float(np.sqrt(np.mean(residual**2))),
      course_control_initial_max_m=float(initial.max()),course_control_candidate_max_m=float(residual.max()),
      sampled_flow_jacobian_min=float(det.min()),sampled_flow_jacobian_max=float(det.max()),
      sampled_flow_minimum_singular_value=float(singular.min()),sampled_flow_samples=len(all_j),
      every_serial_ring_area_and_profile_preserved=all(r['maximum_relative_literal_section_area_error']<1e-10 and r['maximum_pairwise_literal_profile_error_m']<1e-12 for r in metrics),
      fixed_whole_sources_literal=all(r['displacement_max_m']==0. for r in metrics if r['id']in FIXED),
      finite_lingular_assembly_one_proper_rigid_transform=True,
      actual_shared_interfaces_replay_pending=True,actual_mesh_self_crossings_pending=True,
      actual_full_context_contacts_pending=True,actual_wall_lumen_patency_pending=True,render_review_pending=True,
      catalogue_changed=False,installed=False,anatomical_acceptance=False,
      limitation='Coarse source-course initialization only; inherited disconnected/open outer casts and root gaps remain, and nonserial blends still need actual discrete tissue/contact gates.',
      elapsed_seconds=time.monotonic()-began)
    (OUT/'candidate-evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['course_control_initial_rms_m','course_control_candidate_rms_m','sampled_flow_jacobian_min',
      'every_serial_ring_area_and_profile_preserved','fixed_whole_sources_literal','elapsed_seconds']}|dict(output=str(OUT),
      candidate_evidence_sha256=digest_file(OUT/'candidate-evidence.json'))))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['contract','candidate'],required=True);args=parser.parse_args()
    if args.stage=='contract':contract()
    else:candidate()
