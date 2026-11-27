"""Finite urinary outlet from explicit bladder lips and a native course.

The native urethra is a course/profile guide, not evidence of wall thickness.
The illustrative external coat is declared in metres. Only paired wall rims
are connected; no disk spans the fluid or replaces the prostate's real bore.
"""
from pathlib import Path
import json

import numpy as np

from .algorithms.annular_pipe_envelopes import boundary_loops
from .algorithms.tube_continuity import annular_rim
from .algorithms.tube_rim_bridge import bridge_boundaries, end_annulus
from .algorithms.tube_sections import arclength
from .algorithms.tube_station_topology import recover_station_rings
from .algorithms.finite_tissue_boolean import make_solid, provenance, barycentric, classifier_caps_outside
from .algorithms.surface_normals import grouped_normal_words
from .audit import audit
from .bladder_carriers import recover
from .common import atomic_json, atomic_npz, digest_file
from .model import Mesh


def build(joined_dir, outlet_sources, output, *, wall_thickness_m=.00045, collar_length_m=.008):
    joined_dir,outlet_sources,output=map(lambda p:Path(p).resolve(),(joined_dir,outlet_sources,output))
    if output.exists():
        raise FileExistsError(output)
    inventory=json.loads((outlet_sources/'inventory-and-assembly.json').read_text())
    if not inventory['terms'] or len(inventory['terms'])!=49:
        raise ValueError('Complete source-bound outlet inventory is required')
    if not .0001 <= wall_thickness_m <= .002 or not .003 <= collar_length_m <= .015:
        raise ValueError('Explicit bounded finite-coat and local collar construction inputs required')
    report=json.loads((joined_dir/'junctions.json').read_text())
    if report.get('passed') is not True:
        raise ValueError('Begin with the verified continuous bladder/ureter wall')
    wall=Mesh.load(joined_dir/'verified-carriers/wall.npz')
    with np.load(joined_dir/'verified-carriers/carriers.npz',allow_pickle=False)as z:
        if not np.array_equal(wall.faces,z['wall_faces']):
            raise ValueError('Joined wall/carrier bindings changed')
        inner_ids=z['inner_source_face_ids'].copy();inner=z['inner_side_faces'].copy();neck=z['urethral_neck_ring'].copy()
    body_annulus,body_outer=end_annulus(wall.faces,inner,neck)
    native=Mesh.load(outlet_sources/'urethral-course.npz')
    expected=inventory['sources']['ZA:Urethra']
    if native.geometry_sha!=expected['geometry_sha256']:
        raise ValueError('Native urethral course changed')
    v,inverse=np.unique(native.vertices,axis=0,return_inverse=True);f=inverse[native.faces]
    rows=recover_station_rings(v,f)
    if len(rows)!=1:
        raise ValueError('One complete native urethral course is required')
    rings=np.array(rows[0]['rings']);axis=v[rings].mean(1)
    if axis[0,1]<axis[-1,1]:
        rings=rings[::-1];axis=axis[::-1]
    arc=arclength(axis);k=int(np.searchsorted(arc,collar_length_m))
    if not 1<=k<len(rings)-2:
        raise ValueError('The collar requires a source-connected downstream station')
    retained=np.concatenate(rings[k:]);source_faces=np.all(np.isin(f,retained),axis=1)
    used,inv=np.unique(f[source_faces],return_inverse=True);profile=v[used];pf=inv.reshape(-1,3)
    station=np.full(len(v),-1,int)
    for i,ring in enumerate(rings):
        station[ring]=i
    radial=profile-axis[station[used]];length=np.linalg.norm(radial,axis=1)
    if np.any(length<=0):
        raise ValueError('A source profile collapses onto its course')
    # The complete literal native profile is the luminal construction guide.
    # The external wall is offset in each source section, never narrowed.
    outer=profile+wall_thickness_m*radial/length[:,None]
    q=profile[pf];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);r=radial[pf].mean(1)
    if np.median(np.einsum('ij,ij->i',n,r))<0:
        pf=pf[:,[0,2,1]]
    count=len(profile);vertices=np.r_[outer,profile]
    rims=boundary_loops(pf)
    if len(rims)!=2:
        raise ValueError('The retained urethral course does not have two actual mouths')
    proximal=max(rims,key=lambda r:profile[r].mean(0)[1])
    distal=min(rims,key=lambda r:profile[r].mean(0)[1])
    annuli=[annular_rim(r,r+count)for r in rims]
    uf=np.r_[pf,pf[:,[0,2,1]]+count,*annuli]
    urethra=Mesh(vertices,uf,wall.frame)
    ucheck=audit(urethra,require_closed=True,expected_components=1,expected_genus=[1],max_pairs=1000000)
    output.mkdir(parents=True);urethra.save(output/'urethral-wall.npz')
    if not ucheck['passed']:
        atomic_json(output/'urethral-failure.json',ucheck)
        raise ValueError('The finite urethral wall fails actual geometry checks')
    remove=np.concatenate([body_annulus]);keep=np.setdiff1d(np.arange(len(wall.faces)),remove)
    inner_set=set(map(int,inner_ids));inside=[j for j,i in enumerate(keep)if int(i)in inner_set]
    # Keep the physical annulus at the external meatus; remove only the
    # proximal annulus that becomes an internal shared tissue continuation.
    uremove=annuli[0]if set(rims[0])==set(proximal)else annuli[1]
    keys={tuple(sorted(t))for t in uremove}
    ukept=np.array([i for i,t in enumerate(uf)if tuple(sorted(t))not in keys],int)
    offset=len(wall.vertices);all_vertices=np.r_[wall.vertices,vertices]
    pieces=[wall.faces[keep],uf[ukept]+offset]
    inside.extend(len(pieces[0])+j for j,i in enumerate(ukept)if len(pf)<=i<2*len(pf))
    # Current boundaries are read from the actual retained walls; do not
    # infer matching cycles by radius or nearest anatomical neighbour.
    u_outer=min(boundary_loops(uf[ukept]),key=lambda r:abs(vertices[r].mean(0)[1]-profile[proximal].mean(0)[1]) if np.all(r<count)else 1.)
    u_inner=next(r for r in boundary_loops(uf[ukept])if np.all(r>=count))
    outer_bridge,outer_evidence=bridge_boundaries(all_vertices,body_outer,u_outer+offset)
    inner_bridge,inner_evidence=bridge_boundaries(all_vertices,neck,u_inner+offset)
    pieces.append(outer_bridge);inside.extend(range(sum(len(p)for p in pieces),sum(len(p)for p in pieces)+len(inner_bridge)));pieces.append(inner_bridge)
    whole=Mesh(all_vertices,np.concatenate(pieces),wall.frame);whole.save(output/'urinary-wall.npz')
    check=audit(whole,require_closed=True,expected_components=1,expected_genus=[2],max_pairs=1000000)
    atomic_npz(output/'source-lineage.npz',inner_source_face_ids=np.array(inside,int),
        native_course_axis=axis,source_station_arclength_m=arc,source_station_ids=station[used],
        native_luminal_guide=profile,bladder_neck_ring=neck,outer_bridge=outer_bridge,inner_bridge=inner_bridge)
    result=dict(schema='finite-bladder-urethral-outlet/v1',passed=check['passed'],whole_check=check,urethral_check=ucheck,
        source_inventory_sha256=digest_file(outlet_sources/'inventory-and-assembly.json'),
        source_bladder_ureter_join_sha256=digest_file(joined_dir/'junctions.json'),
        wall_thickness_m=wall_thickness_m,clinical_thickness_measured=False,native_section_guide_preserved=True,
        source_first_retained_station=k,source_collar_arclength_m=float(arc[k]),
        outer_bridge=outer_evidence,inner_bridge=inner_evidence,physiological_disk_emitted=False,
        prostate_modified=False,skeleton_modified=False,installed=False,anatomical_acceptance=False)
    atomic_json(output/'outlet.json',result)
    if check['passed']:
        result['carriers']=recover(whole,output/'verified-carriers',inner_face_ids=np.array(inside,int),
            mouth_labels=['external-urethral-orifice','right-renal-root','left-renal-root'],cap_forbid_source_diagonals=True)
        atomic_json(output/'outlet-with-carriers.json',result)
    return result


def prostate_bore(outlet_dir, outlet_sources, output):
    """Remove the actual urinary envelope from the retained prostate volume."""
    outlet_dir,outlet_sources,output=map(lambda p:Path(p).resolve(),(outlet_dir,outlet_sources,output))
    if output.exists():
        raise FileExistsError(output)
    passed=json.loads((outlet_dir/'outlet-with-carriers.json').read_text())
    if passed.get('passed') is not True:
        raise ValueError('Use a passed, open bladder/urethral outlet')
    source=Mesh.load(outlet_sources/'prostate.npz');ps,pid=make_solid(source.vertices,source.faces)
    expected=json.loads((outlet_sources/'inventory-and-assembly.json').read_text())['sources']['ZA:Prostate']
    if source.geometry_sha!=expected['geometry_sha256'] or source.attributes_sha!=expected['attributes_sha256']:
        raise ValueError('The protected prostate source changed')
    carrier=outlet_dir/'verified-carriers/carriers.npz'
    with np.load(carrier,allow_pickle=False)as z:
        ov,of=z['vertices'].copy(),z['outer_faces'].copy();ff=z['fluid_faces'].copy();disks=z['computational_disks'].copy()
    cap_keys={tuple(sorted(t))for t in disks}
    cap_ids=np.array([i for i,t in enumerate(of)if tuple(sorted(t))in cap_keys],int)
    if len(cap_ids)!=len(disks):
        raise ValueError('Computational urinary end disks lost their identity')
    cap_check=classifier_caps_outside(ov,of,cap_ids,source.vertices,source.faces)
    occupied,oid=make_solid(ov,of);fluid,_=make_solid(ov,ff)
    result=ps-occupied;m=result.to_mesh64();v=np.asarray(m.vert_properties)[:,:3].copy();f=np.asarray(m.tri_verts,int).copy()
    owners,parents=provenance(m);old=owners==pid;new=owners==oid
    if not np.all(old|new):
        raise ValueError('The prostate bore lost source-face ancestry')
    words=np.empty((len(f),3,3),dtype='<i2');uv=np.empty((len(f),3,2),dtype='<f4')
    original_words=source.attributes['normals'].astype(np.int32)
    if 'normalSigns'in source.attributes:
        original_words*=source.attributes['normalSigns'][:,None]
    if 'uvs'not in source.attributes:
        raise ValueError('Retain the actual prostate source texture coordinates')
    maximum_plane_error=0.
    for i in np.flatnonzero(old):
        parent=source.faces[parents[i]];triangle=source.vertices[parent];b=barycentric(v[f[i]],triangle)
        reconstruction=b@triangle
        maximum_plane_error=max(maximum_plane_error,float(np.linalg.norm(reconstruction-v[f[i]],axis=1).max()))
        uv[i]=(b@source.attributes['uvs'][parent]).astype('<f4')
        words[i]=np.rint(b@original_words[parent]).astype('<i2')
    bore_faces=f[new]
    if not len(bore_faces)or maximum_plane_error>2e-12:
        raise ValueError('The bore is absent or the retained exterior left its original source planes')
    smooth,normal_evidence=grouped_normal_words(v,{'bore':bore_faces},{'bore':['gland-urinary-contact']*len(bore_faces)})
    words[new]=smooth['bore']
    # New internal contact charts borrow only a nearest-source UV color map;
    # this is not a new histological texture or a clinical tissue measurement.
    import point_cloud_utils as pcu
    _,nearest,bary=pcu.closest_points_on_mesh(v[bore_faces].reshape(-1,3),source.vertices,np.ascontiguousarray(source.faces,dtype=np.int32))
    uv[new]=np.einsum('ij,ijk->ik',bary,source.attributes['uvs'][source.faces[nearest]]).reshape(-1,3,2)
    xyz=v[f];fn=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
    if np.any(np.einsum('tki,ti->tk',words.astype(float),fn)<=0):
        raise ValueError('The retained/new prostate normal field is not outward')
    triangles=np.arange(len(f)*3).reshape(-1,3)
    mesh=Mesh(xyz.reshape(-1,3),triangles,source.frame,
        dict(normals=words.reshape(-1,3),normalSigns=np.ones(len(f)*3,dtype='i1'),uvs=uv.reshape(-1,2)))
    check=audit(mesh,require_closed=True,expected_components=1,expected_genus=[1],max_pairs=1000000)
    escape=float((result^fluid).volume())
    bound=2048*np.finfo(float).eps*float(np.ptp(source.vertices,axis=0).max())**3*len(f)
    output.mkdir(parents=True);mesh.save(output/'prostate.npz')
    atomic_npz(output/'source-provenance.npz',vertices=v,faces=f,prostate_source_faces=np.where(old,parents,-1),
        urinary_outer_source_faces=np.where(new,parents,-1))
    report=dict(schema='source-bound-prostatic-urinary-bore/v1',passed=bool(check['passed']and abs(escape)<=bound),
        prostate_check=check,urinary_fluid_overlap_m3=escape,numerical_volume_bound_m3=bound,
        prostate_source_sha256=digest_file(outlet_sources/'prostate.npz'),urinary_carriers_sha256=digest_file(carrier),
        source_exterior_plane_maximum_error_m=maximum_plane_error,
        source_exterior_uvs_replayed=True,source_exterior_normals_interpolated_only_on_split_faces=True,
        new_internal_uvs_nearest_source_projection=True,new_internal_normal_evidence=normal_evidence,
        computational_end_disks_outside_gland=cap_check,source_volume_m3=float(ps.volume()),
        retained_prostate_volume_m3=float(result.volume()),removed_occupied_urinary_envelope_m3=float((ps^occupied).volume()),
        prostate_exterior_transformed=False,skeleton_modified=False,installed=False,anatomical_acceptance=False,
        remaining=['Reciprocal finite-coat/prostate interface subdivision and ordinary atlas installation.',
                   'Actual surface/section and full pelvic packing review; no complete prostatic anatomy claim.'])
    atomic_json(output/'prostate-bore.json',report)
    return report


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--joined-dir',required=True);p.add_argument('--outlet-sources',required=True);p.add_argument('--output',required=True)
    r=build(**vars(p.parse_args()));print(json.dumps(dict(passed=r['passed'],failures=r['whole_check']['failures']),indent=2))
