"""Source-bound adult lingual-septum proposal; never changes an atlas.

Inspect the literal current midsagittal muscle geometry before constructing one
finite sheet. The adult histological reference constrains relationships and a
0.1--0.2 mm thickness, but cannot supply this donor's exact septal boundary.
All occupied-volume conflicts remain explicit; a closed mesh is not an
anatomically completed connective-tissue cell.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'server'))

import numpy as np
import igl
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import diagnostics
from tongue_whole_input import WORK, SOURCE_AUDIT, fixed_bone_sources

REFERENCE = 'https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/ar.25198'
SELECTED = WORK / 'selected-coupling-audit-v2/coupled-inputs.selected.json'
FRAME = 'model:local-ta98:assembled-metres'


def require(value, message):
    if not value:
        raise ValueError(message)


def bound(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=digest_file(path))


def load_bound(row):
    path = Path(row['path'])
    require(path.is_relative_to(WORK) and digest_file(path) == row['sha256'],
            'A selected geometry source changed')
    return Mesh.load(path)


def inputs():
    manifest = json.loads(SELECTED.read_text())
    require(manifest['packingReady'] and not manifest['reviewPending'],
            'The selected coupled inputs are not ready')
    b = manifest['bindings']
    meshes = {'whole': b['whole']}
    intrinsic = json.loads(Path(b['intrinsics']['path']).read_text())
    for row in intrinsic['meshes']:
        if row['key'] != 'superior-longitudinal':
            meshes[row['key']] = dict(path=row['mesh'], sha256=row['sha256'])
    packed = WORK / 'shared-coverings-exposed-apo-sl-v3/packed'
    for key, filename in [('superior-longitudinal', 'intrinsic_superior-longitudinal'),
                          ('mucosa', 'mucosa'), ('aponeurosis', 'aponeurosis')]:
        meshes[key] = bound(packed / f'TA98TONGUE_male_{filename}.npz')
    extrinsic = json.loads(Path(b['extrinsics']['path']).read_text())
    for row in extrinsic['parts']:
        meshes[row['name']] = dict(path=row['mesh'], sha256=row['sha256'])
    for key, row in b.items():
        require(digest_file(row['path']) == row['sha256'], 'Changed selected binding: '+key)
    return meshes, {key: load_bound(row) for key, row in meshes.items()}


def ray_y(mesh, x, z):
    """All actual facet intersections of a Y-directed line, with witnesses."""
    tri = mesh.triangles
    edge = tri[:, 1:] - tri[:, :1]
    a, b = edge[:, 0, 0], edge[:, 1, 0]
    c, d = edge[:, 0, 2], edge[:, 1, 2]
    det = a*d-b*c
    good = np.abs(det) > 1e-18
    u = np.zeros(len(tri)); v = u.copy()
    qx, qz = x-tri[:, 0, 0], z-tri[:, 0, 2]
    u[good] = (d[good]*qx[good]-b[good]*qz[good])/det[good]
    v[good] = (-c[good]*qx[good]+a[good]*qz[good])/det[good]
    selected = np.flatnonzero(good & (u >= -1e-11) & (v >= -1e-11) & (u+v <= 1+1e-11))
    hits = []
    for f in selected:
        weights = np.array([1-u[f]-v[f], u[f], v[f]])
        point = weights @ tri[f]
        if any(abs(row['y_m']-point[1]) < 1e-10 for row in hits):
            continue
        hits.append(dict(y_m=float(point[1]), face=int(f), barycentric=weights.tolist()))
    return sorted(hits, key=lambda row: row['y_m'])


def plane_segments(mesh, x):
    segments = []
    for tri in mesh.triangles:
        hits = []
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            a, b = tri[i], tri[j]
            if (a[0]-x)*(b[0]-x) < 0:
                hits.append(a+(x-a[0])/(b[0]-a[0])*(b-a))
            elif a[0] == x:
                hits.append(a)
        points = np.unique(np.asarray(hits), axis=0) if hits else []
        if len(points) == 2:
            segments.append(points)
    return np.asarray(segments)


def inspect(output):
    began = time.monotonic()
    require(output.is_relative_to(WORK) and not output.exists(), 'Use a new local stage')
    bindings, meshes = inputs()
    output.mkdir(parents=True)
    positions = np.linspace(.029, .079, 51)
    stations = []
    for z in positions:
        row = dict(z_m=float(z), midline={})
        for key in ['whole', 'superior-longitudinal', 'transverse-vertical-shared']:
            row['midline'][key] = ray_y(meshes[key], 0., z)
        # Medial GG surfaces may stop at X=0. Inspect both actual medial strips
        # rather than infer the GG roof from a bounding box or empty midline ray.
        gg = []
        for key, sign in [('genioglossus-left', 1), ('genioglossus-right', -1)]:
            for x in sign*np.array([.0004,.001,.002,.003]):
                gg.extend(dict(muscle=key, x_m=float(x), **hit)
                          for hit in ray_y(meshes[key], x, z))
        row['genioglossus_medial_hits'] = gg
        tv = row['midline']['transverse-vertical-shared']
        sl = row['midline']['superior-longitudinal']
        if len(tv) == 2 and len(sl) == 2:
            # Reserve the lowermost T/V stratum for actual crossing fibers.
            # One millimetre is an explicit macroscopic display estimate, not a
            # measurement from the reference or a patient-specific landmark.
            floor = max(tv[0]['y_m']+.001,
                        max([h['y_m'] for h in gg], default=-np.inf)+.0003)
            roof = min(tv[-1]['y_m']-.0002, sl[0]['y_m']-.0002)
            row.update(proposed_floor_y_m=float(floor), proposed_roof_y_m=float(roof),
                       available_height_m=float(roof-floor),
                       construction_station_supported=bool(roof-floor > .001))
        stations.append(row)
    atomic_json(output/'extent-constraints.json', dict(
        schema='tongue-septum-extent-constraints/v1',
        source_selected=bound(SELECTED), inputs=bindings,
        authority_inventory=bound(ROOT/'reports/ta98-tongue-inventory.json'),
        whole_organ='Male tongue, source-bound current coupled assembly',
        planned_part=dict(term='A05.1.04.024', parent='A05.1.04.001',
            representation='finite-connective-tissue-sheet', name='lingual septum'),
        primary_reference=REFERENCE,
        reference_scope='Adult histological sections; not donor-specific 3D segmentation',
        reference_constraints=dict(thickness_m=[.0001,.0002], wavy=True,
            below_superior_longitudinal=True, above_lowest_transverse_crossing_stratum=True,
            no_sheet_between_bilateral_genioglossus=True),
        proposed_estimates=dict(thickness_m=.00015, transverse_crossing_stratum_reserve_m=.001,
            intrinsic_boundary_reserve_m=.0002, medial_GG_roof_reserve_m=.0003),
        adjoining_interfaces=[
            dict(second='superior-longitudinal', relation='disjoint', method='actual source ray inferior boundary'),
            dict(second='genioglossus-left/right', relation='disjoint', method='above actual medial roof; no inferior midline wall'),
            dict(second='transverse-vertical-shared', relation='shared-tissue-interface-required',
                 state='pending reciprocal volume ownership; no overlap waiver'),
            dict(second='mucosa/aponeurosis/IL/HG/fixed bones', relation='disjoint')],
        stations=stations, skeleton_modified=False, source_files_modified=False,
        installed=False, anatomical_acceptance=False,
        elapsed_seconds=time.monotonic()-began))
    # Actual facet-plane cuts. This is a source inspection, not an imagined
    # diagram and not a capped solid or a resampled proprietary reference.
    render_sections(output, meshes, stations)
    print(json.dumps(dict(stage='inspected', supported_stations=sum(
        r.get('construction_station_supported',False) for r in stations),
        constraints=str(output/'extent-constraints.json'), image=str(output/'actual-source-sections.png'),
        elapsed_seconds=time.monotonic()-began)), flush=True)


def render_sections(output, meshes, stations):
    """Tiny dependency-free RGB line plot of actual triangle-plane segments."""
    import struct
    import zlib
    width, height, scale = 1740, 640, 7600.
    pixels = np.full((height,width,3), 250, np.uint8)
    def line(first, last, color, weight=1):
        count = max(2,int(np.max(np.abs(last-first)))+1)
        points = np.rint(np.linspace(first,last,count)).astype(int)
        for dx in range(-weight,weight+1):
            for dy in range(-weight,weight+1):
                x,y=points[:,0]+dx,points[:,1]+dy
                good=(x>=0)&(x<width)&(y>=0)&(y<height)
                pixels[y[good],x[good]]=color
    colors = {'whole':[45,45,45], 'superior-longitudinal':[179,71,84],
              'transverse-vertical-shared':[210,164,45],
              'genioglossus-left':[62,139,186], 'genioglossus-right':[62,139,186]}
    for index,x in enumerate([-.002,0.,.002]):
        def point(q):
            return np.array([index*580+30+(q[2]-.017)*scale,
                             height-35-(q[1]-1.587)*scale])
        for z in np.arange(.020,.090,.010):
            line(point([0,1.590,z]),point([0,1.660,z]),[230,230,230],0)
        for y in np.arange(1.590,1.661,.010):
            line(point([0,y,.020]),point([0,y,.088]),[230,230,230],0)
        for key,color in colors.items():
            for a,b in plane_segments(meshes[key],x):line(point(a),point(b),color)
        for key,color in [('proposed_floor_y_m',[88,141,53]),('proposed_roof_y_m',[145,60,147])]:
            supported=[r for r in stations if r.get('construction_station_supported')]
            for a,b in zip(supported,supported[1:]):
                line(point([0,a[key],a['z_m']]),point([0,b[key],b['z_m']]),color)
    raw=b''.join(b'\0'+pixels[y].tobytes()for y in range(height))
    def chunk(kind,data):
        return struct.pack('>I',len(data))+kind+data+struct.pack('>I',zlib.crc32(kind+data)&0xffffffff)
    image=b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',width,height,8,2,0,0,0))
    image+=chunk(b'IDAT',zlib.compress(raw))+chunk(b'IEND',b'')
    (output/'actual-source-sections.png').write_bytes(image)
    atomic_json(output/'source-section-legend.json',dict(section_x_m=[-.002,0,.002],
        source_colors=colors,lower_curve_rgb=[88,141,53],upper_curve_rgb=[145,60,147],
        proposed_curve_not_anatomical_boundary=True,grid_spacing_m=.010,
        plot_scale_pixels_per_m=scale,source_image=bound(output/'actual-source-sections.png')))


def construct(output, extent):
    """One bounded, visibly qualified finite proposal in the supported extent."""
    began=time.monotonic()
    require(output.is_relative_to(WORK) and not output.exists(), 'Use a new stage')
    plan=json.loads(extent.read_text())
    require(plan['schema']=='tongue-septum-extent-constraints/v1' and
            digest_file(SELECTED)==plan['source_selected']['sha256'], 'Extent inputs changed')
    bindings,meshes=inputs()
    require(bindings==plan['inputs'], 'Actual source meshes changed after inspection')
    supported=[r for r in plan['stations']if r.get('construction_station_supported')]
    require(all(b['z_m']-a['z_m']<.00101 for a,b in zip(supported,supported[1:])),
            'Do not bridge an unsupported axial interval')
    require(len(supported)>=10, 'No defensible connected extent')
    # The supported posterior/body extent is a prototype; it is explicitly NOT
    # the entire MLS. Smooth limiting ends terminate inside the current T/V
    # territory until its anterior GG architecture can be resolved.
    z=np.array([r['z_m']for r in supported])
    floor=np.array([r['proposed_floor_y_m']for r in supported])
    roof=np.array([r['proposed_roof_y_m']for r in supported])
    center=(floor+roof)/2
    half=(roof-floor)/2-.00025
    require(np.all(half>.00025), 'The reference-constrained sheet collapsed')
    t=(z-z[0])/(z[-1]-z[0])
    # Reduced Y width at the two finite ends gives a rounded outline; the
    # positive X thickness is retained until the final convex apex caps.
    axial=np.sin(np.pi*t)**.35
    theta=np.arange(12)*2*np.pi/12
    rings=[]
    for i in range(1,len(z)-1):
        yy=center[i]+half[i]*axial[i]*np.sin(theta)
        wave=.00015*np.sin(2*np.pi*(z[i]-z[0])/.010+1.2*np.sin(theta))
        xx=wave+.000075*np.cos(theta)
        rings.append(np.c_[xx,yy,np.full(12,z[i])])
    v=np.r_[np.array([[0.,center[0],z[0]]]),np.concatenate(rings),
            np.array([[0.,center[-1],z[-1]]])]
    faces=[];count=len(rings)
    for j in range(12):faces.append([0,1+(j+1)%12,1+j])
    for i in range(count-1):
        for j in range(12):
            a=1+i*12+j;b=1+i*12+(j+1)%12;c=a+12;d=b+12
            faces.extend([[a,b,c],[b,d,c]])
    last=len(v)-1;start=1+(count-1)*12
    for j in range(12):faces.append([last,start+j,start+(j+1)%12])
    f=np.asarray(faces,int)
    volume=float(np.sum(np.einsum('ij,ij->i',v[f[:,0]],np.cross(v[f[:,1]],v[f[:,2]])))/6)
    if volume<0:f=f[:,[0,2,1]]
    origin=np.array([0.,1.625,0.])
    stored=(v-origin).astype('<f4').astype(float)+origin
    mesh=Mesh(stored,f,FRAME,{'normals':igl.per_vertex_normals(stored,f)},dict(
        kind='original-reference-constrained-finite-MLS-prototype',term='A05.1.04.024',
        source_selected_sha256=digest_file(SELECTED),extent_constraints_sha256=digest_file(extent),
        nominal_central_thickness_m=.00015,primary_reference=REFERENCE,
        source_geometry_copied=False,source_geometry_modified=False,anatomical_acceptance=False,
        anterior_extent_missing=True,reciprocal_transverse_attachment_pending=True))
    output.mkdir(parents=True);path=output/'septum.npz';mesh.save(path)
    selfcheck=diagnostics(mesh);topo=topology(mesh)
    sample=np.r_[mesh.vertices,mesh.triangles.mean(1)]
    winding={key:float(np.abs(igl.winding_number(meshes[key].vertices,meshes[key].faces,sample)).min())
             for key in ['whole','transverse-vertical-shared']}
    pairs=[]
    disjoint=['superior-longitudinal','mucosa','aponeurosis','inferior-longitudinal-left',
              'inferior-longitudinal-right','genioglossus-left','genioglossus-right',
              'hyoglossus-left','hyoglossus-right']
    for key in disjoint:
        _,summary=pair_contacts(mesh,meshes[key],0.,max_pairs=100000,max_tests=1000000)
        contained=np.abs(igl.winding_number(meshes[key].vertices,meshes[key].faces,sample))
        pairs.append(dict(second=key,summary=summary,
            maximum_septum_sample_winding=float(contained.max()),
            sampled_occupied_volume_exclusion=bool(np.all(contained<.01))))
    bones=[]
    for row in fixed_bone_sources():
        bone=Mesh.load(row['meshSavePath'])
        _,summary=pair_contacts(mesh,bone,0.,max_pairs=100000,max_tests=1000000)
        bones.append(dict(part_id=row['partId'],geometry_sha256=row['geometrySha256'],summary=summary))
    clear=lambda row:row['summary']['complete_within_pair_budget']and not(
        row['summary']['crossing_face_pairs']or row['summary']['coplanar_overlap_face_pairs'])
    passed=bool(topo['closed_oriented_manifold'] and len(topo['components'])==1 and
        selfcheck['self_intersection_pairs']==0 and min(winding.values())>.99 and
        all(clear(row)for row in pairs+bones)and
        all(row['sampled_occupied_volume_exclusion']for row in pairs))
    atomic_json(output/'candidate.json',dict(schema='tongue-septum-original-prototype/v1',
        passed_geometry=passed,anatomical_acceptance=False,installed=False,installation_ready=False,
        term='A05.1.04.024',parent='A05.1.04.001',representation='finite-connective-tissue-sheet',
        mesh=bound(path),geometry_sha256=mesh.geometry_sha,inputs=bindings,
        extent_constraints=bound(extent),source_selection=bound(SELECTED),
        inspected_source_sections=bound(extent.parent/'actual-source-sections.png'),
        source_sections_inspected=True,primary_reference=REFERENCE,
        modeled_axial_extent_m=[float(z[0]),float(z[-1])],
        nominal_central_thickness_m=.00015,estimated_lateral_waviness_m=.00015,
        dimensions_are_modeling_estimates=True,finite_outline_is_provisional=True,
        topology=topo,cpp_diagnostics=selfcheck,sample_containment_winding=winding,
        disjoint_pairs=pairs,fixed_bone_pairs=bones,
        ownership=dict(septum='finite tissue prototype',
            transverse_vertical='Currently contains the prototype volume; reciprocal subtraction/attachment chart required',
            occupied_volume_disjointness=False,overlap_permitted=False),
        remaining=['Anterior septal extent cannot be established from current coarse GG/SL carriers',
            'Reciprocal septum/T-V tissue ownership and epimysial attachment not constructed',
            'Actual Blender surface/section anatomical review still required'],
        resource=dict(vertices=len(v),triangles=len(f),estimated_raw_packed_bytes=len(v)*18+len(f)*12,
            image_textures=0),source_files_modified=False,skeleton_modified=False,
        storage_origin_m=origin.tolist(),float32_error_m=float(np.linalg.norm(stored-v,axis=1).max()),
        elapsed_seconds=time.monotonic()-began))
    atomic_json(output/'review-spec.json',dict(frame=FRAME,resolution=640,wall_seconds=50,
        meshes=[dict(mesh=str(path),name='Septum prototype — anterior extent and tissue ownership pending',
                     color=[.78,.72,.65,1.]),
                dict(mesh=bindings['superior-longitudinal']['path'],name='Actual SL context',color=[.7,.28,.3,.15]),
                dict(mesh=bindings['genioglossus-left']['path'],name='Actual left GG context',color=[.36,.48,.66,.15])],
        views=[dict(name='right',direction=[1,0,0]),dict(name='front',direction=[0,0,1]),
               dict(name='coronal-section',direction=[0,0,1],section=dict(normal=[0,0,1],
                    point=[0,1.635,.045],clear_outer=True))],
        qualification='Uninstalled finite MLS prototype; incomplete anterior extent and reciprocal tissue ownership'))
    print(json.dumps(dict(passed_geometry=passed,anatomical_acceptance=False,
        triangles=len(f),output=str(output),elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--extent',type=Path)
    args = parser.parse_args()
    if args.extent:construct(args.output.resolve(),args.extent.resolve())
    else:inspect(args.output.resolve())
