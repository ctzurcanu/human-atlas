"""One source-bound subdorsal-sheet / shared T-V territory proposal.

Offsets are explicit illustrative construction fields, not measured histology.
Self-crossing input offsets stop construction before any Boolean repairs.
"""
from pathlib import Path
import argparse
import json
import sys
import time
from collections import defaultdict
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file,atomic_npz
from geometry_tools.audit import topology
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.contact import pair_contacts
from geometry_tools.tissue import layer_arrangement


def preflight(out,erosion_check=None):
    began=time.monotonic();out.mkdir(parents=True,exist_ok=True)
    work=out.parent
    source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    controls=work/'surface-controls-v1/surface-controls.json'
    data=json.loads(controls.read_text())
    assert digest_file(source)==data['sourceCatalogueSha256']
    native=Mesh.load(work/'source-v1/meshes/native-00.npz')
    normals=native.attributes['normals'].astype(float)
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    origin=np.array([0.,1.625,.055])
    def emit(vertices,faces,name,provenance):
        # All subsequent checks use exactly this bounded Float32 payload.
        emitted=(vertices-origin).astype(np.float32).astype(float)+origin
        mesh=Mesh(emitted,faces,native.frame,{'normals':igl.per_vertex_normals(emitted,faces)},provenance)
        path=out/(name+'.npz');mesh.save(path)
        q=np.concatenate([mesh.vertices,mesh.triangles.mean(1),
            (mesh.triangles[:,0]+mesh.triangles[:,1])/2,(mesh.triangles[:,1]+mesh.triangles[:,2])/2,
            (mesh.triangles[:,2]+mesh.triangles[:,0])/2])
        winding=np.abs(igl.winding_number(native.vertices,native.faces,q))
        _,outer=pair_contacts(native,mesh,0.,max_pairs=100000,max_tests=1000000)
        cpp=diagnostics(mesh);topo=topology(mesh)
        checks=dict(closed_oriented=topo['closed_oriented_manifold'],
            self_crossings_zero=cpp['self_intersection_pairs']==0,contained=bool(np.all(winding>.99)),
            outer_crossings_zero=outer['crossing_face_pairs']==0 and outer['coplanar_overlap_face_pairs']==0,
            contact_budget_complete=outer['complete_within_pair_budget'])
        return mesh,dict(path=str(path),sha256=digest_file(path),geometry_sha256=mesh.geometry_sha,
            topology=topo,cpp_diagnostics=cpp,minimum_absolute_winding=float(winding.min()),
            outer_contacts=outer,checks=checks,passed=all(checks.values()))
    if erosion_check is None:
        inner,inner_check=emit(native.vertices-.0009*normals,native.faces,'native-inner-envelope',dict(
            kind='actual-source-normal-inward-construction',source_geometry_sha256=native.geometry_sha,
            source_catalogue_sha256=digest_file(source),normal_source='literal stored Int16 source normals, normalized',
            inward_offset_m=.0009,approximate_territory=True,anatomical_acceptance=False))
        atomic_json(out/'inner-offset-check.json',inner_check)
    else:
        inner_check=json.loads(erosion_check.read_text())
        assert inner_check['passed'] and inner_check['source_catalogue_sha256']==digest_file(source)
        assert digest_file(inner_check['mesh'])==inner_check['mesh_sha256']
        inner=Mesh.load(inner_check['mesh'])
        assert inner.geometry_sha==inner_check['geometry_sha256']
    if not inner_check['passed']:
        failure=dict(schema='tongue-coherent-territory-preflight/v1',passed=False,
            state='stopped-native-inward-offset-invalid',source_catalogue_sha256=digest_file(source),
            inner=inner_check,source_geometry_modified=False,installed=False,anatomical_acceptance=False,
            elapsed_seconds=time.monotonic()-began)
        atomic_json(out/'preflight.json',failure)
        atomic_json(out/'review-spec.json',dict(frame=native.frame,resolution=640,wall_seconds=60,
            meshes=[dict(mesh=str(out/'native-inner-envelope.npz'),name='Rejected inward offset',color=[.75,.35,.32,1.])],
            views=[dict(name='top',direction=[0,1,0]),dict(name='right',direction=[1,0,0])]))
        print(json.dumps(dict(state=failure['state'],checks=inner_check['checks'],
            self_intersections=inner_check['cpp_diagnostics']['self_intersection_pairs'],elapsed_seconds=failure['elapsed_seconds'])),flush=True)
        return
    dorsal=next(r for r in data['termGroups']if r['key']=='dorsum')
    original_faces=set(dorsal['surface']['faces'])
    incidence=defaultdict(list)
    for face_id,t in enumerate(native.faces):
        for k in range(3):incidence[tuple(sorted((int(t[k]),int(t[(k+1)%3]))))].append(face_id)
    adjacent=defaultdict(set)
    for ids in incidence.values():
        for face_id in ids:adjacent[face_id].update(set(ids)-{face_id})
    # A single, explicit source-mask closure fills only triangular pinholes:
    # every one of the missing triangle's three edge neighbours is already dorsal.
    added=sorted(i for i in range(len(native.faces))if i not in original_faces
                 and len(adjacent[i])==3 and adjacent[i]<=original_faces)
    source_faces=np.asarray(sorted(original_faces|set(added)),int)
    atomic_json(out/'dorsum-mask-regularization.json',dict(
        source_controls_sha256=digest_file(controls),original_source_faces=sorted(original_faces),
        added_source_faces=added,removed_source_faces=[],source_vertices_modified=False,
        method='One pass fills only literal source triangles whose three edge-neighbour triangles are already selected; no new surface or smoothing',
        anatomical_boundary_provisional=True))
    used,lookup=np.unique(native.faces[source_faces],return_inverse=True);faces=lookup.reshape(-1,3)
    fields=np.array([native.vertices[used]-.0007*normals[used],native.vertices[used]-.0025*normals[used]])
    p,surfaces,cells,wall_proof=layer_arrangement(fields,faces,['superior-longitudinal'])
    boundary=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]
                             for u in cells['superior-longitudinal']])
    sheet,sheet_check=emit(p,boundary,'superior-longitudinal-sheet',dict(
        kind='actual-dorsum-source-patch-finite-subdorsal-sheet',source_catalogue_sha256=digest_file(source),
        source_geometry_sha256=native.geometry_sha,source_dorsum_controls_sha256=digest_file(controls),
        source_face_ids=source_faces.tolist(),source_vertex_ids=used.tolist(),
        outer_depth_m=.0007,inner_depth_m=.0025,boundary='Actual open source-patch boundaries closed only by annular thickness rims',
        thickness_qualified='Illustrative source-normal offsets; no measured donor-specific muscle thickness',
        normal_source='literal stored Int16 source normals, normalized',anatomical_acceptance=False))
    atomic_npz(out/'sheet-correspondence.npz',fields=fields,faces=faces,source_face_ids=source_faces,source_vertex_ids=used)
    atomic_json(out/'sheet-check.json',sheet_check)
    passed=inner_check['passed'] and sheet_check['passed'] and len(sheet_check['topology']['components'])==1
    result=dict(schema='tongue-coherent-territory-preflight/v1',passed=passed,
        state='preflight-passed-ready-for-bounded-ownership' if passed else 'stopped-subdorsal-sheet-invalid',
        source_catalogue_sha256=digest_file(source),source_dorsum_controls_sha256=digest_file(controls),
        inner=inner_check,sheet=sheet_check,source_geometry_modified=False,installed=False,
        anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'preflight.json',result)
    review=dict(frame=native.frame,resolution=640,wall_seconds=60,
        meshes=[dict(mesh=str(out/'superior-longitudinal-sheet.npz'),name='Continuous SL source-bound sheet',color=[.72,.30,.29,1.]),
                dict(mesh=str(out/'native-inner-envelope.npz'),name='Candidate inner envelope context',color=[.6,.7,.8,.12])],
        views=[dict(name='top',direction=[0,1,0]),dict(name='right',direction=[1,0,0]),dict(name='front',direction=[0,0,1]),
               dict(name='section',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.625,.055],clear_outer=True))])
    atomic_json(out/'review-spec.json',review)
    print(json.dumps(dict(state=result['state'],sheet_components=len(sheet_check['topology']['components']),
                         elapsed_seconds=result['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--erosion-check',type=Path);args=ap.parse_args()
    out=args.output.resolve()
    if not out.is_relative_to(ROOT/'.local-models/geometry-workbench/tongue-completion'):
        raise ValueError('Keep proposals in local tongue workbench')
    preflight(out,args.erosion_check.resolve()if args.erosion_check else None)
