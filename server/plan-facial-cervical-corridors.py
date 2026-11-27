"""One capped transverse-corridor trial with source mouths/terminals fixed.

Adapted from the cystic v27 two-direction corridor.  Open3D only proposes
distances.  Construction and full actual-triangle context gates remain
mandatory.  No failed/disconnected course is silently replaced or installed.
"""
from pathlib import Path
import json
import sys
import time
import argparse
from copy import deepcopy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import igl
from scipy.optimize import lsq_linear
from scipy.spatial.distance import cdist
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.audit import topology
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.tube_sections import path_tangents, rotation_minimizing_frames
from geometry_tools.algorithms.tube_continuity import boundary_loops


def closed_classifier(mesh):
    gate = topology(mesh)
    if gate['closed_oriented_manifold']:
        return mesh
    if gate['nonmanifold_edges'] or gate['inconsistent_winding_edges'] or gate['invalid_vertex_links'] or gate['duplicate_faces'] or gate['zero_area_faces']:
        return None
    # Named native conduit ends are closed ONLY in this temporary occupancy
    # classifier, never in proposed tissue or viewer geometry.
    vertices, faces = mesh.vertices.copy(), mesh.faces.copy()
    try:
        loops = boundary_loops(faces)
    except ValueError:
        return None
    for loop in loops:
        center = len(vertices)
        vertices = np.r_[vertices, vertices[loop].mean(0)[None]]
        faces = np.r_[faces, np.array([[center,int(loop[(j+1)%len(loop)]),int(loop[j])] for j in range(len(loop))])]
    output = Mesh(vertices, faces, mesh.frame)
    return output if topology(output)['closed_oriented_manifold'] else None


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--version',default='v4')
    args=parser.parse_args()
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/circulatory-completion'
    parent_file = base / 'male-left-facial-cervical-routes-v2/routes.json'
    plan = json.loads(parent_file.read_bytes())
    proof_file = base / 'male-left-facial-cervical-finite-network-v1/feasibility-v1/evidence.json'
    first = json.loads(proof_file.read_bytes())
    directory = base / ('male-left-facial-cervical-routes-'+args.version)
    assert not directory.exists()
    directory.mkdir()
    source = ROOT / plan['source_catalogue']
    assert digest_file(source) == plan['source_catalogue_sha256']
    catalogue = json.loads(source.read_bytes())
    parts = {p['id']:p for p in catalogue['parts']}
    reader = AtlasReader()
    native = Mesh.load(ROOT / plan['fixed_parent']['mesh'])
    parent_classifier = closed_classifier(native)
    assert parent_classifier is not None
    blockers = {}
    for name in sorted(set(r['part_id'] for r in first['forbidden_crossings'])):
        mesh = reader.read_part(catalogue,parts[name],frame=native.frame,catalogue_path=source,catalogue_sha=plan['source_catalogue_sha256'])
        blockers[name] = (mesh, SurfaceIndex(mesh), closed_classifier(mesh))
    values = np.arange(-8,9)*.0005
    grid = np.array([[x,y] for x in values for y in values])
    zero = int(np.argmin(np.linalg.norm(grid,axis=1)))
    transitions = cdist(grid,grid)
    transition_cost = 2*(transitions/.0005)**2
    transition_cost[transitions > .0009] = np.inf
    result = deepcopy(plan)
    result['predecessor_route_contract'] = str(parent_file.relative_to(ROOT))
    result['predecessor_route_contract_sha256'] = digest_file(parent_file)
    summaries = []
    all_passed = True
    for route in result['routes']:
        points = np.load(ROOT / route['guide'])['points']
        n = len(points)
        radius = route['illustrative_radius_m']
        tangent = path_tangents(points)
        frames = rotation_minimizing_frames(points,tangent)
        lateral, transverse = frames[:,:,0], frames[:,:,1]
        candidates = points[:,None] + grid[None,:,0,None]*lateral[:,None] + grid[None,:,1,None]*transverse[:,None]
        flat = candidates.reshape(-1,3)
        parent_distance = np.sqrt(igl.point_mesh_squared_distance(flat,parent_classifier.vertices,parent_classifier.faces)[0]).reshape(n,-1)
        parent_inside = (np.abs(igl.winding_number(parent_classifier.vertices,parent_classifier.faces,flat))>.99).reshape(n,-1)
        complete_inside_parent = parent_inside & (parent_distance > radius+.00002)
        safe = np.ones((n,len(grid)),bool)
        per_blocker = []
        for name,(mesh,index,classifier) in blockers.items():
            if name in {'ZA:Superior pharyngeal constrictor.l'} and route['name']=='tonsillar':
                continue  # Explicit original tissue-perforating interface.
            if time.monotonic()-began > 38:
                atomic_json(directory/'corridor-evidence.json',dict(stopped_before_construction=True,reason='38 s proposal-work deadline',completed_routes=summaries,elapsed_seconds=time.monotonic()-began,installed=False))
                print(json.dumps(dict(stopped_before_construction=True,reason='38 s proposal-work deadline',elapsed_seconds=time.monotonic()-began)))
                return
            closest,ev=index.closest_points(flat,certify=False)
            distance=closest['distance_m'].reshape(n,-1)
            valid=distance > radius+.00012
            if classifier is not None:
                inside=(abs(igl.winding_number(classifier.vertices,classifier.faces,flat))>.99).reshape(n,-1)
                valid &= ~inside
            safe &= valid | complete_inside_parent
            per_blocker.append(dict(part_id=name,valid_candidates=int(valid.sum()),temporary_occupied_classifier=classifier is not None,query=ev))
        safe[:3]=False
        safe[:3,zero]=True
        safe[-2:]=False
        safe[-2:,zero]=True
        # Fixed recipient centers still require explicit obstacle clearance:
        # they are not exempt simply because the initial route used them.
        endpoint_ok=True
        for name,(mesh,index,classifier) in blockers.items():
            if name=='ZA:Superior pharyngeal constrictor.l' and route['name']=='tonsillar':continue
            d=np.sqrt(igl.point_mesh_squared_distance(points[-2:],mesh.vertices,mesh.faces)[0])
            occupied=abs(igl.winding_number(classifier.vertices,classifier.faces,points[-2:]))>.99 if classifier is not None else np.zeros(2,bool)
            if np.any(d <= radius+.00012) or occupied.any():endpoint_ok=False
        cost=np.full(len(grid),np.inf)
        cost[zero]=0
        previous=np.full((n,len(grid)),-1,int)
        failed=None
        for j in range(3,n):
            choices=cost[:,None]+transition_cost
            previous[j]=choices.argmin(0)
            cost=choices.min(0)+.02*np.sum((grid/.001)**2,axis=1)
            cost[~safe[j]]=np.inf
            if not np.isfinite(cost).any():failed=j;break
        if failed is not None or not endpoint_ok:
            all_passed=False
            summaries.append(dict(name=route['name'],feasible=False,first_disconnected_station=failed,fixed_endpoint_avoids_forbidden_occupied_tissue=endpoint_ok,feasible_candidates_per_station=safe.sum(1).tolist(),blockers=per_blocker))
            continue
        selected=np.full(n,zero,int)
        for j in range(n-1,3,-1):selected[j-1]=previous[j,selected[j]]
        target=candidates[np.arange(n),selected]-points
        variable=np.arange(3,n-2)
        second=np.diff(np.eye(n),n=2,axis=0)
        A=np.r_[np.eye(n),8*second]
        displacement=np.zeros_like(points)
        solved=[]
        for axis in range(3):
            rhs=np.r_[target[:,axis]*1000,np.zeros(n-2)]
            goal=target[variable,axis]*1000
            fit=lsq_linear(A[:,variable],rhs,bounds=(goal-.08,goal+.08),max_iter=60,tol=1e-9)
            displacement[variable,axis]=fit.x/1000
            solved.append(bool(fit.success))
        corrected=points+displacement
        path=directory/(route['name']+'-guide.npz')
        atomic_npz(path,points=corrected,outer_radii=np.full(n,radius),inner_radii=np.full(n,radius-route['illustrative_wall_m']),source_guide=points,station_displacement=displacement,frame=np.array(plan['frame']))
        route.update(guide=str(path.relative_to(ROOT)),guide_sha256=digest_file(path),length_m=float(np.linalg.norm(np.diff(corrected,axis=0),axis=1).sum()))
        summaries.append(dict(name=route['name'],feasible=True,maximum_center_shift_m=float(np.linalg.norm(displacement,axis=1).max()),source_first_three_and_final_two_centers_literal=True,calibre_changed=False,smoothing_success=all(solved),blockers=per_blocker))
    proof=dict(schema='facial-cervical-two-direction-corridor/v1',source_catalogue_sha256=plan['source_catalogue_sha256'],previous_feasibility_sha256=digest_file(proof_file),
               lattice_spacing_m=.0005,maximum_offset_per_axis_m=.004,maximum_station_offset_step_m=.0009,completed_routes=summaries,all_routes_feasible=all_passed,
               full_actual_triangles_and_field_prisms_still_require_recheck=True,private_source_copied=False,installed=False,eligible_for_installation=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(directory/'corridor-evidence.json',proof)
    if all_passed:
        result['corridor_evidence_sha256']=digest_file(directory/'corridor-evidence.json')
        atomic_json(directory/'routes.json',result)
    print(json.dumps(dict(routes=[{k:r.get(k) for k in ['name','feasible','first_disconnected_station','fixed_endpoint_avoids_forbidden_occupied_tissue','maximum_center_shift_m']} for r in summaries],all_routes_feasible=all_passed,elapsed_seconds=proof['elapsed_seconds'])))


if __name__=='__main__':main()
