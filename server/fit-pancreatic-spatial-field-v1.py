"""One source-bound convex spatial correction; never installs a candidate."""
from pathlib import Path
import hashlib, json, sys, time
import numpy as np
from scipy import sparse
import clarabel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/pancreatic-public-cohort-v1/finite-network-construction-v1/coupled-spatial-field-v1'

def basis(points, centers, limits):
    p = np.asarray(points, float)
    difference = p[:, None] - centers
    exponent = -np.sum(difference**2, axis=2) / (2*.008**2)
    exponent -= exponent.max(axis=1, keepdims=True)
    phi = np.exp(exponent); phi /= phi.sum(axis=1, keepdims=True)
    grad_log = -difference / .008**2
    grad = phi[:, :, None] * (grad_log - np.sum(phi[:, :, None]*grad_log, axis=1)[:, None])
    lo, hi = limits; width = .008
    u = np.clip((p[:, 0]-lo)/width, 0, 1)
    v = np.clip((hi-p[:, 0])/width, 0, 1)
    left, right = u*u*(3-2*u), v*v*(3-2*v)
    envelope = left*right
    derivative = (np.where((u>0)&(u<1), 6*u*(1-u)/width, 0)*right
                  - left*np.where((v>0)&(v<1), 6*v*(1-v)/width, 0))
    grad = grad*envelope[:, None, None]
    grad[:, :, 0] += phi*derivative[:, None]
    return phi*envelope[:, None], grad

def main():
    began = time.monotonic(); plan_path = BASE/'plan.json'
    plan = json.loads(plan_path.read_text()); source = ROOT/plan['source_halfspaces']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == plan['source_halfspaces_sha256']
    data = json.loads(source.read_text()); rows = data['iterations'][0]['constraints']
    points = np.array([r['point'] for r in rows]); normals = np.array([r['normal'] for r in rows])
    rhs = np.array([r['minimum_dot_m'] for r in rows])*1000
    chosen = [int(np.argmax(rhs))]
    distance = np.linalg.norm(points-points[chosen[0]], axis=1)
    for _ in range(23):
        chosen.append(int(distance.argmax()))
        distance = np.minimum(distance, np.linalg.norm(points-points[chosen[-1]], axis=1))
    centers = points[chosen].copy(); windows = np.array(data['windows']); limits = windows[0, [0,3]]
    weights, _ = basis(points, centers, limits)
    matrix = (weights[:, :, None]*normals[:, None, :]).reshape(len(points), -1)
    count = len(centers); size = count*3
    # A connected, distance-weighted coefficient graph penalizes local strain.
    dd = np.linalg.norm(centers[:, None]-centers[None], axis=2)
    graph = np.exp(-dd**2/(2*.012**2)); np.fill_diagonal(graph, 0)
    laplacian = np.diag(graph.sum(axis=1))-graph
    quadratic = sparse.kron(sparse.eye(count)+.25*sparse.csc_matrix(laplacian), sparse.eye(3), format='csc')
    blocks = [-sparse.csc_matrix(matrix)]; bounds = [-rhs]; cones = [clarabel.NonnegativeConeT(len(rows))]
    for i in range(count):
        block = sparse.lil_matrix((4,size)); block[1:,3*i:3*i+3] = -np.eye(3)
        blocks.append(block.tocsc()); bounds.append(np.array([4.,0,0,0])); cones.append(clarabel.SecondOrderConeT(4))
    settings = clarabel.DefaultSettings(); settings.verbose = False; settings.max_iter = 80; settings.time_limit = 12.
    settings.tol_gap_abs = 1e-8; settings.tol_feas = 1e-8
    solver = clarabel.DefaultSolver(sparse.triu(quadratic).tocsc(), np.zeros(size), sparse.vstack(blocks,format='csc'), np.concatenate(bounds), cones, settings)
    result = solver.solve(); coefficients = np.asarray(result.x).reshape(count,3)
    slack = matrix@coefficients.ravel()-rhs
    passed = str(result.status)=='Solved' and slack.min()>=-1e-6 and np.linalg.norm(coefficients,axis=1).max()<=4.000001
    proof = dict(schema='pancreatic-coupled-spatial-field-fit/v1',plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest(),
        status=str(result.status),iterations=result.iterations,solver_seconds=result.solve_time,minimum_slack_mm=float(slack.min()),
        maximum_control_motion_mm=float(np.linalg.norm(coefficients,axis=1).max()),passed_linearized_problem=passed,
        constraints=len(rows),controls=count,mapped_geometry=False,installed=False,anatomical_acceptance=False)
    atomic_npz(BASE/'field.npz',centers=centers,coefficients_mm=coefficients,support_x=limits)
    if passed:
        anchor = np.load(ROOT/plan['anchor_contract'].replace('contract.json','literal-anchors.npz'))
        ap = np.r_[anchor['coeliac_side_inlet_points'],anchor['fixed_hilar_ownership_cut_points']]
        aw, ag = basis(ap,centers,limits)
        assert np.array_equal(aw,np.zeros_like(aw)) and np.array_equal(ag,np.zeros_like(ag))
        proof['protected_anchor_weights_and_gradients_exactly_zero'] = True
        for name in ['outer','blood','wall']:
            source_mesh = ROOT/plan['source_assembly']/(name+'.npz'); arrays=dict(np.load(source_mesh,allow_pickle=False))
            old = arrays['vertices']; ww, gg=basis(old,centers,limits)
            displacement = ww@coefficients/1000
            jacobian = np.eye(3)[None]+np.einsum('nki,kj->nji',gg,coefficients/1000)
            determinant=np.linalg.det(jacobian); singular=np.linalg.svd(jacobian,compute_uv=False)
            proof[name]=dict(source_sha256=hashlib.sha256(source_mesh.read_bytes()).hexdigest(),
                maximum_motion_mm=float(np.linalg.norm(displacement,axis=1).max()*1000),
                minimum_vertex_jacobian=float(determinant.min()),minimum_vertex_singular=float(singular.min()),maximum_vertex_singular=float(singular.max()),
                calibre_and_triangle_checks_pending=True)
            if determinant.min()<=0: raise ValueError('Nonpositive actual vertex Jacobian')
            arrays['vertices']=old+displacement; arrays['deformation_jacobians']=jacobian
            arrays['undeformed_vertices']=old; atomic_npz(BASE/(name+'-mapped.npz'),**arrays)
        proof['mapped_geometry']=True
    proof['elapsed_seconds']=time.monotonic()-began; atomic_json(BASE/'proof.json',proof)
    print(json.dumps(proof))

if __name__=='__main__': main()
