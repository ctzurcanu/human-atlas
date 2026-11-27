"""Embed bladder ridges in the existing finite wall, with fixed outer tissue.

The source-defined ostia are real boundaries of the retained reservoir lining.
They are not the arbitrary coat-ownership seams in a joining collar. Relief is
a single monotone height deformation of the entire shared tissue bank; no
independent ridge solid or lumen disk is added. Dimensions are review inputs,
not measurements of this donor.
"""
from pathlib import Path
import json

import numpy as np

from .algorithms.annular_pipe_envelopes import boundary_loops
from .algorithms.finite_tissue_boolean import CommonPositionBank
from .algorithms.shared_triangle_edges import bisect_surfaces
from .common import atomic_json, atomic_npz, digest_file
from .distance import point_triangle
from .model import AtlasReader, Mesh
from .spatial import BVH


def smoothstep(t):
    t = np.clip(t, 0., 1.)
    return t*t*t*(10.+t*(-15.+6.*t))


def back_graph(vertices, faces, xy):
    """Exact Float64 triangle witnesses on the posterior graph over X/Y.

    Return NaN outside the supplied physical sheet. In particular an actual
    opening is never completed by an extrapolated plane or a convex hull.
    """
    xy = np.asarray(xy, float)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
        raise ValueError('Finite X/Y points in metres are required')
    triangles = np.asarray(vertices, float)[faces]
    a = triangles[:, 0, :2]
    e, f = triangles[:, 1, :2]-a, triangles[:, 2, :2]-a
    determinant = e[:, 0]*f[:, 1]-e[:, 1]*f[:, 0]
    keep = np.abs(determinant) > 1e-22
    ids = np.flatnonzero(keep)
    a, e, f, determinant, triangles = a[keep], e[keep], f[keep], determinant[keep], triangles[keep]
    height = np.full(len(xy), np.nan); parents = np.full(len(xy), -1, int)
    for start in range(0, len(xy), 256):
        d = xy[start:start+256, None]-a
        u = (d[..., 0]*f[:, 1]-d[..., 1]*f[:, 0])/determinant
        v = (e[:, 0]*d[..., 1]-e[:, 1]*d[..., 0])/determinant
        inside = (u >= -1e-11) & (v >= -1e-11) & (u+v <= 1.+1e-11)
        z = triangles[:, 0, 2]+u*(triangles[:, 1, 2]-triangles[:, 0, 2])+v*(triangles[:, 2, 2]-triangles[:, 0, 2])
        z[~inside] = np.inf
        nearest = z.argmin(1); rows = np.arange(len(d)); values = z[rows, nearest]
        hit = np.isfinite(values)
        height[start+rows[hit]] = values[hit]
        parents[start+rows[hit]] = ids[nearest[hit]]
    return height, parents


def inherited_body_rims(whole, inner_ids, body, body_inner_faces, tolerance_m=2e-12):
    """Recover the actual reservoir/collar seams from literal source facets."""
    lining = whole.faces[inner_ids]
    body_lining = Mesh(body.vertices, body_inner_faces, body.frame)
    tree = BVH(body_lining)
    result, evidence = tree.closest_points(whole.vertices[lining].mean(1), backend='auto')
    parents = result['face_ids']
    distances = [point_triangle(whole.vertices[lining[:, k]], body_lining.triangles[parents])[0] for k in range(3)]
    inherited = np.max(distances, axis=0) <= tolerance_m**2
    rims = boundary_loops(lining[inherited])
    if len(rims) != 3:
        raise ValueError('The actual inherited bladder lining must have three distinct ostial boundaries')
    neck = min(rims, key=lambda q: whole.vertices[q].mean(0)[1])
    remaining = [q for q in rims if q is not neck]
    left = max(remaining, key=lambda q: whole.vertices[q].mean(0)[0])
    right = min(remaining, key=lambda q: whole.vertices[q].mean(0)[0])
    return dict(neck=neck, left=left, right=right), dict(
        source_triangles_replayed=int(inherited.sum()), tolerance_m=tolerance_m,
        actual_triangle_witnesses=evidence, ownership_boundary_used_as_ostium=False)


class TrigonalRelief:
    def __init__(self, body, body_inner_faces, body_outer_faces, rims, *,
                 crest_height_m=.001, crest_half_width_m=.0012,
                 uvula_height_m=.0008, falloff_into_fluid_m=.005):
        if not 0 < crest_height_m <= .002 or not 0 < uvula_height_m <= .0015:
            raise ValueError('Bounded explicit relief heights in metres are required')
        if not .0006 <= crest_half_width_m <= .003:
            raise ValueError('An explicit bounded crest half-width is required')
        maximum = crest_height_m+uvula_height_m
        # The two supports are separate here, but this conservative sum also
        # bounds their possible overlap. It guarantees monotonicity in Z.
        if falloff_into_fluid_m <= 1.875*maximum:
            raise ValueError('Fluid-side relief must remain strictly monotone')
        self.body, self.inner, self.outer = body, body_inner_faces, body_outer_faces
        self.rims = {name:np.asarray(q, float) for name,q in rims.items()}
        left, right, neck = [self.rims[name].mean(0) for name in ('left', 'right', 'neck')]
        self.low, self.high = float(right[0]), float(left[0])
        if self.high-self.low < .015 or min(left[1], right[1])-neck[1] < .015:
            raise ValueError('The explicit source ostia do not define a usable trigone')
        self.neck = neck
        self.crest_y = (float(right[1]), float(left[1]))
        self.crest_height, self.crest_width = float(crest_height_m), float(crest_half_width_m)
        self.uvula_height, self.falloff = float(uvula_height_m), float(falloff_into_fluid_m)
        self.uvula_center = np.array([neck[0], neck[1]+.0065])
        self.uvula_radii = np.array([.003, .0045])
        self.protected = []
        for name in ('left', 'right', 'neck'):
            q = self.rims[name]; center = q[:, :2].mean(0)
            radius = float(np.linalg.norm(q[:, :2]-center, axis=1).max())+.0002
            self.protected.append((center, radius))

    def amplitudes(self, xy):
        xy = np.asarray(xy, float)
        fraction = (xy[:, 0]-self.low)/(self.high-self.low)
        y = self.crest_y[0]+fraction*(self.crest_y[1]-self.crest_y[0])
        crest = self.crest_height*smoothstep(1.-np.abs(xy[:, 1]-y)/self.crest_width)
        crest *= smoothstep((xy[:, 0]-self.low)/.004)*smoothstep((self.high-xy[:, 0])/.004)
        r = (xy-self.uvula_center)/self.uvula_radii
        uvula = self.uvula_height*smoothstep(1.-r[:, 0]**2)*smoothstep(1.-r[:, 1]**2)
        protect = np.ones(len(xy))
        for center, radius in self.protected:
            protect *= smoothstep((np.linalg.norm(xy-center, axis=1)-radius)/.001)
        return np.c_[crest*protect, uvula*protect]

    def displacement(self, positions):
        positions = np.asarray(positions, float)
        amplitude = self.amplitudes(positions[:, :2]).sum(1)
        ids = np.flatnonzero(amplitude > 0)
        delta = np.zeros(len(positions)); derivative_lower_bound = np.ones(len(positions))
        if not len(ids):
            return delta, derivative_lower_bound
        q = positions[ids]
        inner, _ = back_graph(self.body.vertices, self.inner, q[:, :2])
        outer, _ = back_graph(self.body.vertices, self.outer, q[:, :2])
        valid = np.isfinite(inner) & np.isfinite(outer) & (inner-outer > 1e-6)
        if np.any(~valid & (q[:, 2] < self.body.bounds[1, 2])):
            # A feature footprint crossing a hole needs a different reviewed
            # support, not a fabricated surface across that opening.
            raise ValueError('Relief support leaves the actual posterior inner/outer wall graph')
        ids, q, inner, outer = ids[valid], q[valid], inner[valid], outer[valid]
        depth = np.clip((q[:, 2]-outer)/(inner-outer), 0., 1.)
        fluid = smoothstep(1.-(q[:, 2]-inner)/self.falloff)
        delta[ids] = amplitude[ids]*depth*fluid
        derivative_lower_bound[ids] = 1.-1.875*amplitude[ids]/self.falloff
        return delta, derivative_lower_bound

    def candidate_faces(self, positions, faces):
        """Conservative feature support bounds; never miss a narrow ridge."""
        q = positions[faces]; low, high = q.min(1), q.max(1)
        cy = np.array(self.crest_y)
        regions = [(np.array([self.low, cy.min()-self.crest_width]),
                    np.array([self.high, cy.max()+self.crest_width])),
                   (self.uvula_center-self.uvula_radii, self.uvula_center+self.uvula_radii)]
        selected = np.zeros(len(faces), bool)
        for a,b in regions:
            selected |= np.all(high[:, :2] >= a, axis=1) & np.all(low[:, :2] <= b, axis=1)
        # The feature is on the posterior lining; exclude distant front wall.
        selected &= low[:, 2] < -.01
        return selected

    def specification(self):
        return dict(crest_height_m=self.crest_height, crest_half_width_m=self.crest_width,
            crest_endpoint_y_m=list(self.crest_y), uvula_height_m=self.uvula_height,
            uvula_center_xy_m=self.uvula_center.tolist(), uvula_half_extent_xy_m=self.uvula_radii.tolist(),
            fluid_falloff_m=self.falloff, maximum_height_derivative=1.875,
            minimum_global_z_jacobian_bound=1.-1.875*(self.crest_height+self.uvula_height)/self.falloff,
            dimensions_are_review_inputs=True, donor_histometry_measured=False,
            references=['https://commons.wikimedia.org/wiki/File:Gray1140.png',
                        'https://pubmed.ncbi.nlm.nih.gov/32449124/'],
            mri_ap_diameter_is_not_assumed_to_equal_ridge_height=True)


def build(delivery_dir, coordinates_dir, body_carriers, catalogue, output, *, maximum_edge_m=.0009):
    """Repair the complete source envelope before rebuilding its partitions."""
    from .audit import audit
    from .bladder_carriers import recover
    delivery_dir,coordinates_dir,body_carriers,catalogue,output=map(lambda p:Path(p).resolve(),
        (delivery_dir,coordinates_dir,body_carriers,catalogue,output))
    if output.exists():
        raise FileExistsError(output)
    if not .0004 <= maximum_edge_m <= .0015:
        raise ValueError('Use a bounded local refinement spacing in metres')
    d=json.loads((delivery_dir/'delivery.json').read_text())
    info=json.loads((coordinates_dir/'coordinates.json').read_text())
    atlas=json.loads(catalogue.read_text());reader=AtlasReader();byid={p['id']:p for p in atlas['parts']}
    verified=[]
    for group,row in d['groups'].items():
        for pid in row['cells']:
            actual=reader.read_part(atlas,byid[pid],frame=d['frame'])
            source=Mesh.load(delivery_dir/group/(pid.replace(':','_')+'.npz'))
            if actual.surface_sha!=source.surface_sha:
                raise ValueError('The installed source tissue changed: '+pid)
            verified.append(pid)
    whole=Mesh.load(info['wall']);body=Mesh.load(info['body_wall'])
    if digest_file(info['wall'])!=info['wall_sha256']or digest_file(info['body_wall'])!=info['body_wall_sha256']:
        raise ValueError('The pre-partition source wall changed')
    with np.load(body_carriers,allow_pickle=False)as z:
        if not np.array_equal(body.faces,z['wall_faces']):
            raise ValueError('Read landmark coordinates from the fitted body with the original source connectivity')
        bi=z['inner_side_faces'].copy();bo=z['outer_and_rim_faces'].copy()
    with np.load(coordinates_dir/'coordinates.npz',allow_pickle=False)as z:
        inner_ids=z['source_inner_face_ids'].copy()
    clinical_rims,ancestry=inherited_body_rims(whole,inner_ids,body,bi)
    points={name:whole.vertices[ring]for name,ring in clinical_rims.items()}
    relief=TrigonalRelief(body,bi,bo,points)
    output.mkdir(parents=True)
    plan=dict(schema='bladder-source-envelope-relief/v2',catalogue=str(catalogue),catalogue_sha256=digest_file(catalogue),
        installed_source_cells_replayed=verified,source_wall_sha256=digest_file(info['wall']),
        body_sha256=digest_file(info['body_wall']),source_delivery_sha256=digest_file(delivery_dir/'delivery.json'),
        actual_ostia={k:q.tolist()for k,q in points.items()},ostial_ancestry=ancestry,relief=relief.specification(),
        maximum_edge_m=maximum_edge_m,outer_wall_fixed=True,renal_roots_fixed=True,skeleton_fixed=True,
        partitions_rebuilt_after_whole_surface_repair=True,anatomical_acceptance=False,installed=False)
    atomic_json(output/'plan.json',plan)
    rest=np.setdiff1d(np.arange(len(whole.faces)),inner_ids)
    surfaces=dict(inner=whole.faces[inner_ids],outer_and_rims=whole.faces[rest])
    bank=CommonPositionBank(whole.vertices);rounds=[]
    for iteration in range(12):
        v=bank.positions;f=surfaces['inner'];selected=f[relief.candidate_faces(v,f)];edges=set()
        for a,b in ((0,1),(1,2),(2,0)):
            e=selected[:,[a,b]];long=np.linalg.norm(v[e[:,0]]-v[e[:,1]],axis=1)>maximum_edge_m
            edges.update(tuple(sorted(map(int,q)))for q in e[long])
        if not edges:
            break
        midpoints={e:bank.insert((v[e[0]]+v[e[1]])*.5)for e in sorted(edges)}
        surfaces,_,count=bisect_surfaces(surfaces,midpoints)
        rounds.append(dict(iteration=iteration,split_edges=len(edges),changed_faces=count))
    else:
        raise ValueError('Local source-sheet refinement exceeded its work bound')
    v=bank.positions;delta,bound=relief.displacement(v)
    fixed=np.unique(surfaces['outer_and_rims']);delta[fixed]=0.
    moved=v.copy();moved[:,2]+=delta
    for name,q in points.items():
        if np.any(relief.displacement(q)[0]!=0):
            raise ValueError('A source ostial curve moved: '+name)
    f=np.r_[surfaces['inner'],surfaces['outer_and_rims']]
    mesh=Mesh(moved,f,whole.frame);mesh.save(output/'urinary-wall.npz')
    check=audit(mesh,require_closed=True,expected_components=1,expected_genus=[2],max_pairs=1000000)
    atomic_json(output/'whole-check.json',check)
    if not check['passed']:
        raise ValueError('The repaired whole source wall fails the geometry gates')
    carriers=recover(mesh,output/'verified-carriers',inner_face_ids=np.arange(len(surfaces['inner'])),
        mouth_labels=['external-urethral-orifice','right-renal-root','left-renal-root'],cap_forbid_source_diagonals=True)
    result=dict(**plan,passed=True,whole_check=check,carriers=carriers,
        refinement=rounds,source_controls=len(whole.vertices),refined_controls=len(moved),
        moved_controls=int(np.count_nonzero(delta)),maximum_displacement_m=float(delta.max()),
        minimum_continuous_z_jacobian_bound=float(bound.min()),
        ostial_curves_literal=True,outer_source_controls_literal=bool(np.array_equal(moved[fixed],v[fixed])),
        tissue_partitions_pending=True,normal_fields_pending=True,rendered_review_pending=True,
        source_wall_geometry_sha256=mesh.geometry_sha,one_connected_urinary_fluid_domain=True)
    atomic_json(output/'outlet-with-carriers.json',result)
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('delivery-dir','coordinates-dir','body-carriers','catalogue','output'):
        p.add_argument('--'+flag,required=True)
    p.add_argument('--maximum-edge-m',type=float,default=.0009)
    r=build(**vars(p.parse_args()))
    print(json.dumps({k:r[k]for k in ('passed','source_controls','refined_controls','maximum_displacement_m')},indent=2))
