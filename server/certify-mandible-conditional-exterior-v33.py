"""Exact exterior witnesses relative to literal triangles and defect domains.

This constructs no mandible, cap, signed-distance classifier, or anatomy.  The
only topological conclusion is conditional on a future bounded solid whose
entire boundary is contained in the literal source plus eight patch boxes.
Existence and anatomical validity of that solid are deliberately unestablished.
"""
from pathlib import Path
from collections import defaultdict
import itertools
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def sub(a, b):
    return tuple(x-y for x,y in zip(a,b))


def dot(a, b):
    return sum(x*y for x,y in zip(a,b))


def integer_points(points):
    """Lossless shared lattice for the actual finite Float64 coordinates."""
    ratios = [[float(x).as_integer_ratio() for x in p] for p in points]
    denominator = max(d for row in ratios for _,d in row)
    assert denominator & (denominator-1) == 0
    return [tuple(n*(denominator//d) for n,d in row) for row in ratios]


def exact_separation(first, second):
    """A strict separating axis is a sufficient exact disjointness proof.

    No negative result is interpreted as intersection: it is unresolved and
    rejected.  The extra in-plane axes cover disjoint coplanar triangles and
    segment/triangle cases without tolerance-based coplanarity decisions.
    """
    a = integer_points(np.concatenate([first, second]))
    p, q = a[:len(first)], a[len(first):]
    ep = [sub(p[(i+1)%len(p)], p[i]) for i in range(len(p))]
    eq = [sub(q[(i+1)%len(q)], q[i]) for i in range(len(q))]
    np_ = cross(ep[0], ep[1]) if len(p)==3 else (0,0,0)
    nq = cross(eq[0], eq[1]) if len(q)==3 else (0,0,0)
    axes = [(1,0,0),(0,1,0),(0,0,1),np_,nq]
    axes.extend(cross(x,y) for x in ep for y in eq)
    axes.extend(cross(n,x) for n in [np_,nq] for x in ep+eq)
    for k,axis in enumerate(axes):
        if axis == (0,0,0):
            continue
        dp, dq = [dot(axis,x) for x in p], [dot(axis,x) for x in q]
        if max(dp)<min(dq) or max(dq)<min(dp):
            return k
    return None


def overlapping_bvh_pairs(a,b):
    """Comparisons of literal min/max only: no rounded distance pruning."""
    stack = [(0,0)]
    while stack:
        x,y = stack.pop()
        if np.any(a.high[x]<b.low[y]) or np.any(b.high[y]<a.low[x]):
            continue
        if a.leaf(x) and b.leaf(y):
            for i in a.ids(x):
                t = a.triangles[i]
                for j in b.ids(y):
                    u = b.triangles[j]
                    if np.all(t.max(0)>=u.min(0)) and np.all(u.max(0)>=t.min(0)):
                        yield int(i),int(j)
        elif not a.leaf(x) and (b.leaf(y) or a.counts[x]>=b.counts[y]):
            stack.extend((int(c),y) for c in a.children[x])
        else:
            stack.extend((x,int(c)) for c in b.children[y])


def segment_box_separated(segment, box):
    """Exact analytic slab clipping; closed boxes, including their boundary."""
    from fractions import Fraction
    x = integer_points(np.concatenate([segment,box]))
    p,q,lo,hi = x
    begin, end = Fraction(0),Fraction(1)
    for k in range(3):
        d = q[k]-p[k]
        if d==0:
            if p[k]<lo[k] or p[k]>hi[k]:
                return True
            continue
        t0,t1 = Fraction(lo[k]-p[k],d),Fraction(hi[k]-p[k],d)
        if t0>t1:
            t0,t1=t1,t0
        begin,end=max(begin,t0),min(end,t1)
        if begin>end:
            return True
    return False


def main():
    began = time.monotonic()
    work = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17'
    out = work / 'mandible-conditional-exterior-v33'
    assert not out.exists()
    out.mkdir()
    previous_path = work/'mandible-source-defect-localization-v32/evidence.json'
    previous = json.loads(previous_path.read_bytes())
    expected = 'd73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
    binding = json.loads((work/'source-binding-v19/evidence.json').read_bytes())
    assert digest_file(ROOT/binding['current_catalogue']) == expected
    assert digest_file(ROOT/binding['snapshot']) == expected
    assert digest_file(ROOT/previous['artifact']) == previous['artifact_sha256']
    z = np.load(ROOT/previous['artifact'])
    source = Mesh(z['vertices'],z['faces'],previous['frame'])
    assert source.geometry_sha == previous['whole_joined_source_geometry_sha256']
    assert len(source.faces) == 18834
    bank_path = work/'surfaces-and-ancestry.npz'
    assert digest_file(bank_path) == previous['source_ancestry_bank_sha256']
    bank = np.load(bank_path)
    names = json.loads((work/'ports.json').read_bytes())['component_names']
    boxes = np.array([r['bounds_m'] for r in previous['defect_patches']])
    assert boxes.shape == (8,2,3)
    # These are literal patch bounds, conservative supersets of patch convex
    # hulls. No surface or defect alteration is implied by this domain choice.
    whole_bounds = np.array([np.minimum(source.vertices.min(0),boxes[:,0].min(0)),np.maximum(source.vertices.max(0),boxes[:,1].max(0))])
    tree = BVH(source)
    source_triangles = source.triangles
    tlo,thi = source_triangles.min(1),source_triangles.max(1)
    result,arrays=[],dict(defect_domain_boxes=boxes,source_vertices=source.vertices,source_faces=source.faces)
    directions=[(0,0,-1),(0,-1,0),(1,0,0),(-1,0,0),(0,1,0),(0,0,1)]
    directions += [x for x in itertools.product([-1,0,1],repeat=3) if x!=(0,0,0) and x not in directions]
    for ci,name in enumerate(names[1:],1):
        assert time.monotonic()-began<48, '48-second internal deadline'
        faces = bank['outer_faces'][bank['outer_component']==ci]
        used,inv = np.unique(bank['vertices'][faces].reshape(-1,3),axis=0,return_inverse=True)
        f = inv.reshape(-1,3)
        branch = Mesh(used,f,source.frame)
        assert branch.geometry_sha == previous['branch_checks'][ci-1]['actual_emitted_outer_geometry_sha256'] or np.array_equal(branch.triangles,bank['vertices'][faces])
        edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
        graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(used),len(used)))
        components,labels=connected_components(graph,directed=False)
        assert components==1
        branch_bounds=np.array([used.min(0),used.max(0)])
        assert np.all(np.any((branch_bounds[1]<boxes[:,0])|(boxes[:,1]<branch_bounds[0]),axis=1)), 'branch box overlaps defect domain'
        strict_axes=[]
        for box in boxes:
            axes=np.flatnonzero((branch_bounds[1]<box[0])|(box[1]<branch_bounds[0]))
            strict_axes.append(int(axes[0]))
        pairs,axis_ids=[],[]
        branch_triangles=branch.triangles
        for i,j in overlapping_bvh_pairs(BVH(branch),tree):
            assert len(pairs)<250000 and time.monotonic()-began<48
            axis=exact_separation(branch_triangles[i],source_triangles[j])
            assert axis is not None, (name,'source contact/unknown',i,j)
            pairs.append((i,j));axis_ids.append(axis)
        attempts=[];witness=None
        for direction in directions:
            d=np.asarray(direction,float)
            root_index=int(np.argmax(used@d))
            start=used[root_index]
            scale=max((whole_bounds[1,k]-start[k]+.01)/d[k] if d[k]>0 else (whole_bounds[0,k]-start[k]-.01)/d[k] if d[k]<0 else 0 for k in range(3))
            end=start+d*scale
            segment=np.array([start,end])
            box_clear=all(segment_box_separated(segment,b) for b in boxes)
            overlap=np.flatnonzero(np.all((thi>=segment.min(0))&(tlo<=segment.max(0)),axis=1))
            tested=[];axes=[];clear=box_clear
            if clear:
                for i in overlap:
                    axis=exact_separation(segment,source_triangles[i])
                    tested.append(int(i));axes.append(-1 if axis is None else axis)
                    if axis is None:
                        clear=False;break
            attempts.append(dict(direction=direction,source_vertex=root_index,box_clear=box_clear,source_triangle_candidates=int(len(overlap)),source_triangle_tests=len(tested),clear=clear))
            if clear:
                assert np.any((end<whole_bounds[0])|(end>whole_bounds[1]))
                witness=segment
                arrays[name+'-witness']=segment
                arrays[name+'-segment-source-candidates']=overlap
                arrays[name+'-segment-source-axes']=np.array(axes,int)
                break
        assert witness is not None, (name,'no proved escape witness')
        arrays[name+'-source-pairs']=np.array(pairs,int).reshape(-1,2)
        arrays[name+'-source-separating-axes']=np.array(axis_ids,int)
        result.append(dict(branch=name,source_bank_outer_face_count=len(faces),source_bank_outer_geometry_sha256=previous['branch_checks'][ci-1]['actual_emitted_outer_geometry_sha256'],exact_triangle_coordinate_sequence_preserved=True,outer_surface_connected_components=int(components),outer_surface_disjoint_from_all_18834_literal_source_faces=True,triangle_pair_candidates_exactly_separated=len(pairs),outer_bounds_strictly_separate_all_eight_defect_boxes=True,box_separating_axes=strict_axes,witness_m=witness.tolist(),attempts=attempts,conditional_exterior_certificate=True))
        atomic_json(out/'progress.json',dict(branches=result,elapsed_seconds=time.monotonic()-began),overwrite=True)
    artifact=out/'certificate-arrays.npz'
    atomic_npz(artifact,**arrays)
    proof=dict(schema='mandible-conditional-exterior/v33',source_catalogue_sha256=expected,registry_revision=previous['registry_revision'],candidate_geometry_sha256=previous['candidate_geometry_sha256'],source_bank_sha256=digest_file(bank_path),literal_source_geometry_sha256=source.geometry_sha,previous_evidence=str(previous_path.relative_to(ROOT)),previous_evidence_sha256=digest_file(previous_path),frame=source.frame,branches=result,arithmetic='Exact integers representing dyadic Float64 coordinates; strict separating projections and exact rational slab intervals. BVH rejection uses literal bounding comparisons only.',conditional_theorem='Let K be any bounded closed solid with boundary contained in the union of the literal mandible triangles and the eight closed defect-domain boxes. Each connected branch outer surface and its witnessed path is disjoint from that union. The endpoint lies beyond the whole union bounds and has a straight ray to infinity avoiding the union. Thus each entire connected branch outer surface belongs to the unbounded component of R^3 minus boundary(K), and cannot be in the interior of K.',explicit_unproved_premises=['Existence of a bounded clinically appropriate mandible solid K satisfying this boundary-domain restriction is not established.','Anatomical correctness of resolving all source defects within these boxes is not established.','This proof does not classify branch wall material or blood volume; it certifies the connected outer envelope only.','This proof establishes no source normals, wall thickness, or interface repair.'],conditional_domain_justification='Boxes contain all detected local open-edge, nonmanifold-edge, and invalid-link source face patches. This does not establish that every anatomically necessary source alteration or self-intersection correction is confined to them.',actual_mandible_occupied_classification='unknown-no-valid-literal-enclosure',full_anatomical_closure_assumption_justified=False,no_classifier_or_cap_created=True,no_native_changes=True,no_installation=True,complete=True,artifact=str(artifact.relative_to(ROOT)),artifact_sha256=digest_file(artifact),elapsed_seconds=time.monotonic()-began)
    assert digest_file(ROOT/binding['current_catalogue']) == expected
    assert Mesh.load(work/'finite-wall.npz').geometry_sha==previous['candidate_geometry_sha256']
    atomic_json(out/'evidence.json',proof)
    print(json.dumps(dict(output=str((out/'evidence.json').relative_to(ROOT)),branches=[dict(branch=r['branch'],conditional=r['conditional_exterior_certificate'],triangle_pairs=r['triangle_pair_candidates_exactly_separated'],attempts=len(r['attempts'])) for r in result],actual_mandible_occupied_classification=proof['actual_mandible_occupied_classification'],elapsed_seconds=proof['elapsed_seconds'])))


if __name__=='__main__':
    main()
