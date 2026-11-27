"""Independently audit the actual blood classifier and a literal ICA guide."""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.exact_triangle_separation import exact_segment_triangle_intersection

t=time.monotonic()
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1'
D=B/'localized-no-rejoin-wall-v1'
path=D/'inner-union-classifier-diagnostic.npz';m=Mesh.load(path)
check=audit(m,require_closed=True,expected_components=1,expected_genus=[0])
assert check['passed'],check['failures']
profile_path=B/'ica-localized-no-rejoin-v1/world-profiles.npz'
with np.load(profile_path,allow_pickle=False)as z: points=z['world_profiles'].mean(1)
winding=igl.winding_number(m.vertices,m.faces,points)
assert np.all(winding>.999),winding
bvh=BVH(m);intersections=[];tested=0
for i,(a,b)in enumerate(zip(points,points[1:])):
    for face in bvh.query_bounds(np.minimum(a,b),np.maximum(a,b)):
        tested+=1
        hit=exact_segment_triangle_intersection(a,b,m.triangles[face])
        if hit['intersects']:intersections.append(dict(segment=i,face=int(face),predicate=hit))
assert not intersections,intersections
result=dict(schema='carotid-ICA-blood-guide-replay/v1',passed=True,
    inner_classifier_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    actual_profiles_sha256=hashlib.sha256(profile_path.read_bytes()).hexdigest(),
    blood_cast_audit=check,guide_points=len(points),complete_guide_segments=len(points)-1,
    winding_minimum=float(winding.min()),winding_maximum=float(winding.max()),
    exact_segment_triangle_pairs_tested=tested,guide_wall_intersections=intersections,
    scope='One connected genus-zero diagnostic blood cast and the complete literal ICA polygon-centre guide. External guard caps are construction-only. This is not whole circulation or physiological flow acceptance.',
    installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-t)
out=D/'independent-ICA-blood-guide.json';assert not out.exists()
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(passed=True,guide_segments=len(points)-1,exact_pairs=tested,elapsed_seconds=result['elapsed_seconds'])))
