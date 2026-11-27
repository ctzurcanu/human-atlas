"""Check actual emitted root-to-atrium paths, including complete segments."""
from pathlib import Path
from collections import deque
from types import SimpleNamespace
import argparse,json,sys,time
import numpy as np
import igl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.tissue import decode_packed
from geometry_tools.spatial import BVH
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.exact_triangle_separation import exact_segment_triangle_intersection
ap=argparse.ArgumentParser();ap.add_argument('--bank',required=True);args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.bank;e=json.loads((src/'emission.json').read_text());assert digest_file(Path(e['bank']))==e['bank_sha256']
cells=decode_packed(SimpleNamespace(data=Path(e['bank']).read_bytes(),cells=e['parts']),e['frame'],np.array(e['origin']));meshes=list(cells.values());v=meshes[0].vertices
assert all(np.array_equal(v,m.vertices)for m in meshes),'One literal shared position bank required'
f=np.concatenate([m.faces for m in meshes]);_,first,count=np.unique(np.sort(f,axis=1),axis=0,return_index=True,return_counts=True);assert count.max()==2
f=f[first[count==1]];mesh=Mesh(v,f,e['frame']);tree=BVH(mesh);left=Mesh.load(h/'blood-spaces-v130/fluid-boundary-4.npz');courses=json.loads((h/'native-pv-mouths-v66/source-courses.json').read_text())['rings'];rows=[];tests=0
for course in courses:
    p=np.array(course['points_m']);faces=np.array(course['faces'],int);loops=boundary_loops(faces);assert len(loops)==2;adj=[set()for _ in p]
    for tri in faces:
        for a,b in zip(tri,np.roll(tri,1)):adj[a].add(int(b));adj[b].add(int(a))
    level=np.full(len(p),-1);level[loops[0]]=0;queue=deque(map(int,loops[0]))
    while queue:
        a=queue.popleft()
        for b in adj[a]:
            if level[b]<0:level[b]=level[a]+1;queue.append(b)
    assert np.all(level>=0);curve=np.array([p[level==j].mean(0)for j in range(level.max()+1)])
    endpoints=np.abs(igl.winding_number(left.vertices,left.faces,curve[[0,-1]]));inside=np.abs(igl.winding_number(v,f,curve));hits=[]
    for i,(a,b)in enumerate(zip(curve,curve[1:])):
        near=tree.query_bounds(np.minimum(a,b),np.maximum(a,b),distance_m=0.)
        tests+=len(near)
        if tests>100000 or time.monotonic()-start>35:raise TimeoutError('Finite cardiac path query budget')
        for face in near:
            if exact_segment_triangle_intersection(a,b,v[f[face]])['intersects']:hits.append(dict(segment=i,face=int(face)))
    row=dict(part_id=course['part_id'],component=course['component'],stations=len(curve),
        starts_or_ends_in_left_blood_space=bool((endpoints>.5).any()),endpoint_left_blood_winding=endpoints.tolist(),
        all_station_centers_outside_tissue=bool(np.all(inside<.5)),tissue_crossings=hits,
        passed=bool((endpoints>.5).any()and np.all(inside<.5)and not hits),
        curve_m=curve.tolist());rows.append(row)
result=dict(passed=all(r['passed']for r in rows),bank_sha256=e['bank_sha256'],actual_boundary_geometry_sha256=mesh.geometry_sha,
    complete_segment_tests=True,curves=rows,triangle_tests=tests,elapsed_seconds=time.monotonic()-start,
    diagnostic_only=True,anatomical_valve_patency_established=False)
atomic_json(src/'root-fluid-paths.json',result);print(json.dumps(dict(passed=result['passed'],paths=len(rows),crossings=sum(len(r['tissue_crossings'])for r in rows),all_reach_left_blood_space=all(r['starts_or_ends_in_left_blood_space']for r in rows),elapsed_seconds=time.monotonic()-start)),flush=True)
raise SystemExit(0 if result['passed']else 2)
