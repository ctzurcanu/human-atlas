"""Independent emitted-mesh checks; a numerical pass never accepts fragmented anatomy."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import diagnostics


def main(candidate):
    began=time.monotonic()
    report=json.loads((candidate/'construction.json').read_text())
    inputs=json.loads((candidate/'inputs.json').read_text())
    source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(source)==report['source_catalogue_sha256']
    work=candidate.parent
    native=Mesh.load(work/'source-v1/meshes/native-00.npz')
    retained=[]
    for row in inputs['cutters']:
        assert digest_file(row['path'])==row['sha256']
        retained.append((row['key'],Mesh.load(row['path'])))
    rows=[];meshes=[];all_pairs=[]
    for row in report['meshes']:
        assert digest_file(row['mesh'])==row['sha256']
        mesh=Mesh.load(row['mesh']);topo=topology(mesh)
        q=np.concatenate([mesh.vertices,mesh.triangles.mean(1),
            (mesh.triangles[:,0]+mesh.triangles[:,1])/2,(mesh.triangles[:,1]+mesh.triangles[:,2])/2,
            (mesh.triangles[:,2]+mesh.triangles[:,0])/2])
        winding=np.abs(igl.winding_number(native.vertices,native.faces,q))
        _,outer=pair_contacts(native,mesh,0.,max_pairs=100000,max_tests=1000000)
        cpp=diagnostics(mesh)
        checks=dict(closed_oriented=topo['closed_oriented_manifold'],
            all_queries_contained=bool(np.all(winding>.99)),
            envelope_crossings_zero=outer['crossing_face_pairs']==0 and outer['coplanar_overlap_face_pairs']==0,
            source_budget_complete=outer['complete_within_pair_budget'],self_intersections_zero=cpp['self_intersection_pairs']==0)
        rows.append(dict(key=row['key'],geometry_sha256=mesh.geometry_sha,topology=topo,
            independent_queries=len(q),minimum_absolute_winding=float(winding.min()),
            envelope_contacts=outer,cpp_diagnostics=cpp,checks=checks,numerical_checks_passed=all(checks.values())))
        meshes.append((row['key'],mesh))
        for name,other in retained:
            arrays,summary=pair_contacts(mesh,other,0.,max_pairs=50000,max_tests=1000000)
            bad=np.isin(arrays['relation'],['crossing','coplanar-overlap'])
            witness=arrays['pairs'][bad][:10].tolist()
            distances=igl.point_mesh_squared_distance(q,other.vertices,other.faces)[0]
            all_pairs.append(dict(first=row['key'],second=name,kind='retained-occupied',**summary,
                defect_face_pairs=witness,minimum_sampled_distance_m=float(np.sqrt(distances).min())))
    for i,(name,mesh) in enumerate(meshes):
        for other_name,other in meshes[:i]:
            arrays,summary=pair_contacts(mesh,other,0.,max_pairs=50000,max_tests=1000000)
            all_pairs.append(dict(first=name,second=other_name,kind='cross-intrinsic',**summary))
    pair_pass=all(r['complete_within_pair_budget'] and r['crossing_face_pairs']==0 and r['coplanar_overlap_face_pairs']==0 for r in all_pairs)
    passed=all(r['numerical_checks_passed']for r in rows) and pair_pass
    result=dict(schema='tongue-intrinsic-ownership-emitted-audit/v1',passed=passed,
        construction_sha256=digest_file(candidate/'construction.json'),source_catalogue_sha256=digest_file(source),
        rows=rows,contact_pairs=all_pairs,total_pair_count=len(all_pairs),expected_pair_count=63,
        actual_crossing_face_pairs=sum(r['crossing_face_pairs']for r in all_pairs),
        actual_coplanar_overlap_face_pairs=sum(r['coplanar_overlap_face_pairs']for r in all_pairs),
        preservation_gates_passed=report['numerical_preservation_gates_passed'],
        eligible_for_anatomical_review=passed and report['numerical_preservation_gates_passed'],
        anatomical_acceptance=False,installed=False,source_geometry_modified=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(candidate/'checked-geometry.json',result)
    print(json.dumps({k:result[k]for k in ['passed','total_pair_count','actual_crossing_face_pairs',
        'actual_coplanar_overlap_face_pairs','preservation_gates_passed','elapsed_seconds']}),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True,type=Path)
    args=ap.parse_args();main(args.candidate.resolve())
