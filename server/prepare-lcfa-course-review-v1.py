"""Prepare bounded, literal-triangle review inputs; never modify the atlas."""
from pathlib import Path
import json, sys, time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh

B = ROOT / '.local-models/geometry-workbench/circulatory-completion/lateral-circumflex-femoral-native-v1'
INPUT = (ROOT / sys.argv[1]).resolve() if len(sys.argv) > 1 else B / 'bilateral-course-proof-v2'
OUT = (ROOT / sys.argv[2]).resolve() if len(sys.argv) > 2 else B / 'course-render-v1'

def main():
    start = time.monotonic()
    atlas_path = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(atlas_path) == 'b26860231a10b4b4382802751d6a8c15a7491e8baa83d629d73b93f5abb61707'
    atlas = json.loads(atlas_path.read_bytes()); parts = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader(); frame = 'model:local-ta98:assembled-metres'; specs = []
    for side in ['l', 'r']:
        z = np.load(INPUT / (side + '-course.npz'))
        course = Mesh(z['vertices'], z['faces'], frame, dict(normals=z['normal_words'], uvs=z['uvs']))
        used = np.unique(course.faces); low = course.vertices[used].min(0) - .008; high = course.vertices[used].max(0) + .008
        rows = []
        sources = [(course, 'Corrected LCFA course', [.70, .08, .06, 1.])]
        for pid, color in [
            ('TA98NATIVE:male:LCFA:ascending:' + side, [.95,.42,.08,1.]),
            ('ZA:Deep femoral artery.' + side, [.50,.10,.10,1.]),
            ('ZA:Femur.' + side + ':Bone', [.67,.64,.51,1.]),
            ('ZA:Femur.' + side + ':Cartilage', [.55,.64,.66,1.])
        ]:
            if pid not in parts:
                raise KeyError(pid)
            mesh = reader.read_part(atlas, parts[pid], frame=frame)
            tri = mesh.triangles
            keep = np.all(tri.max(1) >= low, axis=1) & np.all(tri.min(1) <= high, axis=1)
            if keep.any():
                mesh = Mesh(mesh.vertices, mesh.faces[keep], frame, mesh.attributes,
                            dict(mesh.source, diagnostic_context_face_subset=True, source_part_id=pid))
                sources.append((mesh, pid, color))
        for i, (mesh, name, color) in enumerate(sources):
            path = OUT / (side + '-' + str(i) + '.npz'); mesh.save(path)
            rows.append(dict(mesh=str(path), name=name, color=color))
        center = course.vertices[used].mean(0)
        spec = dict(frame=frame, meshes=rows, resolution=800, wall_seconds=40,
                    views=[dict(name=side+'-front', direction=[0,0,1]),
                           dict(name=side+'-posterior', direction=[0,0,-1]),
                           dict(name=side+'-oblique', direction=[1 if side=='l' else -1,.15,.4]),
                           dict(name=side+'-section', direction=[0,0,1],
                                section=dict(normal=[0,0,1], point=center.tolist(), clear_outer=True))])
        atomic_json(OUT / (side + '-review-spec.json'), spec)
        specs.append(str(OUT / (side + '-review-spec.json')))
    atomic_json(OUT / 'preparation-proof.json', dict(source_sha256=digest_file(atlas_path),
                actual_triangle_context_subsets=True, geometry_unchanged=True, specs=specs,
                elapsed_seconds=time.monotonic()-start))
    print(json.dumps(dict(specs=specs, elapsed_seconds=time.monotonic()-start)))

if __name__ == '__main__':
    main()
