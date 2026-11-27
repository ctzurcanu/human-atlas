"""Try one witnessed source diagonal repair without moving any control point."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json, digest_file

BASE = ROOT / '.local-models/geometry-workbench/tongue-completion'

def main(output):
    started = time.monotonic()
    output = output.resolve()
    if output.exists():
        raise FileExistsError(output)
    source_audit = json.loads((BASE / 'source-v1/source-audit.json').read_text())
    source = source_audit['sourceCatalogues']['male']
    if digest_file(source['path']) != source['sha256']:
        raise ValueError('Installed source revision changed')
    output.mkdir(parents=True)
    rows = []
    for side, native in [('left', 'native-01.npz'), ('right', 'native-02.npz')]:
        path = BASE / 'source-v1/meshes' / native
        old = Mesh.load(path)
        faces = old.faces.copy()
        first, second = faces[3009].tolist(), faces[3014].tolist()
        shared = set(first) & set(second)
        if len(shared) != 2:
            raise ValueError('Witnessed opposite edge changed')
        boundary = {}
        for face in [first, second]:
            for a, b in zip(face, face[1:] + face[:1]):
                if {a, b} != shared:
                    boundary[a] = b
        start = next(v for v in first if v not in shared)
        cycle = [start]
        for _ in range(3):
            cycle.append(boundary[cycle[-1]])
        a, b, c, d = cycle
        if boundary[d] != a or any({a, c}.issubset(face) for face in faces):
            raise ValueError('Replacement diagonal already belongs to this surface')
        faces[3009], faces[3014] = [a, b, c], [a, c, d]
        candidate = Mesh(old.vertices.copy(), faces, old.frame,
                         {k: v.copy() for k, v in old.attributes.items()},
                         source=dict(kind='witnessed-source-diagonal-repair',
                                     source_path=str(path), source_sha256=digest_file(path),
                                     original_geometry_sha256=old.geometry_sha,
                                     changed_face_ids=[3009, 3014],
                                     coordinates_and_corner_fields_unchanged=True))
        target = output / f'genioglossus-{side}.npz'
        candidate.save(target)
        emitted = Mesh.load(target)
        result = audit(emitted, require_closed=True, max_pairs=250000, max_tests=2000000)
        row = dict(side=side, mesh=str(target), source_sha256=digest_file(path),
                   source_geometry_sha256=old.geometry_sha,
                   geometry_sha256=emitted.geometry_sha, sha256=digest_file(target),
                   changed_faces=[3009, 3014], moved_vertices=0,
                   literal_source_coordinates=np.array_equal(old.vertices, emitted.vertices),
                   literal_source_attributes=all(np.array_equal(v, emitted.attributes[k]) for k, v in old.attributes.items()),
                   numerical=result)
        atomic_json(output / f'{side}-check.json', row)
        rows.append(row)
        if not result['passed']:
            break  # One failed side rejects the bilateral proposal immediately.
    result = dict(schema='tongue-genioglossus-source-diagonal/v1',
                  source_catalogue_sha256=source['sha256'], rows=rows,
                  passed=len(rows) == 2 and all(r['numerical']['passed'] for r in rows), installed=False,
                  evaluated_sides=[r['side'] for r in rows],
                  anatomical_acceptance=False,
                  scope='Only the two witnessed self-crossings; bony attachments still require separate correction',
                  elapsed_seconds=time.monotonic() - started)
    atomic_json(output / 'result.json', result)
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)
