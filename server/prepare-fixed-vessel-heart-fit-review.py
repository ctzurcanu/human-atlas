"""Review a whole-heart similarity with every surrounding vessel held fixed."""
from pathlib import Path
import argparse
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.review import prepare


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--proposal', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    proposal_path = Path(args.proposal).resolve()
    proposal = json.loads(proposal_path.read_text())
    source = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    sha = digest_file(source)
    assert proposal['source_catalogue_sha256'] == sha
    matrix = np.asarray(proposal['matrix'], float)
    scale = float(np.cbrt(np.linalg.det(matrix[:3, :3])))
    rotation = matrix[:3, :3] / scale
    assert scale > 0 and np.max(abs(rotation.T @ rotation - np.eye(3))) < 1e-10
    assert np.linalg.det(rotation) > 0
    out = Path(args.output).resolve()
    out.mkdir(exist_ok=False)
    closure_path = ROOT / '.local-models/geometry-workbench/heart/male-fixed-vessel-heart-closure-v1/closure.json'
    closure = json.loads(closure_path.read_text())
    assert closure['source_catalogue_sha256'] == sha
    moving = {row['part'] for row in closure['physical_heart_tissue']}
    assert len(moving) == 33 and not moving.intersection(closure['fixed_vascular_cuffs'])
    atlas = json.loads(source.read_text())
    parts = {part['id']: part for part in atlas['parts']}
    reader = AtlasReader(ROOT)
    frame = 'model:local-ta98:assembled-metres'
    entries = []
    for ordinal, pid in enumerate(sorted(moving)):
        part = parts[pid]
        mesh = reader.read_part(atlas, part, frame=frame, catalogue_path=source, catalogue_sha=sha)
        # Keep source face order so explicit tissue draw groups remain valid.
        used, inverse = np.unique(mesh.faces, return_inverse=True)
        attrs = {name: values[used] @ rotation.T if name == 'normals' else values[used]
                 for name, values in mesh.attributes.items()}
        placed = Mesh(mesh.vertices[used] @ matrix[:3, :3].T + matrix[:3, 3],
                      inverse.reshape(-1, 3), frame, attrs, mesh.source)
        path = out / f'heart-{ordinal}.npz'
        placed.save(path)
        entries.append(dict(mesh=str(path), name=part['name'], part_id=pid,
                            color=[.67, .34, .31, 1.],
                            tissue_draw_groups=part.get('tissueDrawGroups')))
    baseline = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-heart-placement-source-v1'
    surface_spec = json.loads((baseline / 'review-surface-vessels-spec.json').read_text())
    fixed_vessels = surface_spec['meshes'][1:]
    views = [dict(name='front', direction=[0, 0, 1]),
             dict(name='right', direction=[-1, 0, 0]),
             dict(name='back', direction=[0, 0, -1]),
             dict(name='oblique', direction=[.7, .35, 1])]
    fit_spec = dict(frame=frame, assembled_tissues=True, meshes=entries + fixed_vessels,
                    views=views, resolution=1100, wall_seconds=45)
    atomic_json(out / 'surface-spec.json', fit_spec)
    review = prepare(None, fit_spec)
    context_spec = json.loads((baseline / 'review-fixed-thorax-spec.json').read_text())
    fit_context = dict(fit_spec, meshes=entries + context_spec['meshes'][1:], views=views[:2])
    atomic_json(out / 'context-spec.json', fit_context)
    context_review = prepare(None, fit_context)
    assert digest_file(source) == sha
    atomic_json(out / 'reviews.json', dict(source_catalogue_sha256=sha,
                proposal_sha256=digest_file(proposal_path), matrix=matrix.tolist(),
                surface=review, context=context_review, heart_parts=33,
                all_vessels_fixed=True, source_banks_unchanged=True))
    print(json.dumps(dict(surface=review['request_id'], context=context_review['request_id'],
                          output=str(out), scale=scale)))


if __name__ == '__main__':
    main()
