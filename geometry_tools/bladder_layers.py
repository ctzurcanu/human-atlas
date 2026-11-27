"""First male bladder finite-coat candidate; review only, never installation."""
from pathlib import Path
import json

import numpy as np

from .algorithms.annular_pipe_envelopes import boundary_loops
from .audit import audit
from .common import ROOT, atomic_json, atomic_npz, digest_file
from .db import Database
from .fit import transform_mesh
from .harmonic_wall import layers
from .model import Mesh
from .registry import load_registered
from .review import prepare as prepare_review
from .tissue import pack


def build(source, source_review, output, *, fitted_wall=None):
    source, source_review, output = Path(source).resolve(), Path(source_review).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    evidence = json.loads((source / 'carriers.json').read_text())
    if not evidence['source_boundary_chain_exact'] or evidence['mouth_count'] != 3:
        raise ValueError('Verified three-mouth source carriers are required')
    db = Database(); _, urethra = load_registered(db, 'ZA:Urethra', 'local-ta98')
    uv, inverse = np.unique(urethra.vertices, axis=0, return_inverse=True)
    mouths = boundary_loops(inverse[urethra.faces])
    target_mouth = max(mouths, key=lambda q: uv[q].mean(0)[1])
    target_anchor = uv[target_mouth].mean(0)
    current = Mesh.load(source_review / 'male/current-bladder.npz')
    with np.load(source / 'carriers.npz', allow_pickle=False) as z:
        v, f = z['vertices'].copy(), z['wall_faces'].copy()
        inner, exterior = z['inner_side_faces'].copy(), z['outer_and_rim_faces'].copy()
        source_anchor = v[z['urethral_neck_ring']].mean(0)
    scale = float(np.median(np.ptp(current.bounds, axis=0) / np.ptp(v, axis=0)))
    matrix = np.eye(4); matrix[:3, :3] *= scale; matrix[:3, 3] = target_anchor - scale * source_anchor
    frame = current.frame
    q = transform_mesh(Mesh(v, f, evidence['frame']), matrix, frame).vertices.copy()
    # A common, injective dome field makes a derived construction rather than
    # an unchanged donor copy. Neck/ostial contours and all lower tissue stay fixed.
    height = float(np.ptp(q[:, 1])); start = float(q[:, 1].max() - height * .24)
    u = np.clip((q[:, 1] - start) / (height * .24), 0, 1); w = u * u * (3 - 2 * u)
    radial_center = (q.min(0) + q.max(0)) * .5
    displacement = np.zeros_like(q); displacement[:, [0, 2]] = (q[:, [0, 2]] - radial_center[[0, 2]]) * (.03 * w[:, None])
    q += displacement
    fitted_evidence = None
    if fitted_wall is not None:
        fitted_wall = Path(fitted_wall).resolve()
        corrected = Mesh.load(fitted_wall)
        fitted_evidence = json.loads((fitted_wall.parent / 'fit.json').read_text())
        if corrected.frame != frame or not np.array_equal(corrected.faces, f):
            raise ValueError('The fitted wall must preserve the registered frame and every source face')
        if fitted_evidence.get('passed') is not True:
            raise ValueError('Rejected pelvic fields cannot generate new coats')
        q = corrected.vertices.copy()
    carrier = Mesh(q, f, frame)
    check = audit(carrier, require_closed=True, expected_components=1, expected_genus=[2], max_pairs=1000000)
    output.mkdir(parents=True); carrier.save(output / 'whole-wall.npz')
    fit = dict(source_geometry_sha256=evidence['source_geometry_sha256'], target_geometry_sha256=current.geometry_sha,
               target_urethra_geometry_sha256=urethra.geometry_sha, uniform_scale=scale, matrix=matrix.tolist(),
               target_anchor_m=target_anchor.tolist(), source_anchor_m=source_anchor.tolist(),
               dome_field='x/z expansion up to 3 percent, C1 cubic blend over superior 24 percent; y fixed',
               maximum_dome_displacement_m=float(np.linalg.norm(displacement, axis=1).max()),
               continuous_field_positive_jacobian=True, whole_source_check=check,
               registration_scope='Gross bounds and urethral center only; port profiles and pelvic packing remain unfinished')
    if fitted_evidence is not None:
        fit.update(pelvic_field_sha256=digest_file(fitted_wall.parent / 'fit.json'),
                   fitted_wall_sha256=digest_file(fitted_wall), pelvic_field=fitted_evidence,
                   registration_scope='Bounded common neck field against unchanged pelvic neighbours; joins and final packing remain unfinished')
    atomic_json(output / 'fit.json', fit)
    if not check['passed']:
        raise ValueError('The common source fit introduces a geometry defect')
    names = ['external-connective', 'detrusor', 'submucosa', 'mucosa']
    fractions = [.08, .78, .94]
    inner_lips = np.intersect1d(np.unique(exterior), np.unique(inner))
    arrangement = layers(q, f, exterior, inner, names, method='boundary-harmonic',
                         rim_inner_vertices=inner_lips, interface_fractions=fractions)
    packed = pack(arrangement['positions'], arrangement['surfaces'], arrangement['cells'],
                  {name: 'BLADDER:male:' + name for name in names}, frame, output / 'tissues',
                  origin=np.round(q.mean(0) * 1024) / 1024, max_pairs=1000000,
                  construction=arrangement['evidence'])
    atomic_npz(output / 'source-carriers.npz', vertices=q, wall_faces=f, inner_faces=inner,
               exterior_and_rim_faces=exterior, frame=np.array(frame))
    limits = ['Review candidate; not installed or anatomically accepted.',
              'Relative coat fractions are explicit construction inputs, not histometric measurements.',
              'External connective coat still requires serosal/subserosal/adventitial regional allocation.',
              'Detrusor/trigonal/neck subdivisions, crest/uvula, attachments and mucosal folds remain unfinished.',
              'Ureteric and urethral profile joins and whole pelvic packing are unverified.']
    result = dict(schema='male-bladder-coats/v1', passed=packed['passed'], fit=fit,
                  source_carriers_sha256=digest_file(source / 'carriers.npz'),
                  tissue_manifest=str(output / 'tissues/tissue.json'), total_geometry_bytes=packed['total_bytes'],
                  positions_encoded_once=packed['positions_encoded_once'], new_images=0, new_uvs=0,
                  layers=names, interface_fractions=fractions, remaining=limits,
                  anatomical_acceptance=False, installed=False)
    if packed['passed']:
        colors = [[.78, .58, .42, 1], [.67, .24, .29, 1], [.79, .64, .42, 1], [.88, .48, .52, 1]]
        entries = [dict(mesh=str(output / 'tissues' / ('BLADDER_male_' + name + '.npz')), name=name, color=color)
                   for name, color in zip(names, colors)]
        center = carrier.bounds.mean(0)
        request = prepare_review(db, dict(frame=frame, meshes=entries,
            views=[dict(name='whole-wall', direction=[.5, .2, 1]),
                   dict(name='coronal-four-coat-section', direction=[.1, .2, 1],
                        section=dict(normal=[0, 0, 1], point=center.tolist(), clear_outer=True)),
                   dict(name='sagittal-four-coat-section', direction=[1, .2, .1],
                        section=dict(normal=[1, 0, 0], point=center.tolist(), clear_outer=True))], resolution=900))
        result['review'] = request
    atomic_json(output / 'candidate.json', result)
    return result


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--source', required=True)
    p.add_argument('--source-review', required=True); p.add_argument('--output', required=True)
    p.add_argument('--fitted-wall')
    args = p.parse_args(); result = build(args.source, args.source_review, args.output, fitted_wall=args.fitted_wall)
    print(json.dumps({k: result[k] for k in ['passed', 'layers', 'total_geometry_bytes', 'installed']}, indent=2))
