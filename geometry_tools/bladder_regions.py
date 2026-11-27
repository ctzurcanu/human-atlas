"""Finite bladder coats and muscle regions in one source-bound arrangement.

The three source ostia define the trigonal construction frame. Regional limits
and coat fractions remain explicit review inputs, not measured histometry.
Only existing wall tetrahedra are partitioned; no opening receives a solid disk.
"""
from pathlib import Path
import itertools
import json

import numpy as np

from .common import atomic_json, atomic_npz, digest_file
from .db import Database
from .harmonic_wall import graph_values, tetrahedralize
from .implicit_cells import Arrangement
from .model import Mesh
from .review import prepare
from .tissue import pack


def construct(wall_path, carrier_path, output):
    wall_path, carrier_path, output = Path(wall_path).resolve(), Path(carrier_path).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    fit = json.loads((wall_path.parent / 'fit.json').read_text())
    if fit.get('passed') is not True:
        raise ValueError('Use a passed common-wall candidate')
    wall = Mesh.load(wall_path); v, f = wall.vertices, wall.faces
    with np.load(carrier_path, allow_pickle=False) as z:
        if not np.array_equal(z['wall_faces'], f):
            raise ValueError('The whole source face sequence changed')
        inner, exterior = z['inner_side_faces'], z['outer_and_rim_faces']
        ports = {name: z[name + '_ring'].copy() for name in ['urethral_neck', 'left_ureteric', 'right_ureteric']}
    p, tets, _, _, mesher = tetrahedralize(v, f)
    lookup = {tuple(q): i for i, q in enumerate(p)}
    # The remainder contains the physical mouth annuli as well as exterior
    # skin. Their common controls are INNER lips, not pinched zero-thickness
    # outer/inner seams. Fix those controls to mucosa so it reaches every ostium.
    inner_controls = set(map(int, np.unique(inner)))
    outer = [lookup[tuple(v[i])] for i in np.unique(exterior) if int(i) not in inner_controls]
    lining = [lookup[tuple(v[i])] for i in inner_controls]
    boundary = np.array([[lookup[tuple(v[i])] for i in face] for face in f], int)
    edges = np.concatenate([boundary[:, [i, j]] for i, j in [(0, 1), (1, 2), (2, 0)]])
    surface, surface_evidence = graph_values(p, edges, np.r_[outer, lining], np.r_[np.zeros(len(outer)), np.ones(len(lining))])
    edges = np.concatenate([tets[:, [i, j]] for i, j in itertools.combinations(range(4), 2)])
    fixed = np.unique(boundary)
    phi, harmonic = graph_values(p, edges, fixed, surface[fixed])

    neck, left, right = [v[ports[n]].mean(0) for n in ['urethral_neck', 'left_ureteric', 'right_ureteric']]
    basis = np.column_stack([left - neck, right - neck])
    uv = (p - neck) @ np.linalg.pinv(basis).T
    bary = np.c_[1. - uv.sum(1), uv]
    normal = np.cross(left - neck, right - neck); normal /= np.linalg.norm(normal)
    if normal[2] < 0:
        normal *= -1
    distance = (p - neck) @ normal
    neck_height = float(neck[1] + .008)
    cover_height = float((left[1] + right[1]) / 2 + .004)
    fundus_height = float(max(left[1], right[1]) + .004)
    fractions = [.015, .08, .30, .45, .55, .60, .78, .94]
    columns = [np.full(len(p), -1.), phi, *[phi - c for c in fractions],
               p[:, 1] - neck_height, p[:, 1] - cover_height,
               .010 - distance, fundus_height - p[:, 1],
               *[bary[:, i] + .08 for i in range(3)], .008 - distance, distance + .012]
    fields = np.column_stack(columns)
    neck_index = 2 + len(fractions)
    cover_index, fundus_distance, fundus_top = neck_index + 1, neck_index + 2, neck_index + 3
    trig_start = neck_index + 4
    def label(x):
        value = x[1]; neck_region = x[neck_index] <= 0
        trigone = bool(np.all(x[trig_start:trig_start + 5] >= 0))
        region = 'neck' if neck_region else 'fundus' if trigone or x[fundus_distance] >= 0 and x[fundus_top] >= 0 else 'body'
        if value < .08:
            tissue = ('serosa' if value < .015 else 'subserosa') if x[cover_index] >= 0 else 'adventitia'
        elif value < .78:
            if neck_region:
                tissue = 'neck-external-longitudinal' if value < .30 else 'neck-circular' if value < .55 else 'neck-internal-longitudinal'
            elif trigone and value >= .45:
                tissue = 'trigone-deep' if value < .60 else 'trigone-superficial'
            else:
                tissue = 'detrusor-unstratified'
        elif value < .94:
            tissue = 'submucosa'
        else:
            tissue = 'mucosa-trigone' if trigone else 'mucosa'
        return region + ':' + tissue
    origin = np.round(p.mean(0) * 1024) / 1024
    arrangement = Arrangement(p - origin, fields, label)
    arrangement.facet_labels = {tuple(sorted(face)): i for i, face in enumerate(boundary)}
    positions, surfaces, cells = arrangement.build(tets, skip_fields=[0, 1]); positions += origin
    inner_keys = {tuple(sorted(face)) for face in inner}
    inner_ids = {i for i, face in enumerate(f) if tuple(sorted(face)) in inner_keys}
    covered = set(); wrong_lining = []
    for name, uses in cells.items():
        for use in uses:
            tags = arrangement.surface_source_faces[use['surface']]
            relevant = set(map(int, tags[tags >= 0])) & inner_ids
            if relevant:
                covered |= relevant
                if not name.endswith((':mucosa', ':mucosa-trigone')):
                    wrong_lining.append(name)
    if covered != inner_ids or wrong_lining:
        raise ValueError('Mucosa must own every original inner-wall face, including every ostial lip')
    output.mkdir(parents=True)
    result = pack(positions, surfaces, cells, {name: 'BLADDER:male:' + name for name in cells},
                  wall.frame, output / 'tissues', origin=origin, max_pairs=1000000,
                  construction=dict(**mesher, harmonic=harmonic, boundary_harmonic=surface_evidence,
                                    source_wall_boundary_preserved=True, clinical_thickness_measured=False))
    plan = dict(schema='bladder-landmark-regions/v1', wall=str(wall_path), wall_sha256=digest_file(wall_path),
                carriers_sha256=digest_file(carrier_path), port_centers_m=dict(zip(['neck', 'left', 'right'], [neck.tolist(), left.tolist(), right.tolist()])),
                wall_fractions=fractions, neck_cut_height_m=neck_height, serosal_cut_height_m=cover_height,
                fundus_top_m=fundus_height, trigonal_plane_normal=normal.tolist(), trigonal_barycentric_margin=.08,
                regions=list(cells), geometry_bytes=result['total_bytes'], passed=result['passed'],
                inner_wall_source_faces=len(inner_ids), all_inner_wall_owned_by_mucosa=True,
                every_inner_ostial_lip_is_mucosal=True,
                source_faces_partitioned_once=True, lumen_disks_emitted=False, source_fluid_ports_retained=True,
                anatomical_boundaries_provisional=True, installed=False, anatomical_acceptance=False,
                remaining=['Ostial and urethral shared profiles; pelvic tissue interfaces.',
                           'Interureteric crest, uvula, median umbilical ligament and muscular attachments.',
                           'Full histological and anatomical surface/section review; local fit and source rights.'])
    # Selection sets reference the existing finite tissue cells. They never
    # add a coincident whole-organ copy or turn an opening into a filler mesh.
    groups = {'A08.3.01.001': list(cells)}
    for macro, code in [('body', '003'), ('fundus', '004'), ('neck', '005')]:
        groups['A08.3.01.' + code] = [name for name in cells if name.startswith(macro + ':')]
    for tissue, code in [('serosa', '008'), ('subserosa', '009'), ('detrusor-unstratified', '015'),
                         ('trigone-superficial', '012'), ('trigone-deep', '013'),
                         ('neck-external-longitudinal', '017'), ('neck-circular', '018'),
                         ('neck-internal-longitudinal', '019'), ('submucosa', '022')]:
        groups['A08.3.01.' + code] = [name for name in cells if name.split(':', 1)[1] == tissue]
    groups['A08.3.01.023'] = [name for name in cells if name.split(':', 1)[1] in ['mucosa', 'mucosa-trigone']]
    groups['A08.3.01.010'] = [name for name in cells if any(w in name for w in ['detrusor-', 'trigone-deep', 'trigone-superficial', 'neck-external-', 'neck-circular', 'neck-internal-'])]
    groups['A08.3.01.011'] = [name for name in cells if name.endswith(('trigone-deep', 'trigone-superficial'))]
    groups['A08.3.01.016'] = [name for name in groups['A08.3.01.010'] if name.startswith('neck:')]
    groups['A08.3.01.014'] = [name for name in groups['A08.3.01.010'] if name not in groups['A08.3.01.011']]
    groups['A08.3.01.024'] = [name for name in cells if name.endswith(('mucosa-trigone', 'trigone-deep', 'trigone-superficial'))]
    if any(not members for members in groups.values()):
        raise ValueError('A planned TA98 tissue/region has no actual finite representation')
    plan['ta98_selection_groups'] = {term: [result['part_ids'][name] for name in names] for term, names in groups.items()}
    plan['physiological_boundaries'] = [dict(ta98='A08.3.01.027M' if name == 'urethral_neck' else 'A08.3.01.026',
            side=None if name == 'urethral_neck' else name.split('_')[0],
            source_ring_controls=ids.tolist(), points_m=v[ids].tolist(),
            representation='actual open luminal boundary; not a separate solid', closed_ring=True)
        for name, ids in ports.items()]
    if result['passed']:
        entries = []
        for name, pid in result['part_ids'].items():
            color = ([.82, .57, .49, 1] if 'serosa' in name else [.65, .58, .44, 1] if 'adventitia' in name else
                     [.89, .52, .55, 1] if 'mucosa' in name and 'submucosa' not in name else [.84, .69, .45, 1] if 'submucosa' in name else
                     [.63, .22, .28, 1] if 'circular' in name else [.72, .35, .37, 1])
            entries.append(dict(mesh=str(output / 'tissues' / (pid.replace(':', '_') + '.npz')), name=name, color=color))
        center = wall.bounds.mean(0)
        plan['review'] = prepare(Database(), dict(frame=wall.frame, meshes=entries, resolution=1000,
            views=[dict(name='whole-regional-wall', direction=[.3, .2, 1]),
                   dict(name='posterior-trigone-section', direction=[0, .2, 1],
                        section=dict(normal=[0, 0, 1], point=center.tolist(), clear_outer=True)),
                   dict(name='neck-muscle-section', direction=[1, .1, 0],
                        section=dict(normal=[1, 0, 0], point=neck.tolist(), clear_outer=True))]))
    atomic_json(output / 'regions.json', plan)
    return plan


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('wall'); p.add_argument('--carriers', required=True); p.add_argument('--output', required=True)
    args = p.parse_args(); r = construct(args.wall, args.carriers, args.output)
    print(json.dumps({k: r[k] for k in ['passed', 'regions', 'geometry_bytes']}, indent=2))
