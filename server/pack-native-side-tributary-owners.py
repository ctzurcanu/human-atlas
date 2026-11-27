"""Two disjoint owners of a completed, source-bound three-port finite wall.

Prepare and pack are separate supervised stages; neither installs anatomy.
Spec: schema native-side-tributary-two-owner-pack/v1, whole {directory,
assemble_sha256}, part_ids/terms {collector, tributary}, cut_branch:'tributary',
origin:[x,y,z]. Optional normal_charts[parent|tributary] supplies normal_array,
uv_array and encoding raw-int16 or transported-float64. The latter requires a
hash-bound transport_proof with proper_normal_transport:true,
normal_sign_flips:0 and target_source_sha256 matching the actual source NPZ.
Transported charts may also supply literal_normal_words_array and
literal_normal_mask_array; the same target-bound proof must name both arrays.
"""
from pathlib import Path
from collections import Counter
import argparse
import importlib.util
import json
import re
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file  # noqa: E402
from geometry_tools.model import Mesh  # noqa: E402
from geometry_tools.algorithms.finite_tissue_boolean import oriented_key  # noqa: E402
from geometry_tools.algorithms.tube_continuity import annular_rim, boundary_loops  # noqa: E402
from geometry_tools.tissue import normal_words, pack  # noqa: E402

HELPER = ROOT/'server/pack-superior-gluteal-arterial-owners.py'
loader = importlib.util.spec_from_file_location('native_owner_topology_helpers', HELPER)
helpers = importlib.util.module_from_spec(loader)
loader.loader.exec_module(helpers)
NAMES = ('collector', 'tributary')
GENUS = dict(collector=2, tributary=1)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def path(value, base):
    p = Path(value)
    return (p if p.is_absolute() else base/p).resolve()


def bound(binding, base, files):
    p = path(binding['path'], base)
    require(digest_file(p) == binding['sha256'], f'Changed bound input: {p}')
    files[str(p)] = binding['sha256']
    return p


def load_spec(args):
    p = args.spec.resolve()
    require(digest_file(p) == args.spec_sha256, 'Owner spec SHA256 mismatch')
    s = json.loads(p.read_text())
    require(s['schema'] == 'native-side-tributary-two-owner-pack/v1', 'Unknown owner spec')
    require(s['cut_branch'] == 'tributary', 'Only the explicit single tributary cut is supported')
    require(set(s['part_ids']) == set(s['terms']) == set(NAMES), 'Two explicit IDs and TA98 terms required')
    require(len(set(s['part_ids'].values())) == 2 and all(s['part_ids'].values()), 'Distinct nonempty owner IDs required')
    require(all(re.fullmatch(r'A\d{2}\.\d\.\d{2}\.\d{3}[MF]?', t) for t in s['terms'].values()),
            'Explicit canonical TA98 codes required')
    origin = np.asarray(s['origin'], float)
    require(origin.shape == (3,) and np.isfinite(origin).all(), 'Explicit finite packing origin required')
    return p, s


def source_charts(plan, plan_path, spec, files):
    result = []
    for role in ('parent', 'tributary'):
        p = bound(plan[role], plan_path.parent, files)
        with np.load(p, allow_pickle=False) as z:
            chart = {k: z[k].copy() for k in ('vertices', 'faces', 'profiles', 'profile_vertex_indices')}
            config = spec.get('normal_charts', {}).get(role, {})
            chart['normals'] = z[config.get('normal_array', 'attribute_normals')].copy()
            chart['uvs'] = z[config.get('uv_array', 'attribute_uvs')].copy()
            literal_keys = [config.get('literal_normal_words_array'), config.get('literal_normal_mask_array')]
            require(all(literal_keys) or not any(literal_keys), 'Literal normal words and mask must be declared together')
            if all(literal_keys):
                chart['literal_words'] = z[literal_keys[0]].copy()
                chart['literal_mask'] = z[literal_keys[1]].copy()
        encoding = config.get('encoding', 'raw-int16')
        normal = chart['normals']
        require(normal.shape == chart['vertices'].shape and np.isfinite(normal).all(), 'Native normal chart shape/value failure')
        require(chart['uvs'].shape == (len(normal), 2) and np.isfinite(chart['uvs']).all(), 'Native UV chart missing/nonfinite')
        if encoding == 'raw-int16':
            require(normal.dtype == np.dtype('int16'), 'Raw normal claims require actual Int16 source words')
            require(not any(literal_keys), 'A raw Int16 chart already supplies its literal words')
        elif encoding == 'transported-float64':
            require(normal.dtype == np.dtype('float64'), 'Transported chart requires actual Float64 normals')
            proof_path = bound(config['transport_proof'], plan_path.parent, files)
            proof = json.loads(proof_path.read_text())
            require(proof.get('proper_normal_transport') is True and proof.get('normal_sign_flips') == 0
                    and proof.get('target_source_sha256') == digest_file(p),
                    'Float chart lacks a properly transported, unflipped target-bound proof')
            if all(literal_keys):
                require(proof.get('literal_normal_words_array') == literal_keys[0]
                        and proof.get('literal_normal_mask_array') == literal_keys[1],
                        'Target-bound proof does not name the actual literal normal arrays')
                literal, mask = chart['literal_words'], chart['literal_mask']
                require(literal.dtype == np.dtype('int16') and literal.shape == normal.shape,
                        'Literal normal words require actual Int16[N,3] source data')
                require(mask.dtype == np.bool_ and mask.shape == (len(normal),),
                        'Literal normal mask requires actual bool[N] source data')
                if mask.any():
                    native_length = np.linalg.norm(literal[mask].astype(float), axis=1)
                    target_length = np.linalg.norm(normal[mask], axis=1)
                    require(np.all(native_length > 0) and np.all(target_length > 0), 'Zero masked literal normal')
                    native_direction = literal[mask].astype(float)/native_length[:, None]
                    target_direction = normal[mask]/target_length[:, None]
                    require(np.linalg.norm(native_direction-target_direction, axis=1).max() <= 1e-12,
                            'Masked Float64 normal does not preserve the literal native word direction')
        else:
            raise ValueError(f'Unknown source chart encoding: {encoding}')
        require(np.linalg.norm(normal, axis=1).min() > 0, 'Zero source normal')
        chart.update(encoding=encoding, path=str(p), sha256=digest_file(p))
        result.append(chart)
    return result


def quantized(value):
    length = float(np.linalg.norm(value))
    require(np.isfinite(length) and length > 0, 'Zero/nonfinite interpolated chart normal')
    rounded = np.rint(np.asarray(value, float)/length*32767)
    require(np.isfinite(rounded).all() and rounded.min() >= -32768 and rounded.max() <= 32767,
            'Int16 normal-word range overflow')
    return rounded.astype('<i2')


def charts(positions, encoded, surfaces, masks, data, sources):
    words = normal_words(encoded, surfaces)
    uvs = {sid: np.zeros((len(f), 3, 2), dtype='<f4') for sid, f in surfaces.items()}
    replay = {}
    counts = {role: dict(raw_literal_words=0, raw_barycentric_quantized=0,
                          transported_chart_literal_raw_words=0,
                          float64_literal_chart_quantized_once=0, float64_barycentric_quantized_once=0)
              for role in ('parent', 'tributary')}
    for owner in NAMES:
        sid = owner+':outer'
        components = data['outer_component'][masks[owner]]
        source_faces = data['outer_source_face'][masks[owner]]
        bary = np.empty((len(components), 3, 3))
        for fi, (face, component, source_face) in enumerate(zip(surfaces[sid], components, source_faces)):
            c, sf = int(component), int(source_face)
            require(c in (0, 1), 'Outer facet has unknown native component ancestry')
            src = sources[c]
            require(0 <= sf < len(src['faces']),
                    'Outer facet lacks a native source normal/UV chart; generated extension charts must be supplied explicitly')
            source_indices = src['faces'][sf]
            triangle = src['vertices'][source_indices]
            edge = triangle[1:] - triangle[0]
            gram = edge @ edge.T
            source_face_normal = np.cross(positions[face[1]]-positions[face[0]],
                                          positions[face[2]]-positions[face[0]])
            face_normal = np.cross(encoded[face[1]]-encoded[face[0]], encoded[face[2]]-encoded[face[0]])
            require(np.linalg.norm(face_normal) > 0, 'Position encoding collapsed a surface facet')
            for corner, point in enumerate(positions[face]):
                same = np.flatnonzero(np.all(triangle == point, axis=1))
                if len(same):
                    weights = np.eye(3)[same[0]]
                    value = src['normals'][source_indices[int(same[0])]]
                else:
                    xy = np.linalg.solve(gram, edge @ (point-triangle[0]))
                    weights = np.r_[1-xy.sum(), xy]
                    require(weights.min() >= -1e-8 and weights.max() <= 1+1e-8
                            and np.linalg.norm(weights @ triangle-point) <= 1e-11,
                            'Clipped corner is outside its actual native source triangle chart')
                    value = weights @ src['normals'][source_indices].astype(float)
                raw_literal = bool(len(same) and src['encoding'] == 'raw-int16')
                transported_literal = bool(len(same) and src['encoding'] == 'transported-float64'
                    and 'literal_mask' in src and src['literal_mask'][source_indices[int(same[0])]])
                word = (value.copy() if raw_literal else
                        src['literal_words'][source_indices[int(same[0])]].copy() if transported_literal else
                        quantized(value))
                require(word.dtype == np.dtype('int16') and word.astype(float) @ face_normal > 0
                        and word.astype(float) @ source_face_normal > 0,
                        'Nonoutward source normal after position encoding; no sign flip permitted')
                words[sid][fi, corner] = word
                uvs[sid][fi, corner] = weights @ src['uvs'][source_indices]
                bary[fi, corner] = weights
                key = ('raw_literal_words' if raw_literal else
                       'transported_chart_literal_raw_words' if transported_literal else
                       'raw_barycentric_quantized' if src['encoding'] == 'raw-int16' else
                       'float64_literal_chart_quantized_once' if len(same) else
                       'float64_barycentric_quantized_once')
                counts['parent' if c == 0 else 'tributary'][key] += 1
        replay[owner+'_component'], replay[owner+'_source_face'], replay[owner+'_barycentric'] = components, source_faces, bary
    return words, uvs, replay, counts


def verify_matching_audit(row, mesh, label):
    """A pass for another indexed geometry must never validate this owner bank."""
    require(row.get('passed') is True, f'{label} audit failed')
    require(row.get('geometry_sha256') == mesh.geometry_sha,
            f'{label} audit belongs to a different actual geometry')
    require(row.get('frame') == mesh.frame and row.get('attributes_sha256') == mesh.attributes_sha,
            f'{label} audit frame or attributes do not match the actual mesh')


def prepare(spec_path, spec, output):
    began = time.monotonic()
    require(not output.exists(), 'Prepare output must be new')
    whole_dir = path(spec['whole']['directory'], spec_path.parent)
    files = {str(spec_path): digest_file(spec_path), str(HELPER): digest_file(HELPER)}
    def read(name, sha=None):
        p = whole_dir/name
        value = digest_file(p)
        require(sha is None or value == sha, f'Changed completed whole artifact: {name}')
        files[str(p)] = value
        return json.loads(p.read_text()) if p.suffix == '.json' else p
    receipt = read('assemble.json', spec['whole']['assemble_sha256'])
    require(receipt.get('passed') and not receipt.get('installed') and not receipt.get('anatomical_acceptance'),
            'Completed numerical, uninstalled whole-wall receipt required')
    plan_path = bound(dict(path=receipt['plan'], sha256=receipt['plan_sha256']), whole_dir, files)
    plan = json.loads(plan_path.read_text())
    require(plan['schema'] == 'native-side-tributary-finite-wall/v1', 'Wrong whole producer schema')
    for g in [plan['inventory'], *plan['guards']]:
        bound(g, plan_path.parent, files)
    whole_producer = ROOT/'server/construct-native-side-tributary-wall.py'
    require(digest_file(whole_producer) == receipt['producer_sha256'], 'Whole producer revision changed')
    files[str(whole_producer)] = receipt['producer_sha256']
    for name in ('prepare.json', 'outer.json', 'inner.json'):
        stage = read(name)
        require(stage['plan_sha256'] == receipt['plan_sha256']
                and stage['producer_sha256'] == receipt['producer_sha256'], 'Mixed construction stage revisions')
        if name in ('outer.json', 'inner.json'):
            level = name[:-5]
            for field, suffix in [('sheet_sha256', 'sheet.npz'), ('ancestry_sha256', 'ancestry.npz'),
                                  ('classifier_sha256', 'classifier.npz'),
                                  ('classifier_provenance_sha256', 'classifier-provenance.npz')]:
                read(level+'-'+suffix, stage[field])
        else:
            read('prepared.npz', stage['prepared_sha256'])
    containment = read('inner-containment.json', receipt['inner_containment_sha256'])
    require(containment['passed'] and containment['contained_away_from_external_classifier_caps']
            and containment['outside_actual_wall_source_faces'] == 0,
            'Whole classifier containment failed or actual wall leakage was ignored')
    if not containment['strict_empty']:
        require(containment['nonempty_difference_is_external_cap_roundoff_only']
                and containment['finite_roundoff_volume_bound_passed'], 'Nonempty classifier difference is not proven cap-only roundoff')
    bound(dict(path=containment['outside_artifact'], sha256=containment['outside_artifact_sha256']), whole_dir, files)
    whole = Mesh.load(read('finite-wall.npz', receipt['finite_wall_sha256']))
    require(whole.geometry_sha == receipt['geometry_sha256'] and whole.frame == plan['frame'], 'Whole geometry/frame mismatch')
    verify_matching_audit(read('wall-audit.json'), whole, 'Whole wall')
    blood = Mesh.load(whole_dir/'inner-classifier.npz')
    verify_matching_audit(read('blood-audit.json'), blood, 'Blood classifier')
    with np.load(read('surfaces-and-ancestry.npz'), allow_pickle=False) as z:
        data = {k: z[k].copy() for k in z.files}
    port_rows = read('ports.json')['ports']
    require([(p['source'], p['mouth']) for p in port_rows] == [('parent', 0), ('parent', 1), ('tributary', 0)],
            'Expected exact parent/parent/tributary external port identities')
    positions, outer, inner = data['vertices'], data['outer_faces'], data['inner_faces']
    annuli = [data[f'annulus_{i}'] for i in range(3)]
    full = np.concatenate([outer, inner, *annuli])
    require(np.array_equal(positions, whole.vertices) and np.array_equal(full, whole.faces),
            'Surface ancestry arrangement does not replay actual completed whole wall')
    edges, whole_topology = helpers.closed_edge_check(full)
    require(whole_topology['genus'] == 2, 'Whole finite wall genus must be two')
    lookup = {tuple(positions[i]): int(i) for i in np.unique(full)}
    require(len(lookup) == len(np.unique(full)), 'Shared whole controls are not canonical')
    sources = source_charts(plan, plan_path, spec, files)
    profiles = sources[1]['profiles']
    centers = profiles.mean(1)
    alpha = float(plan['inner_radial_fraction'])
    upstream = int(plan['branch_mouth'])
    candidates = range(1, len(profiles)-1) if upstream == 0 else range(len(profiles)-2, 0, -1)
    cut = None
    for station in candidates:
        pair = np.array([profiles[station], centers[station]+alpha*(profiles[station]-centers[station])])
        try:
            ids = np.array([[lookup[tuple(p)] for p in ring] for ring in pair])
        except KeyError:
            continue
        blocked = {tuple(sorted((int(a), int(b)))) for ring in ids
                   for a, b in zip(ring, np.roll(ring, -1))}
        if len(blocked) == 2*profiles.shape[1] and all(e in edges for e in blocked):
            cut = dict(source_station=int(station), upstream_mouth=upstream,
                       outer=ids[0], inner=ids[1], blocked=blocked)
            break
    require(cut is not None, 'No intact paired native tributary ring beyond the local Boolean')
    rows, cols = [], []
    for e, uses in edges.items():
        if e not in cut['blocked']:
            a, b = uses[0][0], uses[1][0]
            rows.extend([a, b]); cols.extend([b, a])
    count, labels = connected_components(coo_matrix((np.ones(len(rows)), (rows, cols)),
                                         shape=(len(full), len(full))).tocsr(), directed=False)
    require(count == 2, 'The paired cut must make exactly two externally rimmed tissue regions')
    offsets = np.cumsum([0, len(outer), len(inner), *map(len, annuli)])
    owner_by_label = {}
    for i, owner in enumerate(('collector', 'collector', 'tributary')):
        unique = np.unique(labels[offsets[2+i]:offsets[3+i]])
        require(len(unique) == 1, 'An external rim crosses an ownership region')
        label = int(unique[0])
        require(label not in owner_by_label or owner_by_label[label] == owner, 'Distinct external blood mouths share the wrong tissue region')
        owner_by_label[label] = owner
    require(set(owner_by_label.values()) == set(NAMES), 'Both finite owners must retain their external mouth')
    face_owner = np.array([owner_by_label[int(label)] for label in labels])
    surfaces, cells, masks = {}, {name: [] for name in NAMES}, {}
    for owner in NAMES:
        masks[owner] = face_owner[:len(outer)] == owner
        for sheet_name, faces, start in [('outer', outer, 0), ('inner', inner, len(outer))]:
            mask = face_owner[start:start+len(faces)] == owner
            require(mask.any(), 'Finite owner lacks an inner or outer sheet')
            if owner == 'tributary':
                require(set(map(int, data[sheet_name+'_component'][mask])) == {1}, 'Distal tributary contains collector source faces')
            sid = owner+':'+sheet_name
            surfaces[sid] = faces[mask]
            cells[owner].append(dict(surface=sid, side=1, neighbor=None))
    for i, owner in enumerate(('collector', 'collector', 'tributary')):
        sid = f'external-mouth:{i}'
        surfaces[sid] = annuli[i]
        cells[owner].append(dict(surface=sid, side=1, neighbor=None))
    outer_loop = next((q for q in boundary_loops(surfaces['tributary:outer']) if set(q) == set(cut['outer'])), None)
    inner_loop = next((q for q in boundary_loops(surfaces['tributary:inner']) if set(q) == set(cut['inner'])), None)
    require(outer_loop is not None and inner_loop is not None, 'Paired cut is not the actual distal owner boundary')
    permutation = [int(np.flatnonzero(cut['outer'] == i)[0]) for i in outer_loop]
    surfaces['shared-wall-ring'] = annular_rim(outer_loop, cut['inner'][permutation])
    cells['tributary'].append(dict(surface='shared-wall-ring', side=1, neighbor='collector'))
    cells['collector'].append(dict(surface='shared-wall-ring', side=-1, neighbor='tributary'))
    combined = Counter()
    topologies = {}
    for owner in NAMES:
        faces = np.concatenate([surfaces[u['surface']] if u['side'] == 1 else surfaces[u['surface']][:, [0, 2, 1]]
                                for u in cells[owner]])
        _, topologies[owner] = helpers.closed_edge_check(faces)
        require(topologies[owner]['genus'] == GENUS[owner], 'Unexpected finite owner genus')
        combined.update(helpers.oriented_counter(faces))
    cancelled = Counter({f: n-combined[oriented_key([f[0], f[2], f[1]])] for f, n in combined.items()
                         if n > combined[oriented_key([f[0], f[2], f[1]])]})
    require(cancelled == helpers.oriented_counter(full), 'Disjoint ownership and reciprocal interface fail exact whole replay')
    origin = np.asarray(spec['origin'], float)
    encoded = (positions-origin).astype('<f4').astype(float)+origin
    words, uvs, chart_replay, counts = charts(positions, encoded, surfaces, masks, data, sources)
    arrays, sm, wm, um = dict(positions=positions), {}, {}, {}
    for i, (sid, faces) in enumerate(surfaces.items()):
        sm[sid], wm[sid], um[sid] = f'faces_{i}', f'words_{i}', f'uvs_{i}'
        arrays[sm[sid]], arrays[wm[sid]], arrays[um[sid]] = faces, words[sid], uvs[sid]
    output.mkdir(parents=True)
    arrangement = output/'prepared-arrangement.npz'
    atomic_npz(arrangement, **arrays)
    atomic_npz(output/'outer-chart-ancestry.npz', **chart_replay)
    pack_spec = dict(source=str(arrangement), frame=whole.frame, origin=origin.tolist(), surfaces=sm,
                     corner_words=wm, corner_uvs=um, cells=cells, part_ids=spec['part_ids'],
                     max_pairs=250000, max_tests=2000000)
    atomic_json(output/'pack-spec.json', pack_spec)
    receipt = dict(schema='native-side-tributary-two-owner-prepared/v1', preparation_passed=True,
        whole_geometry_sha256=whole.geometry_sha, terms=spec['terms'], part_ids=spec['part_ids'],
        source_directory=str(whole_dir), cells=topologies, whole_wall_cancellation_exact=True,
        all_original_whole_faces_owned_once=True, reciprocal_annular_interface=True, lumen_disks=0,
        cut=dict(source_station=cut['source_station'], upstream_mouth=upstream,
                 controls_per_ring=profiles.shape[1], first_intact_paired_source_ring=True),
        native_chart_counts=counts, normal_sign_flips=0,
        chart_claim='Raw Int16 words remain literal at raw-source corners and explicitly proved unchanged transported-chart corners. Other transported Float64 and barycentric normals are quantized once. New inner/rim charts are generated geometric fields.',
        source_files=files,
        prepared_files={str(p): digest_file(p) for p in
                        [arrangement, output/'outer-chart-ancestry.npz', output/'pack-spec.json']},
        producer_sha256=digest_file(Path(__file__)), anatomical_acceptance=False, installed=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'prepared.json', receipt)
    return receipt


def pack_prepared(spec_path, output):
    began = time.monotonic()
    prepared = json.loads((output/'prepared.json').read_text())
    require(prepared['preparation_passed'] and prepared['producer_sha256'] == digest_file(Path(__file__)), 'Prepared producer changed')
    require(str(spec_path) in prepared['source_files'], 'Prepared arrangement belongs to another spec')
    for p, sha in (prepared['source_files'] | prepared['prepared_files']).items():
        require(digest_file(Path(p)) == sha, f'Changed saved pack input: {p}')
    spec = json.loads((output/'pack-spec.json').read_text())
    with np.load(spec['source'], allow_pickle=False) as data:
        result = pack(data['positions'], {k: data[v] for k, v in spec['surfaces'].items()}, spec['cells'],
            spec['part_ids'], spec['frame'], output/'finite-cells', origin=spec['origin'],
            corner_words={k: data[v] for k, v in spec['corner_words'].items()},
            corner_uvs={k: data[v] for k, v in spec['corner_uvs'].items()},
            max_pairs=spec['max_pairs'], max_tests=spec['max_tests'], checkpoint_directory=output/'validation-checkpoints')
    require(result['passed'], 'Decoded finite owner checks failed')
    for owner, row in result['checks']['cells'].items():
        components = row['topology']['components']
        require(len(components) == 1 and components[0]['closed_orientable_genus'] == GENUS[owner], 'Decoded owner topology changed')
    evidence = dict(schema='native-side-tributary-two-owner-packed/v1', finite_owner_checks_passed=True,
        finite_cells=2, terms=prepared['terms'], part_ids=prepared['part_ids'], cut=prepared['cut'],
        whole_geometry_sha256=prepared['whole_geometry_sha256'], all_original_whole_faces_owned_once=True,
        reciprocal_annular_interface=True, whole_wall_cancellation_exact=True, lumen_disks=0,
        blood_classifier_unchanged=True, native_chart_counts=prepared['native_chart_counts'],
        normal_sign_flips=0, bank_sha256=result['bank_sha256'], bank_bytes=result['total_bytes'],
        anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'evidence.json', evidence)
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec', type=Path, required=True)
    parser.add_argument('--spec-sha256', required=True)
    parser.add_argument('--output-directory', type=Path, required=True)
    parser.add_argument('--stage', choices=('prepare', 'pack'), required=True)
    args = parser.parse_args()
    spec_path, spec = load_spec(args)
    output = args.output_directory.resolve()
    if args.stage == 'pack':
        require(not (output/'evidence.json').exists() and not (output/'finite-cells').exists(), 'No silent packing retry or overwrite')
    try:
        result = prepare(spec_path, spec, output) if args.stage == 'prepare' else pack_prepared(spec_path, output)
        print(json.dumps(result))
    except Exception as error:
        if output.exists() and not (output/(args.stage+'-failure.json')).exists():
            atomic_json(output/(args.stage+'-failure.json'), dict(error_type=type(error).__name__, error=str(error),
                stage=args.stage, installed=False, anatomical_acceptance=False))
        raise


if __name__ == '__main__':
    main()
