"""Read-only source-support audit; never creates or installs geometry."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.model import AtlasReader

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion'
OUT = BASE / 'male-coeliac-cystic-composition-v1'
CO = BASE / 'male-coeliac-four-port-v1'
MASTER = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
PROPOSAL = ROOT / '.local-models/geometry-workbench/jobs/4f6750028a234cccfb8bd48f92a16f19e9976546d154a9690da7db57b883326a/prepared-1/proposal.json'
PID = 'TA98HEPATIC:male:artery'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def face_keys(mesh):
    """Exact oriented triangles and corresponding literal corner attributes."""
    keys, attributes = [], []
    for face, triangle in zip(mesh.faces, mesh.triangles):
        start = min(range(3), key=lambda i: tuple(triangle[i]))
        order = np.roll(np.arange(3), -start)
        keys.append(tuple(triangle[order].ravel()))
        attributes.append(tuple((name, tuple(np.asarray(value[face][order]).ravel()))
                                for name, value in sorted(mesh.attributes.items())))
    return keys, attributes


def face_owners(record, count):
    labels = [()] * count
    for group in record.get('tissueDrawGroups', []):
        start, n = group['start'] // 3, group['count'] // 3
        tag = tuple(sorted((k, str(v)) for k, v in group.items()
                           if k not in ('start', 'count')))
        for i in range(start, start + n):
            labels[i] = tag
    return labels


def main():
    started = time.monotonic()
    reader = AtlasReader(ROOT)
    source_path = CO / 'literal-source-records-before-lingular.json'
    source = json.loads(source_path.read_text())
    source['chunks'] = {int(k): v for k, v in source['chunks'].items()}
    current = json.loads(MASTER.read_text())
    coeliac = json.loads(PROPOSAL.read_text())
    datasets = {}
    for name, atlas in [('source', source), ('coeliac', coeliac), ('current', current)]:
        record = next(p for p in atlas['parts'] if p['id'] == PID)
        mesh = reader.read_part(atlas, record, frame='model:local-ta98:assembled-metres')
        keys, attrs = face_keys(mesh)
        assert len(keys) == len(set(keys)), f'{name}: duplicate oriented face identity'
        datasets[name] = dict(record=record, mesh=mesh, keys=keys, attrs=attrs,
                              lookup={k: i for i, k in enumerate(keys)},
                              owners=face_owners(record, len(mesh.faces)))
    source_keys = datasets['source']['keys']
    source_lookup = datasets['source']['lookup']
    arrays, reports = {}, {}
    for name in ['coeliac', 'current']:
        data = datasets[name]
        mapped = np.asarray([source_lookup.get(key, -1) for key in data['keys']], dtype=np.int64)
        retained = set(mapped[mapped >= 0].tolist())
        removed = np.asarray(sorted(set(range(len(source_keys))) - retained), dtype=np.int64)
        generated = np.flatnonzero(mapped < 0)
        attr_changed = np.asarray([i for i, j in enumerate(mapped) if j >= 0 and
                                  data['attrs'][i] != datasets['source']['attrs'][j]], dtype=np.int64)
        ownership_changed = np.asarray([i for i, j in enumerate(mapped) if j >= 0 and
                                       data['owners'][i] != datasets['source']['owners'][j]], dtype=np.int64)
        arrays[name + '_source_face_ancestry'] = mapped
        arrays[name + '_removed_source_faces'] = removed
        arrays[name + '_generated_faces'] = generated
        arrays[name + '_retained_attribute_changed_faces'] = attr_changed
        arrays[name + '_retained_ownership_changed_faces'] = ownership_changed
        support = datasets['source']['mesh'].triangles[removed].reshape(-1, 3)
        gen = data['mesh'].triangles[generated].reshape(-1, 3)
        reports[name] = dict(triangles=len(data['mesh'].faces),
                             literal_source_faces=len(retained),
                             removed_source_faces=len(removed), generated_faces=len(generated),
                             literal_face_attributes_changed=len(attr_changed),
                             literal_face_ownership_changed=len(ownership_changed),
                             removed_source_bounds=[support.min(0).tolist(), support.max(0).tolist()],
                             generated_bounds=[gen.min(0).tolist(), gen.max(0).tolist()])
    co_removed = set(arrays['coeliac_removed_source_faces'].tolist())
    cy_removed = set(arrays['current_removed_source_faces'].tolist())
    intersection = sorted(co_removed & cy_removed)
    common = sorted(set(range(len(source_keys))) - co_removed - cy_removed)
    arrays['changed_source_face_intersection'] = np.asarray(intersection, np.int64)
    arrays['both_literal_source_faces'] = np.asarray(common, np.int64)
    co_support = datasets['source']['mesh'].triangles[sorted(co_removed)].reshape(-1, 3)
    cy_support = datasets['source']['mesh'].triangles[sorted(cy_removed)].reshape(-1, 3)
    d, _ = cKDTree(cy_support).query(co_support)
    OUT.mkdir(parents=True, exist_ok=True)
    arrays_path = OUT / 'source-support-ancestry-v1.npz'
    np.savez_compressed(arrays_path, **arrays)
    evidence = dict(schema='coeliac-cystic-common-source-support/v1',
                    source_catalogue_record_archive=dict(path=str(source_path), sha256=sha(source_path)),
                    current_catalogue=dict(path=str(MASTER), sha256=sha(MASTER)),
                    coeliac_proposal=dict(path=str(PROPOSAL), sha256=sha(PROPOSAL)),
                    datasets={name: dict(bank_path=d['mesh'].source['bank_path'],
                                         bank_sha256=d['mesh'].source['bank_sha256'],
                                         geometry_sha256=d['mesh'].geometry_sha,
                                         attributes_sha256=d['mesh'].attributes_sha,
                                         triangles=len(d['mesh'].faces))
                              for name, d in datasets.items()},
                    changes=reports,
                    common_literal_source_faces=len(common),
                    changed_source_face_intersection_count=len(intersection),
                    changed_source_support_minimum_vertex_separation_m=float(d.min()),
                    ancestry_arrays=dict(path=str(arrays_path), sha256=sha(arrays_path)),
                    geometry_emitted=False, installation_allowed=False,
                    qualification='Exact source triangles, corner attributes and ownership only. Disjoint face support does not alone establish a valid compositional seam or tissue/lumen acceptance.',
                    elapsed_seconds=time.monotonic() - started)
    output = OUT / 'source-support-evidence-v1.json'
    output.write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps(dict(output=str(output), sha256=sha(output), changes=reports,
                          intersection=len(intersection), common=len(common),
                          minimum_support_vertex_separation_m=float(d.min()),
                          elapsed_seconds=evidence['elapsed_seconds']), indent=2))


if __name__ == '__main__':
    main()
