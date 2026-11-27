"""Compose two verified disjoint hepatic edits; emit an unapplied proposal."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import audit
from geometry_tools.constraints import check_constraint

helper_spec = importlib.util.spec_from_file_location('composition_support_audit', Path(__file__).with_name('audit-coeliac-cystic-composition-v1.py'))
helper = importlib.util.module_from_spec(helper_spec)
helper_spec.loader.exec_module(helper)
face_keys, face_owners = helper.face_keys, helper.face_owners

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion'
WORK = BASE / 'male-coeliac-cystic-composition-v1'
OUT = WORK / 'candidate-v4'
MASTER = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
OLD = BASE / 'male-coeliac-four-port-v1/binary-candidate-v2/candidate-atlas.json'
HEPATIC = 'TA98HEPATIC:male:artery'
CYSTIC = 'TA98CYSTICBLOOD:male:artery'
COMMON = 'TA98HEPATIC:male:region:A12.2.12.015'
FRAME = 'model:local-ta98:assembled-metres'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def rle(tags):
    result = []
    for i, tag in enumerate(tags):
        values = dict(tag)
        if result and {k: v for k, v in result[-1].items() if k not in ('start', 'count')} == values:
            result[-1]['count'] += 3
        else:
            result.append(dict(start=i * 3, count=3, **values))
    return result


def main():
    started = time.monotonic()
    assert sha(MASTER) == 'd73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
    seam = json.loads((WORK / 'literal-splice-seam-v2.json').read_text())
    assert seam['source_to_replacement_boundary_coordinates_and_orientations_literal']
    assert seam['current_cystic_generated_faces_removed'] == 0
    support = np.load(WORK / 'source-support-ancestry-v1.npz')
    current = json.loads(MASTER.read_text())
    old = json.loads(OLD.read_text())
    cmap = {p['id']: p for p in current['parts']}
    omap = {p['id']: p for p in old['parts']}
    reader = AtlasReader(ROOT)
    def read(atlas, records, pid):
        return reader.read_part(atlas, records[pid], frame=FRAME)
    cm, om = read(current, cmap, HEPATIC), read(old, omap, HEPATIC)
    ck, ca = face_keys(cm)
    ok, oa = face_keys(om)
    clook, olook = {k: i for i, k in enumerate(ck)}, {k: i for i, k in enumerate(ok)}
    regions = [p['id'] for p in current['parts']
               if p['id'].startswith('TA98HEPATIC:male:region:A12.2.12.')]
    assert len(regions) == 10
    def actual_owners(atlas, records, lookup, count):
        owners = np.full(count, '', dtype='<U80')
        for pid in regions:
            mesh = read(atlas, records, pid)
            for key in face_keys(mesh)[0]:
                i = lookup[key]
                assert not owners[i], 'Overlapping regional ownership'
                owners[i] = pid
        return owners
    current_owner = actual_owners(current, cmap, clook, len(cm.faces))
    old_owner = actual_owners(old, omap, olook, len(om.faces))
    assert np.all(current_owner != '')
    missing = np.flatnonzero(old_owner == '')
    assert np.array_equal(missing, np.arange(40235, 40290))
    # These are the common-hepatic reciprocal annulus, omitted by the old view emitter.
    old_owner[missing] = COMMON
    remove = np.asarray(seam['current_v43_faces_to_remove'], dtype=int)
    assert np.all(current_owner[remove] == COMMON)
    keep = np.flatnonzero(~np.isin(np.arange(len(cm.faces)), remove))
    add = support['coeliac_generated_faces']
    assert len(add) == 138 and np.all(old_owner[add] == COMMON)
    assert set(cm.attributes) == set(om.attributes) == {'normals', 'normalSigns'}
    names = sorted(cm.attributes)
    vertices = [v.copy() for v in cm.vertices]
    attributes = {name: [v.copy() for v in cm.attributes[name]] for name in names}
    def vertex_key(mesh, i):
        return (tuple(mesh.vertices[i]), tuple((name, tuple(np.asarray(mesh.attributes[name][i]).reshape(-1)))
                                               for name in names))
    vertex_lookup = {vertex_key(cm, i): i for i in range(len(cm.vertices))}
    appended_faces = []
    for face in om.faces[add]:
        new = []
        for i in face:
            key = vertex_key(om, i)
            if key not in vertex_lookup:
                vertex_lookup[key] = len(vertices)
                vertices.append(om.vertices[i].copy())
                for name in names:
                    attributes[name].append(om.attributes[name][i].copy())
            new.append(vertex_lookup[key])
        appended_faces.append(new)
    faces = np.r_[cm.faces[keep], np.asarray(appended_faces, np.int64)]
    attrs = {name: np.asarray(values, dtype=cm.attributes[name].dtype)
             for name, values in attributes.items()}
    used, inverse = np.unique(faces, return_inverse=True)
    composed = Mesh(np.asarray(vertices)[used], inverse.reshape(-1, 3), FRAME,
                    {name: value[used] for name, value in attrs.items()})
    vertices, faces, attrs = composed.vertices, composed.faces, composed.attributes
    assert len(faces) == 40884
    owners = np.r_[current_owner[keep], old_owner[add]]
    current_tags = face_owners(cmap[HEPATIC], len(cm.faces))
    old_tags = face_owners(omap[HEPATIC], len(om.faces))
    tags = [current_tags[i] for i in keep] + [old_tags[i] for i in add]
    # Both interfaces retain literal source facets, including their opposite corner words.
    before_to_after = np.full(len(cm.faces), -1, dtype=np.int64)
    before_to_after[keep] = np.arange(len(keep))
    old_to_after = np.full(len(om.faces), -1, dtype=np.int64)
    old_to_after[add] = np.arange(len(keep), len(faces))
    co = read(old, omap, 'ZA:Coeliac trunk')
    cy = read(current, cmap, CYSTIC)
    co_spec = dict(kind='shared-interface', first='ZA:Coeliac trunk', second=HEPATIC,
                   parameters=dict(first_faces=list(range(3715, 3770)),
                                   second_faces=old_to_after[np.arange(40235, 40290)].tolist()))
    cy_spec = dict(kind='shared-interface', first=HEPATIC, second=CYSTIC,
                   parameters=dict(first_faces=before_to_after[np.arange(41512, 41548)].tolist(),
                                   second_faces=list(range(4644, 4680))))
    interfaces = [('coeliac', co_spec, check_constraint(co_spec, co, composed)),
                  ('cystic', cy_spec, check_constraint(cy_spec, composed, cy))]
    assert all(result['passed'] for _, _, result in interfaces)
    OUT.mkdir(exist_ok=False)
    bank = bytearray()
    def append(value):
        bank.extend(b'\0' * (-len(bank) % 4))
        offset = len(bank)
        bank.extend(np.ascontiguousarray(value).tobytes())
        return offset
    origin = np.asarray(cmap[HEPATIC]['geometryOrigin'])
    local = (composed.vertices - origin).astype('<f4')
    assert np.array_equal(local.astype(float) + origin, composed.vertices)
    fields = dict(positions=append(local), normals=append(attrs['normals']),
                  normalSigns=append(attrs['normalSigns']), vertexCount=len(vertices),
                  geometryOrigin=origin.tolist())
    out = copy.deepcopy(current)
    rows = {p['id']: p for p in out['parts']}
    chunk = len(out['chunks'])
    whole = rows[HEPATIC]
    whole.update(fields, chunk=chunk, indices=append(faces.astype('<u4')),
                 indexCount=faces.size, triangles=len(faces), bounds=composed.bounds.tolist(),
                 tissueDrawGroups=rle(tags))
    regional_ids = {}
    for pid in regions:
        take = np.flatnonzero(owners == pid)
        selected = faces[take].astype('<u4')
        region_tags = [tuple((k, v) for k, v in tags[i] if k != 'replacement') for i in take]
        rows[pid].update(fields, chunk=chunk, indices=append(selected), indexCount=selected.size,
                         triangles=len(selected), bounds=Mesh(vertices, selected, FRAME).bounds.tolist(),
                         tissueDrawGroups=rle(region_tags))
        regional_ids[pid] = take
    assert sum(len(v) for v in regional_ids.values()) == len(faces)
    bank_path = OUT / 'hepatic-composed.bin'
    bank_path.write_bytes(bank)
    bank_sha = sha(bank_path)
    out['chunks'].append(dict(url='/local-models/' + str(bank_path.relative_to(ROOT / '.local-models')),
                              bytes=len(bank)))
    for pid in [HEPATIC] + regions:
        rows[pid]['hepaticAttributeSignature'] = bank_sha
    for pid in ['ZA:Coeliac trunk', 'ZA:Left gastric artery', 'ZA:Splenic artery']:
        index = next(i for i, p in enumerate(out['parts']) if p['id'] == pid)
        record = copy.deepcopy(omap[pid])
        record['chunk'] = len(out['chunks'])
        out['chunks'].append(copy.deepcopy(old['chunks'][omap[pid]['chunk']]))
        out['parts'][index] = record
    candidate = OUT / 'candidate-atlas.json'
    write(candidate, out)
    decoded = reader.read_part(out, whole, frame=FRAME)
    assert decoded.geometry_sha == composed.geometry_sha and decoded.attributes_sha == composed.attributes_sha
    # No current cystic object, recipient object, skeleton, heart or unrelated part is changed.
    changed = [HEPATIC] + regions + ['ZA:Coeliac trunk', 'ZA:Left gastric artery', 'ZA:Splenic artery']
    result_rows = {p['id']: p for p in out['parts']}
    assert all(p == result_rows[p['id']] for p in current['parts'] if p['id'] not in changed)
    composed.save(OUT / 'hepatic-composed.npz')
    current_partition = np.load(BASE / 'male-cystic-blood-donor-v1/artery-final-source-frame-v43/candidate/owned/artery/partition.npz')
    old_parent_faces = np.load(BASE / 'male-coeliac-four-port-v1/source-owners-v3/hepatic-parent-faces.npz')['parents']
    composed_keys, composed_attrs = face_keys(composed)
    selectors = dict(current_face_ids=np.r_[keep, np.full(len(add), -1)],
                     coeliac_face_ids=np.r_[np.full(len(keep), -1), add],
                     current_native_face_ids=np.r_[support['current_source_face_ancestry'][keep],
                                                   np.full(len(add), -1)],
                     inherited_source_domain_ids=np.r_[current_partition['parent_source_face_ids'][keep], old_parent_faces[add]],
                     owner_ids=owners)
    for pid, take in regional_ids.items():
        selectors['region_' + pid.rsplit('.', 1)[-1]] = take
        record = rows[pid]
        emitted = reader.read_part(out, record, frame=FRAME)
        # Reader compacts attributes; compare literal corners, independent of unused carrier vertices.
        ek, ea = face_keys(emitted)
        assert ek == [composed_keys[i] for i in take]
        assert ea == [composed_attrs[i] for i in take]
        if pid != COMMON:
            original = read(current, cmap, pid)
            assert emitted.geometry_sha == original.geometry_sha
            assert emitted.attributes_sha == original.attributes_sha
    np.savez_compressed(OUT / 'composition-ancestry.npz', **selectors)
    checks = audit(decoded, self_intersections=False, expected_components=1, require_closed=True)
    write(OUT / 'hepatic-topology.json', checks)
    assert checks['passed'], checks['failures']
    proof = dict(schema='coeliac-cystic-exact-composition/v1', current_catalogue_sha256=sha(MASTER),
                 old_coeliac_catalogue_sha256=sha(OLD), proposed_catalogue=str(candidate),
                 proposed_sha256=sha(candidate), changed_parts=changed,
                 physical_parts=['ZA:Coeliac trunk', 'ZA:Left gastric artery', 'ZA:Splenic artery', HEPATIC],
                 bank=dict(path=str(bank_path), sha256=bank_sha, bytes=len(bank)),
                 hepatic=dict(geometry_sha256=decoded.geometry_sha, attributes_sha256=decoded.attributes_sha,
                              faces=len(faces), retained_current_faces=len(keep), appended_old_coeliac_faces=len(add)),
                 region_counts={pid: len(take) for pid, take in regional_ids.items()},
                 regional_faces_cover_whole_exactly_once=True,
                 corrected_old_omitted_coeliac_annulus_owner=COMMON,
                 interfaces=[dict(name=name, spec=spec, result=result) for name, spec, result in interfaces],
                 all_unrelated_records_literal=True, current_cystic_artery_record_literal=True,
                 topology_passed=True, self_contact_pending=True, context_pending=True,
                 actual_render_review_pending=True, anatomical_acceptance=False,
                 install_allowed=False, installed=False, elapsed_seconds=time.monotonic()-started)
    write(OUT / 'emission.json', proof)
    print(json.dumps(dict(output=str(candidate), emission_sha256=sha(OUT / 'emission.json'),
                          hepatic_faces=len(faces), region_counts=proof['region_counts'],
                          elapsed_seconds=proof['elapsed_seconds']), indent=2))


if __name__ == '__main__':
    main()
