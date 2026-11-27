"""Pack the reviewed palatoglossus pair as one economical, uninstalled addon.

The addon is data for a subsequent whole-assembly candidate.  It cannot replace
the coupled geometry/installation audit, and does not mark tongue or mucosal
pillar construction complete.  No catalogue, original mesh, or bank is edited.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'server'))
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.model import AtlasReader, Mesh
from tongue_whole_input import WORK, CATALOGUE, whole_input


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def bound(path):
    path = Path(path).resolve()
    require(path.is_file() and path.is_relative_to(ROOT), 'Expected a local bound input')
    return dict(path=str(path), sha256=digest_file(path))


def load(row):
    path = Path(row['path']).resolve()
    require(path.is_relative_to(ROOT/'.local-models') and digest_file(path) == row['sha256'],
            'Source geometry changed: '+str(path))
    mesh = Mesh.load(path)
    require(not row.get('geometry_sha256') or row['geometry_sha256'] == mesh.geometry_sha,
            'Source geometry identity changed')
    return mesh


def main(source, output, coupling_report=None):
    began = time.monotonic()
    def remaining():
        if time.monotonic()-began > 40:
            raise TimeoutError('Palatoglossus addon exceeded its 40-second worker limit')
    source, output = source.resolve(), output.resolve()
    require(source.is_relative_to(WORK) and source.is_file(), 'Select an explicit local candidate')
    require(output.is_relative_to(WORK) and output != WORK and not output.exists(),
            'Use a new output directory')
    candidate = read(source)
    require(candidate.get('schema') == 'tongue-palatoglossus-candidate/v1'
            and candidate.get('numerical_candidate_pass') is True
            and candidate.get('installed') is False and candidate.get('anatomical_acceptance') is False,
            'Use a passed, uninstalled numerical candidate; no anatomical acceptance is inferred')
    plan_path = Path(candidate['plan']['path'])
    require(digest_file(plan_path) == candidate['plan']['sha256'], 'Construction plan changed')
    plan = read(plan_path)
    require(plan['source_catalogue']['sha256'] == digest_file(CATALOGUE), 'Current catalogue changed')
    whole, _, whole_binding = whole_input(Path(plan['whole']['path']),
        WORK/'basal-whole-protected-v1/fit-result.json')
    palate = load(plan['soft_palate'])
    context = [(r['part_id'], load(r)) for r in plan['fixed_bones']+plan['fixed_adjoining']]
    require(len(context) == 12, 'All eight bones and four adjoining tissues are required')
    spec = importlib.util.spec_from_file_location('tongue_addon_packer',
        ROOT/'server/build-tongue-coupled-candidate.py')
    packer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(packer)
    bank, records, meshes, checks = bytearray(), [], {}, []
    material_id = 'TA98TONGUE:male:material:palatoglossus'
    concept_id = 'TA98TONGUE:male:muscle:palatoglossus'
    limitation = ('Original macroscopic palatoglossus envelope. The complete palatal fan and '
        'lingual termination lie inside the actual native tissues. Histological fibre ownership '
        'and faucial-pillar mucosal continuity remain unfinished; not accepted complete anatomy.')
    input_bindings = {str(source): digest_file(source), str(plan_path): digest_file(plan_path)}
    for side in ('left', 'right'):
        remaining()
        row = candidate['sides'][side]
        mesh = load(dict(path=row['mesh'], sha256=row['sha256'], geometry_sha256=row['geometry_sha256']))
        require(mesh.frame == whole.frame == palate.frame, 'The actual source frames must match')
        numerical = audit(mesh, expected_components=1, expected_genus=[0], require_closed=True,
                          max_pairs=30000, max_tests=200000)
        require(numerical['passed'], 'Selected actual muscle fails its finite geometry audit')
        for name, peer in context:
            remaining()
            _, contact = pair_contacts(mesh, peer, 0., max_pairs=30000, max_tests=200000)
            require(contact['complete_within_pair_budget'] and not contact['crossing_face_pairs']
                    and not contact['coplanar_overlap_face_pairs'], 'Actual fixed context crossing: '+name)
            checks.append(dict(side=side, counterpart=name, summary=contact))
        pid = 'TA98TONGUE:male:palatoglossus:'+side
        template = packer.old.append_mesh(bytearray(), mesh, part_id=pid, concept_id=concept_id,
            name='Palatoglossus ('+side+')', material=material_id, terms=['A05.1.04.110'],
            chunk=0, limitation=limitation)
        template.update(side=side, canonicalTa98Parent='A04.1.06.001',
            histologicalEnthesisComplete=False, mucosalPillarContinuityComplete=False)
        records.append(packer.pack_record(bank, mesh, template, np.array([0., 1.625, 0.])))
        meshes[pid] = mesh
        input_bindings[row['mesh']] = row['sha256']
    _, paired = pair_contacts(*meshes.values(), 0., max_pairs=30000, max_tests=200000)
    require(paired['complete_within_pair_budget'] and not paired['crossing_face_pairs']
            and not paired['coplanar_overlap_face_pairs'], 'The actual pair intersects')
    require(len(bank) <= 40000 and sum(p['indexCount']//3 for p in records) == 1200,
            'The addon must retain its two compact, 600-triangle muscle surfaces')
    remaining()
    output.mkdir(parents=True)
    bank_path = output/'palatoglossus.bin'
    bank_path.write_bytes(bank)
    url = '/local-models/'+str(bank_path.relative_to(ROOT/'.local-models'))
    catalogue = dict(parts=records, chunks=[dict(url=url, bytes=len(bank))])
    reader = AtlasReader()
    decoded = []
    for part in records:
        actual = reader.read_part(catalogue, part, frame=whole.frame)
        expected = meshes[part['id']]
        require(np.array_equal(actual.triangles, expected.triangles), 'Emitted physical triangles changed')
        ids = np.asarray(actual.source['source_vertex_ids'])
        words, signs = packer.old.normal_words(expected)
        require(np.array_equal(actual.attributes['normals'], words[ids])
                and np.array_equal(actual.attributes['normalSigns'], signs[ids]),
                'Actual emitted normal words/signs changed')
        decoded.append(dict(part_id=part['id'], geometry_sha256=actual.geometry_sha,
                            triangles_literal=True, normals_literal=True))
    for row in plan['fixed_bones']+plan['fixed_adjoining']+[plan['whole'], plan['soft_palate']]:
        input_bindings[row['path']] = row['sha256']
    require(all(digest_file(p) == sha for p, sha in input_bindings.items()), 'An input changed during packing')
    coupling = None
    if coupling_report is not None:
        coupling_report = coupling_report.resolve()
        require(coupling_report.is_relative_to(WORK), 'Coupling evidence must stay local')
        coupling = read(coupling_report)
        require(coupling.get('schema') == 'palatoglossus-selected-covering-extrinsic-coupling/v1'
                and coupling.get('measurement_complete') is True
                and coupling.get('complete_within_pair_budget') is True
                and coupling.get('actual_pairs') == len(coupling.get('rows', [])) == 16,
                'Use complete actual covering/extrinsic measurements')
        for row in coupling['rows']:
            expected = meshes.get(row['palatoglossus_part_id'])
            require(expected is not None and expected.geometry_sha == row['palatoglossus_geometry_sha256'],
                    'Coupling evidence measures another PG surface')
        require(all(digest_file(p) == sha for p, sha in coupling['source_bindings'].items()),
                'A measured coupling source changed')
        input_bindings.update(coupling['source_bindings'])
        input_bindings[str(coupling_report)] = digest_file(coupling_report)
    input_bindings[str(Path(__file__).resolve())] = digest_file(__file__)
    addon = dict(schema='tongue-palatoglossus-addon/v1', parts=records, chunks=catalogue['chunks'],
        materials={material_id: packer.old.material([.74,.32,.34])},
        termGroup=dict(key='palatoglossus', parent='muscles', name='Palatoglossus',
            term='A05.1.04.110', conceptId=concept_id, partIds=list(meshes), kind='exact'),
        canonical_ta98_parent='A04.1.06.001', input_bindings=input_bindings,
        whole_input=whole_binding, output_bank=bound(bank_path), actual_decoded=decoded,
        checks=checks, paired=paired, physical_meshes=2, triangles=1200, packed_bytes=len(bank),
        new_images=0, original_catalogue_modified=False, installed=False, packing_passed=True,
        assembly_ready=False, covering_entry_gate_passed=(coupling or {}).get('covering_entry_gate_passed'),
        coupling_report=bound(coupling_report) if coupling_report else None,
        unresolved_covering_pairs=(coupling or {}).get('unresolved_covering_pairs'),
        anatomical_acceptance=False, whole_assembly_audit_required=True,
        remaining_limits=[limitation, 'The addon is not a guarded installation proposal.',
            'Covering entry remains unresolved; a whole-organ containment check is insufficient.'],
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'addon.json', addon)
    print(json.dumps({k:addon[k] for k in ('packing_passed','assembly_ready','physical_meshes','triangles',
        'packed_bytes','unresolved_covering_pairs','installed','whole_assembly_audit_required','elapsed_seconds')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--coupling-report', type=Path)
    args = parser.parse_args()
    main(args.source, args.output, args.coupling_report)
