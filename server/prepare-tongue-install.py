"""Prepare a source-bound local male tongue review/install job; never apply it.

The catalogue-candidate job independently audits the actual emitted bank and
the declared interfaces. A separately recorded review of that emitted bank is
still required before the ordinary guarded installer can publish it locally.
"""
from __future__ import annotations
import argparse
import copy
from itertools import combinations
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.constraints import check_constraint
from geometry_tools.jobs import file_bindings
from geometry_tools.model import AtlasReader

WORK = ROOT / '.local-models/geometry-workbench/tongue-completion'
CURRENT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME = 'model:local-ta98:assembled-metres'
MODE = 'macroscopic-muscle-extent-reference-envelope'
WHOLE = 'TA98TONGUE:male:whole'
LEGACY_STYLO = ['LOCAL:male:l_styloglossus:1030', 'LOCAL:male:r_styloglossus:981']
REMAINING = [
    'Partial local review only; no complete or anatomically accepted tongue is claimed.',
    'Intrinsic carriers are approximate macroscopic muscle-extent reference envelopes, not exclusive histological tissue volumes.',
    'Transverse/vertical controls share one physical carrier with illustrative directional fibre fields; it is counted once.',
    'Measured intrinsic/native GG/HG/styloglossus envelope crossings remain explicitly reported; they are not a zero-overlap tissue claim.',
    'Native GG/HG attachment footprints and exact interdigitation remain pending; this job neither modifies nor certifies them.',
    'The styloid proof covers the literal three-face origin patches only; a complete tongue insertion is not certified.',
    'Surface-region boundaries are approximate. Mucosa, body/root partitions, septum, aponeurosis, papillae and other absent terms are not completed.',
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def checked_path(path, sha):
    path = Path(path).resolve()
    require(path.is_relative_to(ROOT) and path.is_file() and digest_file(path) == sha, 'Missing or changed source: '+str(path))
    return path


def interface_plan(intrinsics, new_ids, bone_ids, bone_geometry, origins):
    """Only explicit native origin patches permit coplanar bone contact."""
    checks = [dict(kind='fixed', first=pid, parameters=dict(geometry_sha256=bone_geometry[pid])) for pid in bone_ids]
    checks.extend(dict(kind='nested', first=pid, second='ZA:Tongue') for pid in intrinsics)
    checks.extend(dict(kind='disjoint', first=a, second=b, parameters=dict(occupancy='solid')) for a, b in combinations(intrinsics, 2))
    for pid in new_ids:
        host = origins.get(pid, {}).get('host_part_id')
        for bone in bone_ids:
            parameters = dict(allow_coplanar_contact=True) if bone == host else dict(occupancy='solid')
            checks.append(dict(kind='disjoint', first=pid, second=bone, parameters=parameters))
    for pid, origin in origins.items():
        require(origin['constraint']['kind'] == 'shared-interface', 'Only a literal opposed origin seam is supported')
        require(origin['constraint']['parameters']['second_faces'] == [0, 1, 2], 'Styloid origin must retain its three emitted cap faces')
        checks.append(dict(kind='shared-interface', first=origin['host_part_id'], second=pid, parameters=origin['constraint']['parameters']))
    stylo_ids = sorted(origins)
    checks.extend(dict(kind='disjoint', first=a, second=b, parameters=dict(occupancy='solid')) for a, b in combinations(stylo_ids, 2))
    return checks


def main(candidate_path, inventory_path, output):
    began = time.monotonic()
    deadline = began + 30
    def remaining():
        require(time.monotonic() < deadline, 'Tongue job preparation exceeded its 30 second worker budget')
    candidate_path = candidate_path.resolve()
    if candidate_path.is_dir():
        candidate_path /= 'atlas.json'
    inventory_path, output = inventory_path.resolve(), output.resolve()
    require(candidate_path.is_relative_to(WORK) and inventory_path.is_relative_to(WORK), 'Use explicit local tongue candidate and inventory files')
    require(output.is_relative_to(WORK) and output != WORK and not output.exists(), 'Use a new local tongue job stage')
    candidate, packed_inventory = read(candidate_path), read(inventory_path)
    evidence_path = candidate_path.with_name('packing-evidence.json')
    evidence = read(evidence_path)
    require(evidence.get('schema') == 'tongue-local-packed-candidate/v1' and evidence.get('all_gates_passed') is True and evidence.get('installed') is False, 'A complete, uninstalled packed candidate is required')
    require(evidence.get('representation_mode') == MODE and evidence.get('tissue_ownership_claim') is False and evidence.get('anatomical_acceptance') is False, 'This preparation supports the explicit non-owning local envelope proposal only')
    require(checked_path(evidence['source_catalogue'], evidence['source_catalogue_sha256']) == CURRENT, 'Candidate was not built from the current male catalogue')
    before = read(CURRENT)
    before_sha = digest_file(CURRENT)
    require(before_sha == evidence['source_catalogue_sha256'], 'Current catalogue changed after packing')
    bp, ap = {p['id']:p for p in before['parts']}, {p['id']:p for p in candidate['parts']}
    require(len(bp) == len(before['parts']) and len(ap) == len(candidate['parts']), 'Duplicate physical part IDs')
    require(all(ap.get(pid) == part for pid, part in bp.items()) and set(bp) <= set(ap), 'Every existing physical source record must remain literal')
    require(candidate['chunks'][:len(before['chunks'])] == before['chunks'] and len(candidate['chunks']) == len(before['chunks'])+1, 'Append exactly one bank; preserve every native bank reference')
    require(all(candidate['materials'].get(k) == v for k, v in before['materials'].items()), 'Native material records changed')
    require(candidate.get('defaultHidden') == before.get('defaultHidden') and all(pid in candidate.get('defaultHidden', []) for pid in LEGACY_STYLO), 'Legacy source visibility changed')
    native_concepts = {c['id']:c for c in before['concepts']}
    emitted_concepts = {c['id']:c for c in candidate['concepts']}
    require(all(emitted_concepts.get(cid) == c for cid, c in native_concepts.items()), 'Native concept identities or selection changed')
    marker = candidate.get('finiteTongueIntegration', {})
    require(marker.get('schema') == 'tongue-local/v1' and marker.get('wholeConceptId') == WHOLE and marker.get('representationMode') == MODE and marker.get('tissueOwnershipClaim') is False, 'Candidate lacks the explicit reference-envelope marker')
    groups = marker.get('termGroups', [])
    require(len(groups) == 15 and len({g['term'] for g in groups}) == 15 and len({g['conceptId'] for g in groups}) == 15, 'Expected fifteen identity-distinct controls')
    require(len(set(marker['partIds'])) == len(marker['partIds']) == 11 and emitted_concepts[WHOLE]['elements'] == marker['partIds'], 'Shared carriers must be counted once in the eleven-mesh assembly')
    new_ids = sorted(set(ap)-set(bp))
    intrinsics = sorted(pid for pid in new_ids if ap[pid].get('representationMode') == MODE)
    stylo_ids = sorted(pid for pid in new_ids if pid.startswith('TA98TONGUE:male:styloglossus:'))
    require(len(new_ids) == 6 and len(intrinsics) == 4 and len(stylo_ids) == 2 and set(new_ids) == set(intrinsics+stylo_ids), 'Expected four ROI carriers and two new styloglossus sources')
    require(all(ap[pid].get('license') and ap[pid].get('provenance') for pid in new_ids), 'New source licence/provenance metadata must be retained')
    require(all(ap[pid].get('roiBoundaries') == 'approximate' and ap[pid].get('tissueOwnershipClaim') is False for pid in intrinsics), 'Intrinsic controls must retain approximate, non-owning identity')
    bank = AtlasReader().asset(candidate['chunks'][-1]['url'])
    require(digest_file(bank) == evidence['output_bank_sha256'] == marker['bankSha256'], 'Actual emitted bank changed after packing')
    require(packed_inventory.get('schema') == 'ta98-tongue-local-control-inventory/v1' and packed_inventory.get('term_count') == len(packed_inventory.get('terms', [])) == 36, 'The complete classified 36-term packed inventory is required')
    authority = packed_inventory['authority']
    checked_path(authority['path'], authority['sha256'])
    checked_path(authority['corrections_path'], authority['corrections_sha256'])
    require(all(r.get('representation_class') and r['representation_class'] != 'unclassified' for r in packed_inventory['terms']), 'No terminology item may remain unclassified')
    source_audit_path = WORK/'source-v1/source-audit.json'
    checked_path(source_audit_path, evidence['source_audit_sha256'])
    source_audit = read(source_audit_path)
    bone_rows = {r['partId']:r for r in source_audit['contextParts'] if r['sourcePart'].get('system') == 'skeletal'}
    require(len(bone_rows) == 8 and all(bp[pid] == row['sourcePart'] for pid, row in bone_rows.items()), 'All eight declared skeletal sources must stay fixed')
    reader, meshes, actual_reads, emitted_rows = AtlasReader(), {}, {}, {r['part_id']:r for r in evidence['records']}
    read_ids = set(new_ids) | set(bone_rows) | {'ZA:Tongue'}
    for pid in sorted(read_ids):
        remaining()
        mesh = reader.read_part(candidate, ap[pid], frame=FRAME, catalogue_path=candidate_path)
        meshes[pid] = mesh
        actual_reads[mesh.source['bank_path']] = mesh.source['bank_sha256']
        if pid in new_ids:
            require(mesh.geometry_sha == emitted_rows[pid]['decoded_geometry_sha256'], 'Actual emitted geometry differs from packing evidence')
            topo = topology(mesh)
            require(topo['closed_oriented_manifold'] and len(topo['components']) == 1 and topo['components'][0]['closed_orientable_genus'] == 0, 'Every new carrier must be a closed, single-component source')
        elif pid in bone_rows:
            require(mesh.geometry_sha == bone_rows[pid]['geometrySha256'], 'Native skeleton geometry changed')
    origins = {pid:emitted_rows[pid]['styloid_interface'] for pid in stylo_ids}
    for pid, origin in origins.items():
        proof = check_constraint(origin['constraint'], meshes[origin['host_part_id']], meshes[pid])
        require(proof.get('passed') is True and proof.get('literal_opposed_normal_words') is True, 'Actual styloid origin is not literally opposed in coordinates and normal words')
    constraints = interface_plan(intrinsics, new_ids, sorted(bone_rows), {pid:meshes[pid].geometry_sha for pid in bone_rows}, origins)
    context = evidence.get('retained_context_audit', {})
    require(context.get('passed') is True and context.get('representation_mode') == MODE and context.get('tissue_ownership_claim') is False, 'Complete source-bound macro-envelope context evidence is required')
    context_path = checked_path(context['path'], context['sha256'])
    require(len(context.get('classifications', [])) == 24, 'All four intrinsic × six retained muscle counterpart pairs must be classified')
    for row in context['classifications']:
        require(row['classification'] in ('allowed-macroscopic-envelope-overlap', 'no-measured-overlap'), 'Unsupported overlap classification')
        require(row['classification'] != 'allowed-macroscopic-envelope-overlap' or row.get('primary_reference_url'), 'An allowed reference-envelope crossing needs its primary anatomical reference')
    intrinsic_path = checked_path(evidence['intrinsic_construction'], evidence['intrinsic_construction_sha256'])
    intrinsic_check_path = checked_path(intrinsic_path.with_name('checked-geometry.json'), evidence['intrinsic_checks_sha256'])
    input_paths = {candidate_path, inventory_path, evidence_path, source_audit_path, context_path, intrinsic_path, intrinsic_check_path, ROOT/'server/build-tongue-local-candidate.py', Path(__file__).resolve(), ROOT/'ATTRIBUTION.md'}
    if evidence.get('packing_metadata'):
        input_paths.add(checked_path(evidence['packing_metadata']['path'], evidence['packing_metadata']['sha256']))
    for path, sha in [(WORK/'styloglossus-v2/candidate.json', evidence['styloglossus_candidate_sha256']), (WORK/'styloglossus-v2/context-checks.json', evidence['styloglossus_checks_sha256']), (WORK/'surface-controls-v1/surface-controls.json', evidence['surface_controls_sha256'])]:
        input_paths.add(checked_path(path, sha))
    for row in evidence['records']:
        input_paths.add(checked_path(row['input_npz'], row['input_sha256']))
    bindings = {str(p):digest_file(p) for p in input_paths}
    for p in list(input_paths):
        if p.suffix == '.json':
            bindings.update(file_bindings(read(p)))
    bindings.update(actual_reads)
    bindings[str(bank)] = digest_file(bank)
    bindings[str(CURRENT)] = before_sha
    bindings[authority['path']] = authority['sha256']
    bindings[authority['corrections_path']] = authority['corrections_sha256']
    inventory = copy.deepcopy(packed_inventory)
    inventory.update(schema='ta98-geometry-inventory/v1', model='local-ta98', model_sex='male', root_terms=['A05.1.04.001', 'A04.1.05.001'], source_catalogue=str(CURRENT), source_catalogue_sha256=before_sha, authority_path=authority['path'], authority_sha256=authority['sha256'], authority_corrections_path=authority['corrections_path'], authority_corrections_sha256=authority['corrections_sha256'])
    for row in inventory['terms']:
        row['source_meshes'] = row['source_meshes']['male']
        row['applicable_to_model'] = row.get('required_in_typical_adult', True)
        row['candidate_part_ids'] = (row.get('candidate_control') or {}).get('partIds', [])
        row['anatomical_acceptance'] = False
    semantics = [dict(kind='shared-interface', first=o['host_part_id'], second=pid, scope='Only the literal opposed three-face styloid cap; the complete insertion remains pending') for pid, o in origins.items()]
    semantics.extend(dict(kind='nonexclusive-macroscopic-reference-envelope', **row, tissue_ownership_claim=False, anatomical_acceptance=False) for row in context['classifications'])
    semantics.extend([dict(kind='adjacent-unjoined', scope='Native genioglossus–mandible and hyoglossus–hyoid attachments', verification='pending; retained unchanged and not certified by this proposal'), dict(kind='source-bound-surface-reference', scope='Six whole-tongue facet, point and boundary-curve controls; approximate ROI boundaries, no new solids')])
    spec = dict(model='local-ta98', proposed_catalogue=str(candidate_path), expected_catalogue_sha256=before_sha, changed_parts=new_ids, closed_parts=new_ids, require_closed=True, max_pairs=250000, max_tests=2000000, inventory=str(output/'inventory.json'), interface_checks=constraints,
        assembly=dict(whole_organ=WHOLE, part_ids=marker['partIds'], fixed=['Entire skeleton', 'Every original physical source record, concept and material', 'Native GG/HG positions and pending attachments', 'Entire female model'], adjoining_interfaces=semantics, remaining_limits=REMAINING, representation_mode=MODE, tissue_ownership_claim=False, anatomical_acceptance=False),
        evidence_files=sorted(bindings), explicit_input_bindings=bindings, stage='partial-local-review', external_deployment=False)
    remaining()
    require(digest_file(CURRENT) == before_sha, 'Current catalogue changed during preparation')
    require(all(digest_file(path) == sha for path, sha in bindings.items()), 'A bound input changed during preparation')
    output.mkdir(parents=True)
    atomic_json(output/'inventory.json', inventory)
    atomic_json(output/'source-bindings.json', dict(schema='tongue-install-input-bindings/v1', bindings=bindings, actual_read_banks=actual_reads, candidate_sha256=digest_file(candidate_path), source_catalogue_sha256=before_sha, installed=False, elapsed_seconds=time.monotonic()-began))
    spec['evidence_files'].append(str(output/'source-bindings.json'))
    atomic_json(output/'install-spec.json', spec)
    atomic_json(output/'job-spec.json', dict(kind='catalogue-candidate', parameters=spec))
    print(json.dumps(dict(output=str(output), new_parts=len(new_ids), assembly_meshes=len(marker['partIds']), terms=36, interface_checks=len(constraints), bound_files=len(bindings), installed=False, job_started=False, elapsed_seconds=time.monotonic()-began)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True, help='Explicit packed atlas.json or its directory')
    parser.add_argument('--inventory', type=Path, required=True, help='That candidate’s complete 36-term inventory.json')
    parser.add_argument('--output', type=Path, required=True, help='New local tongue job-stage directory')
    args = parser.parse_args()
    main(args.candidate, args.inventory, args.output)
