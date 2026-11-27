"""Prepare an ordinary guarded install job for an emitted coupled tongue.

No job is started and nothing is applied.  The installer still independently
audits the emitted banks/interfaces and requires a separate actual render
review.  Current unresolved coverings cannot produce a qualifying candidate.
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
sys.path.insert(0, str(ROOT/'server'))
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.constraints import check_constraint
from geometry_tools.jobs import file_bindings
from geometry_tools.model import AtlasReader
from tongue_basal_contract import APO, FRAME, MODE, MUCO, checked_binding, require
from tongue_packing_stages import require_finalized_stage

WORK = ROOT/'.local-models/geometry-workbench/tongue-completion'
CURRENT = ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
WHOLE = 'TA98TONGUE:male:whole'
GG = {'ZA:Genioglossus muscle.l','ZA:Genioglossus muscle.r'}
HG = {'ZA:Hyoglossus muscle.l','ZA:Hyoglossus muscle.r'}
EXTRINSICS = GG|HG


def read(path):
    return json.loads(Path(path).read_text())


def source_preservation(before, after, new_ids):
    bp, ap = ({p['id']:p for p in c['parts']} for c in (before, after))
    require(len(ap) == len(after['parts']) and set(ap)-set(bp) == set(new_ids),
            'Duplicate or unexpected physical part identities')
    require(after['chunks'][:len(before['chunks'])] == before['chunks']
            and len(after['chunks']) == len(before['chunks'])+2,
            'Preserve every native bank and append one muscle bank plus one literal covering bank')
    geometry_keys = {'positions','normals','normalSigns','uvs','indices','vertexCount','indexCount',
        'chunk','bounds','geometryOrigin','sourceOffset','geometryScale','geometryInstance',
        'geometricAttachment','geometricAttachmentToleranceM','histologicalEnthesisComplete','reviewStatus',
        'attachmentEvidence','storedSourceReplay','material','materialCorrection'}
    for pid, original in bp.items():
        if pid == 'ZA:Tongue':
            require(ap[pid] == dict(original, suppressed=True), 'Original tongue arrays/attributes must remain literal and suppressed')
        elif pid in EXTRINSICS:
            require({k:v for k,v in ap[pid].items() if k not in geometry_keys}
                    == {k:v for k,v in original.items() if k not in geometry_keys},
                    'GG/HG replacement changed identity, licence or provenance')
            family = 'genioglossus' if pid in GG else 'hyoglossus'
            material_id = 'TA98TONGUE:male:material:selected-'+family
            expected_material = dict(color=[189,82,87],roughness=.72,pairedNormalSigns=True,
                tissueAppearance=dict(color=[.74,.32,.34],pattern='fibrous',scale=.0004,relief=.000015))
            require('uvs' not in ap[pid] and ap[pid].get('material') == material_id
                    and after['materials'].get(material_id) == expected_material
                    and ap[pid].get('materialCorrection') == dict(
                        reason='Selected corrected GG/HG producer has no UV buffer; use existing UV-free tissue shader',
                        previousMaterial=original['material'],procedural=True,sampledTextures=False,
                        shader='app/tissue-wall-material.ts',pattern='fibrous'),
                    'Only the explicit UV-free procedural replacement is allowed on four corrected GG/HG parts')
            require(ap[pid].get('geometryOrigin') == [0.,1.625,0.]
                    and ap[pid].get('geometricAttachment')
                    and ap[pid].get('histologicalEnthesisComplete') is False
                    and ap[pid].get('attachmentEvidence'),
                    'GG/HG must retain reviewed origin, actual producer evidence and honest attachment scope')
        else:
            require(ap[pid] == original, 'Unrelated source or skeleton changed: '+pid)
    for key, material in before['materials'].items():
        require(after['materials'].get(key) == material, 'A native material changed')
    bc, ac = ({c['id']:c for c in a['concepts']} for a in (before, after))
    require(all(ac.get(cid) == value for cid, value in bc.items() if cid != 'ZA:Tongue'),
            'An unrelated native concept changed')
    native = copy.deepcopy(bc['ZA:Tongue'])
    native.update(elements=after['finiteTongueIntegration']['partIds'],
                  coverageLimitation=ac['ZA:Tongue'].get('coverageLimitation'))
    require(ac['ZA:Tongue'] == native, 'Legacy tongue alias changed more than its active assembly elements/limitation')
    require(set(after.get('defaultHidden', [])) == set(before.get('defaultHidden', []))|{'ZA:Tongue'},
            'Unrelated default visibility changed')


def main(candidate, output):
    began = time.monotonic()
    def remaining():
        if time.monotonic()-began >= 30:
            raise TimeoutError('Coupled tongue job preparation exceeded its 30 second worker budget')
    candidate, output = candidate.resolve(), output.resolve()
    if candidate.is_dir():
        candidate /= 'atlas.json'
    require(candidate.is_relative_to(WORK) and candidate.is_file(), 'Use an explicit local emitted candidate')
    require(output.is_relative_to(WORK) and output != WORK and not output.exists(), 'Use a new local job stage')
    require_finalized_stage(candidate)
    evidence_path = candidate.with_name('packing-evidence.json')
    evidence, after, before = read(evidence_path), read(candidate), read(CURRENT)
    require(evidence.get('schema') == 'tongue-coupled-packed-candidate/v1'
            and evidence.get('all_input_gates_passed') is True and evidence.get('installed') is False
            and evidence.get('emitted_render_review_required') is True,
            'A qualifying source-bound emitted candidate and separate pending review are required')
    require(evidence.get('source_catalogue_sha256') == digest_file(CURRENT)
            and Path(evidence['source_catalogue']).resolve() == CURRENT,
            'Source catalogue changed since packing')
    require(evidence.get('representation_mode') == MODE and evidence.get('tissue_ownership_claim') is False
            and evidence.get('anatomical_acceptance') is False, 'Do not promote approximate muscle ROIs to accepted tissue')
    bindings = evidence.get('input_bindings', {})
    require(bindings and all(Path(p).is_file() and digest_file(p) == sha for p, sha in bindings.items()),
            'An actual selected input or emitted bank has changed')
    marker = after.get('finiteTongueIntegration', {})
    groups = marker.get('termGroups', [])
    require(marker.get('schema') == 'tongue-local/v1' and marker.get('wholeConceptId') == WHOLE
            and marker.get('legacyWholeAlias') is True and len(groups) == 17
            and len({g['term'] for g in groups}) == 17
            and len(marker.get('partIds', [])) == len(set(marker['partIds'])) == 12
            and 'ZA:Tongue' not in marker['partIds'], 'The seventeen controls must own twelve active, unique physical meshes')
    bp, ap = ({p['id']:p for p in a['parts']} for a in (before, after))
    new_ids = sorted(set(ap)-set(bp))
    intrinsics = evidence.get('intrinsic_part_ids', [])
    origins = evidence.get('styloid_interfaces', {})
    require(len(new_ids) == 8 and len(intrinsics) == 4 and len(origins) == 2
            and set(new_ids) == set(intrinsics)|set(origins)|{MUCO,APO},
            'Expected four ROI, two SG and two finite covering sources')
    source_preservation(before, after, new_ids)
    reader, meshes, actual_reads = AtlasReader(), {}, {}
    source_audit = read(WORK/'source-v1/source-audit.json')
    bones = {r['partId']:r for r in source_audit['contextParts'] if r['sourcePart'].get('system') == 'skeletal'}
    require(len(bones) == 8, 'All eight fixed native bone sources must be explicit')
    rows = {r['part_id']:r for r in evidence['records']}
    for pid in sorted(set(new_ids)|EXTRINSICS|set(bones)):
        remaining()
        mesh = reader.read_part(after, ap[pid], frame=FRAME)
        meshes[pid] = mesh
        actual_reads[mesh.source['bank_path']] = mesh.source['bank_sha256']
        if pid in rows:
            require(mesh.geometry_sha == rows[pid]['decoded_geometry_sha256'], 'Actual emitted geometry changed')
            topo = topology(mesh)
            require(topo['closed_oriented_manifold']
                    and topo['components'] == rows[pid]['topology']['components'],
                    'Actual finite-cell components changed; a valid hollow cell may have two boundary components')
        elif pid in bones:
            require(mesh.geometry_sha == bones[pid]['geometrySha256'], 'Fixed skeleton changed')
    constraints = [dict(kind='fixed',first=pid,parameters=dict(geometry_sha256=meshes[pid].geometry_sha)) for pid in bones]
    for pid in sorted(set(new_ids)|EXTRINSICS):
        for bone in bones:
            if origins.get(pid, {}).get('host_part_id') == bone:
                parameters = dict(allow_coplanar_contact=True)
            elif topology(meshes[pid])['closed_oriented_manifold'] and topology(meshes[bone])['closed_oriented_manifold']:
                parameters = dict(occupancy='solid')
            else:
                # A literal open hyoid partition is a surface, never an invented capped solid.
                parameters = {}
            constraints.append(dict(kind='disjoint',first=pid,second=bone,parameters=parameters))
    for a,b in combinations(intrinsics,2):
        constraints.append(dict(kind='disjoint',first=a,second=b,parameters=dict(occupancy='solid')))
    for pid, origin in origins.items():
        proof = check_constraint(origin['constraint'],meshes[origin['host_part_id']],meshes[pid])
        require(proof.get('passed') and proof.get('literal_opposed_normal_words'), 'Actual SG origin lost literal opposed cortex')
        constraints.append(dict(kind='shared-interface',first=origin['host_part_id'],second=pid,
                                parameters=origin['constraint']['parameters']))
    for a,b in combinations(sorted(origins),2):
        constraints.append(dict(kind='disjoint',first=a,second=b,parameters=dict(occupancy='solid')))
    for pid in set(new_ids)|EXTRINSICS:
        for group in ap[pid].get('tissueDrawGroups', []):
            peer=group.get('neighbor')
            if not peer or pid > peer:
                continue
            require(peer in meshes, 'A reciprocal tissue neighbor is absent')
            faces=lambda first,second:[f for g in ap[first].get('tissueDrawGroups',[]) if g.get('neighbor')==second
                for f in range(g['start']//3,(g['start']+g['count'])//3)]
            constraints.append(dict(kind='shared-interface',first=pid,second=peer,
                parameters=dict(first_faces=faces(pid,peer),second_faces=faces(peer,pid))))
    inventory_path = candidate.with_name('inventory.json')
    inventory = read(inventory_path)
    require(inventory.get('term_count') == len(inventory.get('terms', [])) == 36,
            'Complete classified thirty-six-term inventory is required')
    authority=inventory['authority']
    for name in ('path','corrections_path'):
        checked_binding(dict(path=authority[name],sha256=authority['sha256' if name=='path' else 'corrections_sha256']),ROOT)
    inventory.update(schema='ta98-geometry-inventory/v1',model='local-ta98',model_sex='male',
        root_terms=['A05.1.04.001','A04.1.05.001'],source_catalogue=str(CURRENT),
        source_catalogue_sha256=digest_file(CURRENT),authority_path=authority['path'],authority_sha256=authority['sha256'],
        authority_corrections_path=authority['corrections_path'],authority_corrections_sha256=authority['corrections_sha256'])
    for row in inventory['terms']:
        require(row.get('representation_class') not in (None,'unclassified'), 'An inventory term is unclassified')
        row.update(source_meshes=row['source_meshes']['male'],applicable_to_model=row.get('required_in_typical_adult',True),
                   candidate_part_ids=(row.get('candidate_control') or {}).get('partIds',[]),anatomical_acceptance=False)
    bindings.update(actual_reads)
    for path in (candidate,evidence_path,inventory_path,Path(__file__).resolve(),ROOT/'server/build-tongue-coupled-candidate.py',
                 ROOT/'server/tongue_basal_contract.py',ROOT/'ATTRIBUTION.md'):
        bindings[str(path)]=digest_file(path)
    bindings.update(file_bindings(inventory))
    remaining()
    require(all(digest_file(path)==sha for path,sha in bindings.items()), 'A bound input changed during job preparation')
    output.mkdir(parents=True)
    atomic_json(output/'inventory.json',inventory)
    spec=dict(model='local-ta98',proposed_catalogue=str(candidate),expected_catalogue_sha256=digest_file(CURRENT),
        changed_parts=sorted(set(new_ids)|EXTRINSICS|{'ZA:Tongue'}),closed_parts=sorted(set(new_ids)|EXTRINSICS),require_closed=False,
        max_pairs=250000,max_tests=2000000,inventory=str(output/'inventory.json'),interface_checks=constraints,
        evidence_files=sorted(bindings),explicit_input_bindings=bindings,stage='partial-local-review',external_deployment=False,
        assembly=dict(whole_organ=WHOLE,part_ids=marker['partIds'],fixed=['Entire skeleton','All unrelated source anatomy','Entire female model'],
            representation_mode=MODE,tissue_ownership_claim=False,anatomical_acceptance=False,
            adjoining_interfaces=[dict(kind='source-bound-external-corrected-whole-containment',scope='Actual emitted ROI checks bound in packing evidence; retained suppressed ZA:Tongue is the unmodified original, not the corrected container'),
                dict(kind='source-bound-origin-insertion-proposal',scope='All four selected GG/HG arrays and actual passed producer evidence are bound; histological entheses and complete individual attachment footprints remain unaccepted')],
            remaining_limits=[marker['coverageLimitation'],'An actual emitted-bank surface/section review is required before applying this job.']))
    atomic_json(output/'source-bindings.json',dict(schema='tongue-coupled-install-bindings/v1',bindings=bindings,actual_read_banks=actual_reads,installed=False))
    spec['evidence_files'].append(str(output/'source-bindings.json'))
    atomic_json(output/'install-spec.json',spec)
    atomic_json(output/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec))
    print(json.dumps(dict(output=str(output),interface_checks=len(constraints),new_parts=8,actual_meshes=12,
        controls=17,installed=False,job_started=False,elapsed_seconds=time.monotonic()-began)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    main(args.candidate,args.output)
