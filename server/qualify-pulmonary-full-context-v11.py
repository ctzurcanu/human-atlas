"""Classify the literal current pulmonary context without claiming solid anatomy."""
from pathlib import Path
from collections import Counter
import json
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file, identity
from geometry_tools.model import AtlasReader


def main():
    started = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1'
    context_path = base / 'actual-runtime-context-v9/runtime-context.json'
    context = json.loads(context_path.read_text())
    master_path = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    body_path = ROOT / '.local-models/organ-models/spinal-cord/male-body.json'
    assert digest_file(master_path) == context['sourceCatalogueSha256']
    assert digest_file(body_path) == context['ordinaryBodySha256']
    for filename, expected in context['sourceBindings'].items():
        assert digest_file(ROOT / filename) == expected, filename
    master = json.loads(master_path.read_text())
    body = json.loads(body_path.read_text())
    raw = {p['id']: p for p in master['parts']}
    runtime = {p['id']: p for p in body['parts']}
    family = {p['id'] for p in context['sourceRows']}
    guides = {p['id'] for p in context['fixedGuideRows']}
    reader = AtlasReader()
    banks = {}
    rows = []
    geometry_keys = ('vertexCount', 'indexCount', 'positions', 'indices',
                     'normals', 'normalSigns', 'uvs', 'colors', 'geometryInstance')

    def geometry_record(atlas, part):
        value = {k: part.get(k) for k in geometry_keys}
        value.update(bank_url=atlas['chunks'][part['chunk']]['url'],
                     geometryOrigin=part.get('geometryOrigin', [0, 0, 0]),
                     geometryScale=part.get('geometryScale', 1),
                     sourceOffset=part.get('sourceOffset', [0, 0, 0]))
        return value

    for record in context['context']:
        pid = record['id']
        part = raw[pid]
        row = dict(id=pid, retained=record['retained'],
                   source_record_identity=identity(part),
                   runtime_exclusion_basis='executed ordinary-body and preferred-source pipeline',
                   occupied_volume_established=False)
        if not record['retained']:
            row['role'] = 'excluded-by-executed-runtime-pipeline'
            row['runtime_suppressed'] = record['normalRuntimeSuppressed']
            row['replacement_parts'] = record['replacementParts']
            rows.append(row)
            continue
        assert pid in runtime
        layout = geometry_record(master, part)
        assert layout == geometry_record(body, runtime[pid]), pid
        bank = reader.asset(layout['bank_url'])
        stat = bank.stat()
        banks[str(bank)] = dict(path=str(bank), url=layout['bank_url'],
                                bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
        row.update(bank_path=str(bank), geometry_layout_identity=identity(layout),
                   same_literal_layout_bank_and_world_pose=True,
                   geometry_instance=part.get('geometryInstance'),
                   tissue_declarations={k: part[k] for k in (
                       'pairedTissue', 'sectionSolid', 'tissueWallAssembly',
                       'tissueWallLayer', 'tissueGeometrySchema', 'reviewStatus',
                       'sourceOuterSurfaceOnly', 'sourceLumenPatencyUnverified') if k in part})
        if pid in family:
            row['role'] = 'pulmonary-family-source-requires-moving-or-root-ownership-contract'
        elif pid in guides:
            row['role'] = ('fixed-pulmonary-territory-envelope-for-containment' if
                           part.get('tissueWallLayer') == 'segment-envelope' else
                           'fixed-airway-outer-reference-keep-blood-outside')
        elif (':attachments:' in part['material'] and
              'Muscular insertions' in str(part.get('groups', []))):
            row.update(role='nonowning-attachment-footprint-reference',
                       role_basis=dict(material=part['material'], groups=part['groups']),
                       no_automatic_parent_tissue_exclusion=True)
        elif (part.get('pairedTissue') is True and part.get('tissueDrawGroups') and
              part.get('tissueWallAssembly')):
            row['role'] = 'declared-owned-paired-tissue-cell-requires-actual-cell-audit'
        elif part.get('system') in ('arterial', 'venous'):
            row['role'] = 'fixed-vascular-surface-wall-and-lumen-identity-unverified'
        else:
            row['role'] = 'fixed-anatomical-surface-occupied-volume-unverified'
        row['surface_contact_participation_retained'] = True
        rows.append(row)
    out = base / 'full-context-participation-v11'
    manifest = dict(schema='pulmonary-context-literal-bank-manifest/v11',
                    context_sha256=digest_file(context_path), banks=list(banks.values()),
                    bank_contents_hash_pending=True)
    atomic_json(out / 'bank-manifest.json', manifest)
    result = dict(schema='pulmonary-full-context-physical-participation/v11',
                  context_sha256=digest_file(context_path),
                  source_catalogue_sha256=digest_file(master_path),
                  ordinary_body_sha256=digest_file(body_path),
                  loader_bindings=context['sourceBindings'],
                  counts=dict(Counter(r['role'] for r in rows)), records=rows,
                  retained_records=sum(r['retained'] for r in rows),
                  no_system_or_provider_prefix_exclusion=True,
                  all_retained_surfaces_remain_in_contact_diagnostics=True,
                  metadata_flags_do_not_prove_tissue_volume_or_valid_topology=True,
                  anatomical_registration_not_established_by_shared_world_frame=True,
                  geometry_or_visibility_changed=False, anatomical_acceptance=False,
                  bank_manifest_sha256=digest_file(out / 'bank-manifest.json'),
                  elapsed_seconds=time.monotonic() - started)
    atomic_json(out / 'evidence.json', result)
    print(json.dumps(dict(counts=result['counts'], banks=len(banks),
                          total_bank_bytes=sum(r['bytes'] for r in banks.values()),
                          elapsed_seconds=result['elapsed_seconds'])))


if __name__ == '__main__':
    main()
