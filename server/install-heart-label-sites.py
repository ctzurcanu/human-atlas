"""Install tested local cardiac references while preserving anatomical banks."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import time

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--stage', required=True)
args = parser.parse_args()
started = time.monotonic()
stage = ROOT / '.local-models/geometry-workbench/heart' / args.stage
rollback = ROOT / '.local-models/ta98-review' / args.stage
assert stage.is_dir() and not rollback.exists()
sha = lambda data: hashlib.sha256(data).hexdigest()
rows = []
for sex in ['male', 'female']:
    target = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas' / sex / 'cardiac-label-sites.json'
    catalogue = target.with_name('atlas.json')
    proof = json.loads((stage / sex / 'label-proof.json').read_text())
    assert proof['source_catalogue_sha256'] == sha(catalogue.read_bytes()), 'Carrier catalogue changed'
    assert proof['all_parts_literal'] and proof['all_chunks_literal'] and proof['fixed_skeleton']
    assert proof['new_textures'] == 0
    if proof.get('typed_geometry_count'):
        assert proof['typed_geometry_count']==71 and proof['new_atlas_part_meshes']==0
        assert proof['new_reference_meshes']==proof['new_meshes']==36
        plan=json.loads((stage/'reference-plan.json').read_text())
        evidence_path=stage.parent/plan['source_geometry_stage']/'float32-audit.json'
        evidence=json.loads(evidence_path.read_text())
        assert evidence['complete']and all(row['passed']for row in evidence['geometries'].values())
        assert proof['float32_audit_complete']and proof['float32_audit_passed']
        assert sha(evidence_path.read_bytes())==proof['float32_audit_sha256']==plan['float32_audit_sha256']
    else:assert proof['new_meshes']==0
    source = (stage / sex / 'cardiac-label-sites.json').read_bytes()
    bundle = json.loads(source)
    registry = bundle['integration']
    assert registry['sex'] == sex and registry['labeledTermCount'] == 113
    assert len(registry['termGroups']) == len(bundle['concepts']) == 71
    assert registry['geometryBackedTermCount'] == 42 and registry['distinctTissueCoverageUnchanged']
    rows.append(dict(model=sex, path=str(target), sha256=sha(source), bytes=len(source),
                     geometry_modified=False,reference_geometry_modified=proof.get('typed_geometry_count',0)==71,typed_reference_count=proof.get('typed_geometry_count',0),geometry_backed_terms=42, qualified_reference_terms=71,
                     labeled_terms=113, source_catalogue_sha256=proof['source_catalogue_sha256'],
                     previous_backup=str(rollback / (sex + '-previous.json')) if target.exists() else None))

rollback.mkdir()
for row in rows:
    if row['previous_backup']:
        Path(row['previous_backup']).write_bytes(Path(row['path']).read_bytes())
manifest = rollback / 'manifest.json'
manifest.write_text(json.dumps(dict(status='prepared', entries=rows), indent=2) + '\n')
installed = []
try:
    for row in rows:
        target = Path(row['path'])
        temp = target.with_name(target.name + '.pending')
        temp.write_bytes((stage / row['model'] / 'cardiac-label-sites.json').read_bytes())
        os.replace(temp, target)
        installed.append(row)
        assert sha(target.read_bytes()) == row['sha256']
except BaseException:
    for row in installed:
        if row['previous_backup']:
            Path(row['path']).write_bytes(Path(row['previous_backup']).read_bytes())
        else:
            Path(row['path']).unlink()
    manifest.write_text(json.dumps(dict(status='rolled-back', entries=rows), indent=2) + '\n')
    raise
manifest.write_text(json.dumps(dict(status='installed', entries=rows,
                                    elapsed_seconds=time.monotonic()-started), indent=2) + '\n')
print(f'Installed 113 cardiac term controls per sex. Body geometry unchanged. Rollback: {manifest}. '
      f'Elapsed {time.monotonic()-started:.2f}s.', flush=True)
