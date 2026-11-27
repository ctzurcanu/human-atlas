"""Bounded, fail-closed checkpoints for the coupled tongue packer.

File discovery never hashes data. Geometry snapshots every discovered input
before using it; separate invocations hash each deterministic batch exactly
once. Final publication requires completed supervisor reports, the independent
decoded-bank read, every batch proof, and unchanged file identities. Pending
files are deliberately not a catalogue-candidate accepted by the installer.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import stat
import time

from geometry_tools.common import atomic_json, digest_file

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = 'tongue-packing-checkpoint/v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def deadline(seconds=40):
    end = time.monotonic()+seconds
    def remaining():
        if time.monotonic() >= end:
            raise TimeoutError('Tongue packing stage exceeded its finite worker budget')
    return remaining


def fingerprint(path):
    """Include ctime/inode so restoring mtime cannot hide a replacement."""
    s = Path(path).stat()
    require(stat.S_ISREG(s.st_mode), 'Binding is not a regular file: '+str(path))
    return dict(size=s.st_size, mtime_ns=s.st_mtime_ns, ctime_ns=s.st_ctime_ns,
                device=s.st_dev, inode=s.st_ino, mode=s.st_mode)


def discover_files(value):
    """The same string-path coverage as jobs.file_bindings, without hashing."""
    found = set()
    def walk(item):
        if isinstance(item, dict):
            for value in item.values():
                walk(value)
        elif isinstance(item, (list, tuple)):
            for value in item:
                walk(value)
        elif isinstance(item, str) and len(item) < 4096:
            try:
                path = Path(item)
                if path.is_file():
                    found.add(str(path.resolve()))
            except OSError:
                pass
    walk(value)
    return found


def unchanged(rows):
    for row in rows:
        require(fingerprint(row['path']) == row['stat'],
                'A snapshotted input/artifact changed: '+row['path'])


def runtime_identity():
    report = os.environ.get('ATLAS_RUNTIME_REPORT')
    require(report, 'Run each packing stage through geometry_tools.runtime')
    return dict(path=str(Path(report).resolve()), worker_pid=os.getpid())


def completed_runtime(binding):
    report = read(binding['path'])
    require(report.get('schema') == 'bounded-command/v1'
            and report.get('pid') == binding['worker_pid']
            and report.get('state') == 'completed'
            and report.get('returncode') == report.get('worker_returncode') == 0,
            'A packing stage did not finish successfully under its supervisor: '+binding['path'])
    return report


def code_files():
    paths = {p.resolve() for p in (ROOT/'geometry_tools').rglob('*.py')
             if 'tests' not in p.parts}
    paths.update(ROOT/'server'/name for name in (
        'build-tongue-coupled-candidate.py', 'build-tongue-local-candidate.py',
        'tongue_packing_stages.py', 'tongue_basal_contract.py',
        'tongue_whole_input.py', 'prepare-tongue-coupled-install.py'))
    paths.add(ROOT/'app/tissue-wall-material.ts')
    return paths


def begin_geometry(inputs):
    """Snapshot before preflight/geometry; only small code is hashed here."""
    inputs = Path(inputs).resolve()
    initial = fingerprint(inputs)
    plan = read(inputs)
    values = [plan]
    expected = {}
    for row in plan.get('bindings', {}).values():
        path = Path(row['path']).resolve()
        expected[str(path)] = row['sha256']
        if path.suffix == '.json':
            before = fingerprint(path)
            values.append(read(path))
            require(fingerprint(path) == before, 'Input changed while discovering its dependencies')
    for row in plan.get('covering_meshes', {}).values():
        expected[str(Path(row['path']).resolve())] = row['sha256']
    paths = {str(inputs)} | set().union(*(discover_files(v) for v in values))
    codes = {str(p.resolve()) for p in code_files()}
    paths |= codes
    rows = [dict(path=p, stat=fingerprint(p), expected_sha256=expected.get(p), role='code' if p in codes else 'input')
            for p in sorted(paths)]
    require(fingerprint(inputs) == initial, 'Packing plan changed during discovery')
    # Code loaded for this geometry pass is bound now; large data waits for the
    # supervised verification batches instead of being hashed repeatedly.
    for row in rows:
        if row['role'] == 'code':
            row['expected_sha256'] = digest_file(row['path'])
    unchanged(rows)
    return dict(rows=rows, runtime=runtime_identity(), inputs=str(inputs))


def batches(rows, max_files=64, max_bytes=128*1024*1024):
    result, current, size = [], [], 0
    for index, row in enumerate(rows):
        if current and (len(current) >= max_files or size+row['stat']['size'] > max_bytes):
            result.append(current)
            current, size = [], 0
        current.append(index)
        size += row['stat']['size']
    if current:
        result.append(current)
    return result


def seal_geometry(output, snapshot, remaining):
    output = Path(output)
    unchanged(snapshot['rows'])
    rows = list(snapshot['rows'])
    for name in ('geometry.bin', 'pending-atlas.json', 'pending-inventory.json', 'decoded-packing-proof.json'):
        remaining()
        path = output/name
        rows.append(dict(path=str(path.resolve()), stat=fingerprint(path),
                         expected_sha256=digest_file(path), role='pending-output'))
    rows.sort(key=lambda row: row['path'])
    require(len(rows) == len({r['path'] for r in rows}), 'Duplicate checkpoint identities')
    remaining()
    manifest = dict(schema='tongue-packing-file-manifest/v1', rows=rows,
                    batches=batches(rows), installed=False)
    manifest_path = output/'file-manifest.json'
    atomic_json(manifest_path, manifest)
    checkpoint = dict(schema=SCHEMA, state='geometry-passed-awaiting-verification',
        runtime=snapshot['runtime'], inputs=snapshot['inputs'],
        manifest_path=str(manifest_path.resolve()), manifest_sha256=digest_file(manifest_path),
        verification_batches=len(manifest['batches']), candidate_published=False, installed=False)
    atomic_json(output/'checkpoint.json', checkpoint)
    return checkpoint


def load_checkpoint(output):
    output = Path(output).resolve()
    checkpoint_path = output/'checkpoint.json'
    checkpoint = read(checkpoint_path)
    require(checkpoint.get('schema') == SCHEMA
            and checkpoint.get('state') == 'geometry-passed-awaiting-verification'
            and checkpoint.get('candidate_published') is False,
            'A complete geometry checkpoint is required; partial old output cannot be resumed')
    completed_runtime(checkpoint['runtime'])
    manifest_path = Path(checkpoint['manifest_path']).resolve()
    require(manifest_path == output/'file-manifest.json'
            and digest_file(manifest_path) == checkpoint['manifest_sha256'], 'File manifest changed')
    manifest = read(manifest_path)
    require(manifest.get('schema') == 'tongue-packing-file-manifest/v1'
            and manifest.get('batches') == batches(manifest['rows'])
            and len(manifest['rows']) == len({r['path'] for r in manifest['rows']})
            and len(manifest['batches']) == checkpoint['verification_batches'],
            'Incomplete or altered checkpoint manifest')
    unchanged(manifest['rows'])
    return checkpoint, manifest, digest_file(checkpoint_path)


def timed_digest(path, remaining):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while True:
            remaining()
            block = stream.read(1024*1024)
            if not block:
                break
            value.update(block)
    return value.hexdigest()


def independent_geometry_read(output, remaining):
    """A new worker decodes every pending physical part, not a copied pass flag."""
    from geometry_tools.audit import topology
    from geometry_tools.constraints import check_constraint
    from geometry_tools.model import AtlasReader
    candidate = read(output/'pending-atlas.json')
    evidence = read(output/'decoded-packing-proof.json')
    inventory = read(output/'pending-inventory.json')
    require(evidence.get('schema') == 'tongue-coupled-pending-geometry/v1'
            and evidence.get('all_input_gates_passed') is True
            and evidence.get('installed') is False, 'Unqualified pending geometry proof')
    marker = candidate.get('finiteTongueIntegration', {})
    require(marker.get('activePhysicalMeshes') == 12 and len(set(marker['partIds'])) == 12
            and len(marker.get('termGroups', [])) == 17
            and inventory.get('term_count') == 36, 'Pending assembly identities changed')
    parts = {p['id']:p for p in candidate['parts']}
    reader, meshes, proofs = AtlasReader(), {}, []
    frame = 'model:local-ta98:assembled-metres'
    for row in evidence['records']:
        remaining()
        pid = row['part_id']
        require(pid not in meshes, 'Duplicate decoded proof')
        actual = reader.read_part(candidate, parts[pid], frame=frame)
        topo = topology(actual)
        require(actual.geometry_sha == row['decoded_geometry_sha256']
                and topo['closed_oriented_manifold']
                and topo['components'] == row['topology']['components'],
                'Independent decoded geometry/topology differs: '+pid)
        meshes[pid] = actual
        proofs.append(dict(part_id=pid, geometry_sha256=actual.geometry_sha, topology=topo))
    require(set(meshes) == set(marker['partIds']), 'Independent read did not cover all active meshes')
    interfaces = []
    for pid in sorted(meshes):
        for group in parts[pid].get('tissueDrawGroups', []):
            peer = group.get('neighbor')
            if not peer or pid > peer:
                continue
            remaining()
            require(peer in meshes, 'Missing shared-tissue counterpart')
            faces = lambda a,b: [f for g in parts[a].get('tissueDrawGroups', []) if g.get('neighbor') == b
                for f in range(g['start']//3, (g['start']+g['count'])//3)]
            proof = check_constraint(dict(kind='shared-interface', parameters=dict(
                first_faces=faces(pid,peer), second_faces=faces(peer,pid))), meshes[pid], meshes[peer])
            require(proof['passed'], 'Independent reciprocal chart read failed')
            interfaces.append(dict(first=pid, second=peer, actual_check=proof))
    return dict(schema='tongue-independent-packing-read/v1', passed=True,
                actual_meshes=proofs, actual_shared_interfaces=interfaces, installed=False)


def verify_stage(output, batch):
    began, remaining = time.monotonic(), deadline()
    output = Path(output).resolve()
    checkpoint, manifest, checkpoint_sha = load_checkpoint(output)
    require(0 <= batch < len(manifest['batches']), 'Unknown verification batch')
    proof_path = output/('verified-batch-'+str(batch)+'.json')
    require(not proof_path.exists(), 'A batch is immutable; do not overwrite an earlier attempt')
    rows = []
    for index in manifest['batches'][batch]:
        remaining()
        row = manifest['rows'][index]
        actual_sha = timed_digest(row['path'], remaining)
        require(not row.get('expected_sha256') or actual_sha == row['expected_sha256'],
                'Producer/code/output hash differs: '+row['path'])
        require(fingerprint(row['path']) == row['stat'], 'File changed while hashing: '+row['path'])
        rows.append(dict(index=index, path=row['path'], sha256=actual_sha, stat=row['stat']))
    decoded = None
    if batch == 0:
        decoded = independent_geometry_read(output, remaining)
    unchanged(manifest['rows'])
    remaining()
    proof = dict(schema='tongue-packing-batch-verification/v1', passed=True, batch=batch,
        checkpoint_sha256=checkpoint_sha, manifest_sha256=checkpoint['manifest_sha256'],
        files=rows, independent_geometry_read=decoded, runtime=runtime_identity(),
        elapsed_seconds=time.monotonic()-began, installed=False)
    atomic_json(proof_path, proof)
    return dict(verified_batch=batch, files=len(rows), total_batches=len(manifest['batches']),
                independent_geometry_read=decoded is not None, installed=False)


def finalized_files(output, checkpoint, manifest, checkpoint_sha):
    bindings, proofs = {}, []
    for batch, indices in enumerate(manifest['batches']):
        path = output/('verified-batch-'+str(batch)+'.json')
        proof = read(path)
        require(proof.get('schema') == 'tongue-packing-batch-verification/v1'
                and proof.get('passed') is True and proof.get('batch') == batch
                and proof.get('checkpoint_sha256') == checkpoint_sha
                and proof.get('manifest_sha256') == checkpoint['manifest_sha256']
                and [r['index'] for r in proof.get('files', [])] == indices,
                'Missing or altered verification batch: '+str(batch))
        completed_runtime(proof['runtime'])
        for actual, index in zip(proof['files'], indices):
            row = manifest['rows'][index]
            require(actual['path'] == row['path'] and actual['stat'] == row['stat']
                    and len(actual.get('sha256', '')) == 64
                    and (not row.get('expected_sha256') or actual['sha256'] == row['expected_sha256']),
                    'Batch does not verify its actual manifest member')
            bindings[row['path']] = actual['sha256']
        if batch == 0:
            require(proof.get('independent_geometry_read', {}).get('passed') is True,
                    'Missing independent decoded-bank verification')
        proofs.append(dict(path=str(path), sha256=digest_file(path)))
    require(len(bindings) == len(manifest['rows']), 'Not every input/code/artifact was hashed')
    unchanged(manifest['rows'])
    return bindings, proofs


def finalize_stage(output):
    began, remaining = time.monotonic(), deadline()
    output = Path(output).resolve()
    checkpoint, manifest, checkpoint_sha = load_checkpoint(output)
    require(not any((output/name).exists() for name in ('atlas.json', 'inventory.json', 'packing-evidence.json', 'finalization.json')),
            'Final publication is immutable; partial publication is not reusable')
    bindings, proofs = finalized_files(output, checkpoint, manifest, checkpoint_sha)
    evidence = read(output/'decoded-packing-proof.json')
    evidence.update(schema='tongue-coupled-packed-candidate/v1', input_bindings=bindings,
        staged_packing=dict(schema=SCHEMA, checkpoint_sha256=checkpoint_sha,
            manifest_sha256=checkpoint['manifest_sha256'], verification_proofs=proofs,
            independent_geometry_read_passed=True, runtime=runtime_identity()),
        elapsed_seconds=time.monotonic()-began)
    remaining()
    # Publish only after every gate above. A caller must additionally check the
    # final supervisor record; a timeout after a write can never qualify it.
    atomic_json(output/'packing-evidence.json', evidence)
    for source, target in (('pending-inventory.json','inventory.json'), ('pending-atlas.json','atlas.json')):
        (output/target).write_bytes((output/source).read_bytes())
    atomic_json(output/'finalization.json', dict(schema='tongue-packing-finalization/v1',
        passed=True, checkpoint_sha256=checkpoint_sha, runtime=evidence['staged_packing']['runtime'],
        files={name:digest_file(output/name) for name in ('atlas.json','inventory.json','packing-evidence.json')},
        installed=False))
    return dict(output=str(output), candidate_published=True, verification_batches=len(proofs),
                bound_files=len(bindings), installed=False, elapsed_seconds=time.monotonic()-began)


def require_finalized_stage(candidate):
    output = Path(candidate).resolve().parent
    final = read(output/'finalization.json')
    require(final.get('schema') == 'tongue-packing-finalization/v1' and final.get('passed') is True,
            'A fully finalized staged packing candidate is required')
    completed_runtime(final['runtime'])
    checkpoint, manifest, checkpoint_sha = load_checkpoint(output)
    require(final['checkpoint_sha256'] == checkpoint_sha, 'Finalization belongs to another checkpoint')
    finalized_files(output, checkpoint, manifest, checkpoint_sha)
    require(set(final.get('files', {})) == {'atlas.json','inventory.json','packing-evidence.json'}, 'Incomplete publication')
    for name, sha in final['files'].items():
        require(digest_file(output/name) == sha, 'Final publication changed: '+name)
    return final
