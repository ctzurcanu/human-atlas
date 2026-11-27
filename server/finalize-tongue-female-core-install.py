"""Finish the journal of the reviewed female core after a stopped derived rebuild.

This never applies, refits, rolls back, or edits anatomy. A complete source-bound
forest and literal rollback backups are mandatory before any installed DB state.
"""
from pathlib import Path
import argparse, fcntl, gzip, hashlib, importlib.util, json, sys, time
import numpy as np
from psycopg.types.json import Jsonb

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from geometry_tools.common import CACHE, atomic_json, digest_file, identity
from geometry_tools.db import Database
from geometry_tools.install import catalogue_diff
from geometry_tools.jobs import artifacts, implementation_signature, verify_bindings
from geometry_tools.model import AtlasReader
from structural_relation_cache import resolver_inputs_match

FEMALE = '8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'
PREVIOUS = '770f47f802c79cf6f3e4336bf2fa7c42cf710998069221f533a1ea3f8322e470'
MALE = '1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39'
JOB = 'f2a24e1c787acae83b18fcf290af92bfd855e52589b069f54def5213188aa346'
NAME = 'female-tongue-core-partial-v1'
STAGE = CACHE / 'tongue-completion/female-core-only-staged-validation-v1'


def load(path):
    return json.loads(Path(path).read_bytes())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main(finalize, output):
    began = time.monotonic()
    proposal_path = CACHE / 'jobs' / JOB / 'prepared-1/install.json'
    proposal = load(proposal_path)
    manifest = ROOT / '.local-models/ta98-review' / NAME / 'manifest.json'
    records = load(manifest)
    require(len(records) == 17 and len({r['path'] for r in records}) == 17,
            'Incomplete or repeated rollback manifest')
    backups = {}
    for row in records:
        backup = ROOT / row['backup']
        require(digest_file(backup) == row['sha256'], 'Rollback backup changed')
        backups[row['path']] = gzip.decompress(backup.read_bytes())
    require(proposal['job_id'] == JOB and proposal['model'] == 'local-female-ta98'
            and proposal['proposed_sha256'] == FEMALE
            and proposal['expected_catalogue_sha256'] == PREVIOUS,
            'Wrong reviewed female proposal')
    proof_path = STAGE / 'independent-final-proof-v1.json'
    proof = load(proof_path)
    require(proof['passed'] and proof['final_job_id'] == JOB
            and proof['corrected_unique_constraints92'] == 92,
            'Independent complete validation proof missing')
    verify_bindings(proof['bindings'])
    interrupted_path = CACHE / 'runtime/a32c8264ac644270884c139d019acf53.json'
    interrupted = load(interrupted_path)
    require(interrupted['state'] == 'timed_out' and interrupted['returncode'] == 124
            and str(proposal_path.relative_to(ROOT)) in interrupted['command'] and NAME in interrupted['command'],
            'Exact interrupted apply evidence is missing')
    verify_bindings({proposal['proposed_catalogue']: FEMALE,
                     proposal['inventory']: proposal['inventory_sha256'],
                     **proposal['bank_bindings']})
    require(implementation_signature() == proposal['validation_implementation'],
            'Geometry validation implementation changed')
    db = Database()
    with db.connect() as c:
        c.execute("SET LOCAL statement_timeout = '10000ms'")
        job = c.execute('SELECT * FROM atlas_geom.jobs WHERE id=%s', (JOB,)).fetchone()
    require(job and job['state'] == 'reviewed'
            and job['input_signature'] == proof['final_input_signature']
            and job['result']['passed'], 'Exact reviewed job is required')
    actual_artifacts = artifacts(db, JOB)
    require(set(proposal['bank_bindings'].values()) <= {r['sha256'] for r in actual_artifacts},
            'Reviewed artifacts do not bind every physical bank')
    parent = load(ROOT / 'app/data/ta98-parents.json')
    require(parent['hierarchySha256'] == digest_file(ROOT / 'public/assets/ta98.json'),
            'TA98 parent index does not bind the current hierarchy')
    helper_spec = importlib.util.spec_from_file_location('local_revision_recovery', ROOT / 'scripts/local_catalogue_revision.py')
    helper = importlib.util.module_from_spec(helper_spec)
    helper_spec.loader.exec_module(helper)
    old_routes = json.loads(backups['reports/alimentary-visit-paths.json'])
    new_routes = load(ROOT / 'reports/alimentary-visit-paths.json')
    model_proof, derived_proof = [], []
    reader = AtlasReader()
    for sex, expected in [('male', MALE), ('female', FEMALE)]:
        relative = f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json'
        current = ROOT / relative
        before_raw = backups[relative]
        after_raw = current.read_bytes()
        before, after = json.loads(before_raw), json.loads(after_raw)
        require(sha(after_raw) == expected, 'Installed catalogue differs from reviewed bytes')
        require(gzip.decompress(Path(str(current) + '.gz').read_bytes()) == after_raw,
                'Compressed catalogue is stale')
        if sex == 'male':
            require(before_raw == after_raw, 'Male catalogue changed during female recovery')
        else:
            require(sha(before_raw) == PREVIOUS, 'Female rollback does not contain the reviewed source')
            catalogue_diff(before, after, proposal['changed_parts'])
        visit_key = f'.local-models/ta98-review/visits/{sex}.json'
        old_visit_raw, visit_raw = backups[visit_key], (ROOT / visit_key).read_bytes()
        old_visit, visit = json.loads(old_visit_raw), json.loads(visit_raw)
        fresh = old_visit['sourceCatalogueSha256'] == sha(before_raw)
        if fresh:
            require(visit['sourceCatalogueSha256'] == expected, 'Fresh visit binding is stale')
            old_visit.pop('sourceCatalogueSha256'); visit.pop('sourceCatalogueSha256')
            require(old_visit == visit, 'Visit coordinates or controls changed')
        else:
            require(old_visit_raw == visit_raw
                    and helper.stale_visit_geometry_unchanged(before, after, visit),
                    'Stale visit geometry is not literally preserved')
        old_route = next(r for r in old_routes if r['sex'] == sex)
        route = next(r for r in new_routes if r['sex'] == sex)
        require(old_route['routeSha256'] == sha(old_visit_raw)
                and route['routeSha256'] == sha(visit_raw), 'Visit report hash mismatch')
        if not fresh:
            require(old_route == route, 'Stale visit report was rebound')
        folder = ROOT / 'reports/structural-dependencies/current'
        forest_path = folder / f'{sex}-forest.json'
        forest = load(forest_path)
        ids = {p['id'] for p in after['parts']}
        leaves = [p['id'] for r in forest['forest'] for g in r['children'] for p in g['children']]
        require(forest['sourceSha256'] == expected and forest['ownershipComplete']
                and len(leaves) == len(set(leaves)) == len(ids)
                and set(leaves) == set(forest['partToGroup']) == ids,
                'Current complete ownership forest is missing')
        for r in forest['forest']:
            for g in r['children']:
                require(all(forest['partToGroup'][p['id']] == g['id'] for p in g['children']),
                        'Forest ownership leaf/group mismatch')
        rel_path = folder / f'{sex}-anatomical-edges.json'
        relations = load(rel_path)
        require(relations['sourceSha256'] == expected
                and digest_file(rel_path) == forest['anatomicalRelationEvidence']['sha256']
                and len(relations['edges']) == forest['anatomicalRelationEvidence']['edges']
                and all(e['from'] in ids and e['to'] in ids and not e['coMove'] for e in relations['edges']),
                'Complete relation export does not bind the current forest')
        require(resolver_inputs_match(ROOT, relations), 'Relationship resolver input bindings are stale')
        if sex == 'female':
            staged = relations.get('stagedProducer', {})
            manifest_path = Path(staged.get('manifest', ''))
            require(staged.get('schema') == 'exact-anatomical-relation-stages/v1'
                    and digest_file(manifest_path) == staged.get('manifestSha256'),
                    'Complete staged relationship manifest is missing')
            relation_manifest = load(manifest_path)
            require(relation_manifest['sourceSha256'] == FEMALE
                    and relation_manifest['partCount'] == staged['partCount'] == len(ids)
                    and len(staged['batchFiles']) == staged['completeOrdinalBatches'] == relation_manifest['batchCount'],
                    'Relationship batches do not cover the entire catalogue')
            for name, expected_hash in staged['producerInputs'].items():
                require(digest_file(ROOT / name) == expected_hash, 'Staged producer changed')
            cursor, exported_edges = 0, []
            part_order = [p['id'] for p in after['parts']]
            for index, binding in enumerate(staged['batchFiles']):
                batch_path = Path(binding['path'])
                require(digest_file(batch_path) == binding['sha256'], 'Relationship batch changed')
                batch = load(batch_path)
                require(batch['complete'] and batch['index'] == index and batch['start'] == cursor
                        and batch['end'] == batch['completedThrough']
                        and batch['manifestSha256'] == staged['manifestSha256'],
                        'Relationship batch traversal is incomplete or noncontiguous')
                part_hash = sha(json.dumps(part_order[cursor:batch['end']], ensure_ascii=False, separators=(',', ':')).encode())
                require(part_hash == batch['partIdsSha256'] == relation_manifest['batchPartIdsSha256'][index],
                        'Relationship batch does not name its literal source interval')
                cursor = batch['end']
                exported_edges.extend(batch['edges'])
            require(cursor == len(ids) and exported_edges == relations['edges'],
                    'Published relation export omits or changes a completed batch')
        spatial = forest['spatialIndex']
        spatial_path = ROOT / spatial['file']
        require(digest_file(spatial_path) == spatial['sha256'], 'Spatial index changed')
        with np.load(spatial_path, allow_pickle=False) as z:
            ordered = z['partIds'].tolist()
            require(ordered == [p['id'] for p in after['parts']] and z['bounds'].shape == (len(ids), 2, 3),
                    'Spatial index has an incomplete source order')
            for pid in after.get('finiteTongueIntegration', {}).get('partIds', []):
                part = next(p for p in after['parts'] if p['id'] == pid)
                mesh = reader.read_part(after, part, frame=f'model:local-{sex}-ta98:assembled-metres', catalogue_path=current)
                require(np.max(np.abs(z['bounds'][ordered.index(pid)] - mesh.bounds)) <= 1e-12,
                        'Tongue spatial bounds do not match actual decoded world geometry: ' + pid)
        for ci, expected_hash in spatial['chunkSha256'].items():
            require(reader.file_sha(reader.asset(after['chunks'][int(ci)]['url'])) == expected_hash,
                    'Forest source bank changed')
        model_proof.append(dict(sex=sex, sourceCatalogueSha256=sha(before_raw),
                                catalogueSha256=expected, visitCoordinatesLiteral=True,
                                staleVisitPreservedWithoutRebinding=not fresh,
                                edit=dict(unchanged=True) if sex == 'male' else
                                dict(proposal_sha256=digest_file(proposal_path), changed_parts=proposal['changed_parts'])))
        derived_proof.append(dict(sex=sex, part_count=len(ids), forest=str(forest_path),
                                  forest_sha256=digest_file(forest_path), relations_sha256=digest_file(rel_path),
                                  spatial_sha256=digest_file(spatial_path)))
    installed_files = {str(ROOT / r['path']): digest_file(ROOT / r['path']) for r in records}
    iid = identity('installation', proposal['model'], FEMALE, JOB)
    result = dict(schema='interrupted-local-install-recovery/v1', id=iid, model=proposal['model'],
                  installed=False, verified=True, catalogue_sha256=FEMALE, job_id=JOB,
                  rollback_manifest=str(manifest), models=model_proof, derived=derived_proof,
                  interrupted_apply=str(interrupted_path), interrupted_apply_sha256=digest_file(interrupted_path),
                  installed_files=installed_files, external_deployment=False,
                  elapsed_seconds=time.monotonic() - began)
    if finalize:
        lock = CACHE / 'install.lock'
        with lock.open('a') as stream:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            verify_bindings(installed_files)
            published = manifest.parent / 'published.json'
            journal = dict(installation_id=iid, files=installed_files, job_id=JOB)
            if published.exists():
                require(load(published) == journal, 'A different published recovery journal already exists')
            else:
                atomic_json(published, journal)
            report = ROOT / 'reports' / (NAME + '.json')
            report_data = dict(schema=NAME, models=model_proof, installed=True,
                               backupManifest=str(manifest.relative_to(ROOT)),
                               recoveredAfterInterruptedDerivedRebuild=True)
            if report.exists():
                require(load(report) == report_data, 'A different recovery install report already exists')
            else:
                atomic_json(report, report_data)
            evidence = dict(proposal=str(proposal_path), proof=model_proof,
                            installed_files=installed_files, recovered_after_interrupted_derived_rebuild=True,
                            validation_proof=str(proof_path), validation_proof_sha256=digest_file(proof_path),
                            interrupted_apply=str(interrupted_path), interrupted_apply_sha256=digest_file(interrupted_path),
                            derived_completion=derived_proof)
            with db.connect() as c:
                c.execute("SET LOCAL statement_timeout = '10000ms'")
                status = c.execute('SELECT state FROM atlas_geom.jobs WHERE id=%s FOR UPDATE', (JOB,)).fetchone()
                require(status['state'] == 'reviewed', 'Reviewed status changed before journal finalization')
                c.execute('INSERT INTO atlas_geom.installations(id,job_id,model_id,previous_catalogue_sha256,catalogue_sha256,manifest_path,evidence)VALUES(%s,%s,%s,%s,%s,%s,%s)',
                          (iid, JOB, proposal['model'], PREVIOUS, FEMALE, str(manifest), Jsonb(evidence)))
                c.execute("UPDATE atlas_geom.jobs SET state='installed',updated_at=now()WHERE id=%s", (JOB,))
            result['installed'] = True
    result['elapsed_seconds'] = time.monotonic() - began
    atomic_json(output, result)
    print(json.dumps(dict(installed=result['installed'], verified=True, id=iid,
                          report=str(output), elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--finalize', action='store_true')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    main(args.finalize, args.output.resolve())
