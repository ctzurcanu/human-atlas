"""Stage/reproduce the iliac native-wall batch with iliac metadata labels.

Run each stage under the workspace runtime supervisor. ``construct`` uses the
reviewed native-wall geometry implementation, then finalizes local proposal
metadata without rewriting its binary bank or changing any geometry arrays.
"""
import argparse
import hashlib
import importlib.util
import json
import runpy
import signal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '.local-models/geometry-workbench/circulatory-completion/native-iliac-walls-v1'
PARTS = [('ZA:Common iliac artery.l', 'common-left', 'A12.2.14.001'),
         ('ZA:Common iliac artery.r', 'common-right', 'A12.2.14.001'),
         ('ZA:External iliac artery.l', 'external-left', 'A12.2.16.002'),
         ('ZA:External iliac artery.r', 'external-right', 'A12.2.16.002')]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finalize_metadata():
    path, proof_path = OUT / 'proposed-atlas.json', OUT / 'packing-proof.json'
    atlas = json.loads(path.read_text())
    proof = json.loads(proof_path.read_text())
    before_hash = sha(path)
    bank_hash = sha(OUT / 'native-walls.bin')
    assert bank_hash == proof['bank_sha256']
    changes = []
    for pid, short, _ in PARTS:
        part = next(p for p in atlas['parts'] if p['id'] == pid)
        if 'nativeForegutWall' not in part:
            assert 'nativeIliacWall' in part
            assert part['tissueWallAssembly'] == 'iliac:male:' + short
            continue
        part['nativeIliacWall'] = part.pop('nativeForegutWall')
        assert part['tissueWallAssembly'] == 'foregut:male:' + short
        part['tissueWallAssembly'] = 'iliac:male:' + short
        changes.append(dict(part_id=pid, renamed_field=['nativeForegutWall', 'nativeIliacWall'],
                            assembly_before='foregut:male:' + short,
                            assembly_after='iliac:male:' + short))
    if changes:
        path.write_text(json.dumps(atlas, indent=2) + '\n')
    after_hash = sha(path)
    current_inventory = OUT / 'inventory-current.json'
    authoritative_report = ROOT / 'reports/ta98-iliac-wall-inventory-20261024.json'
    bindings = dict(inventory_current_path=str(current_inventory),
                    inventory_current_sha256=sha(current_inventory),
                    authoritative_inventory_path=str(authoritative_report),
                    authoritative_inventory_sha256=sha(authoritative_report))
    if changes:
        proof['pre_metadata_catalogue_sha256'] = before_hash
    proof['proposed_catalogue_sha256'] = after_hash
    proof.update(bindings)
    proof['construction_helper_path'] = str(Path(__file__).resolve())
    proof['construction_helper_sha256'] = sha(__file__)
    proof['metadata_namespace'] = 'iliac'
    proof_path.write_text(json.dumps(proof, indent=2) + '\n')
    provenance = dict(schema='native-iliac-local-metadata-finalization/v1',
                      pre_metadata_catalogue_sha256=proof.get('pre_metadata_catalogue_sha256', before_hash),
                      proposed_catalogue_sha256=after_hash, changes=changes,
                      bank_sha256=bank_hash, geometry_and_bank_rewritten=False,
                      live_catalogue_changed=False, construction_helper_path=str(Path(__file__).resolve()),
                      construction_helper_sha256=sha(__file__), **bindings)
    (OUT / 'metadata-finalization.json').write_text(json.dumps(provenance, indent=2) + '\n')
    assert sha(OUT / 'native-walls.bin') == bank_hash
    print(json.dumps(provenance))


def construct():
    path = ROOT / 'server/build-native-foregut-walls-v1.py'
    spec = importlib.util.spec_from_file_location('native_wall_geometry_implementation', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUT, module.PARTS = OUT, PARTS
    module.construct()
    finalize_metadata()


if __name__ == '__main__':
    signal.alarm(35)
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['source', 'construct', 'finalize-metadata'])
    args = parser.parse_args()
    if args.stage == 'source':
        runpy.run_path(str(ROOT / '.local-models/geometry-workbench/circulatory-completion/prepare-native-iliac-walls-v1.py'), run_name='__main__')
    elif args.stage == 'construct':
        construct()
    else:
        finalize_metadata()
