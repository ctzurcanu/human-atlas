"""Pack audited, render-only section classifiers without closing visible hila.

The existing liver surface and every adjoining organ stay byte-identical.
Closed source lobe volumes classify the cut face; separately closed conduit
envelopes subtract their union. Computational mouth disks never enter the
display mesh. This is not a new tissue partition or additional TA98 coverage.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    start = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/liver-completion'
    output = Path(args.output).resolve()
    output.mkdir(exist_ok=False)
    catalogue = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    raw = catalogue.read_bytes()
    atlas = json.loads(raw)
    old = base / 'shared-lobes-v54/liver.bin'
    assert hashlib.sha256(old.read_bytes()).hexdigest() == atlas['liverIntegration']['bankSha256']
    packed = bytearray(old.read_bytes())
    origin = np.array([0., 1.25, 0.])
    results = {}
    arrays = {}

    def append(a):
        while len(packed) % 4:
            packed.append(0)
        key = (a.dtype.str, a.shape, hashlib.sha256(a.tobytes()).hexdigest())
        if key in arrays:
            return arrays[key]
        offset = len(packed)
        packed.extend(a.tobytes())
        arrays[key] = offset
        return offset

    def volume(name, path):
        source = Mesh.load(path)
        positions = (source.vertices - origin).astype('<f4')
        mesh = Mesh(positions.astype(float) + origin, source.faces, 'model:local-ta98:assembled-metres')
        # Audit the exact Float32 positions the GPU receives, including source
        # computational disks. A passed pre-quantization source is insufficient.
        check = audit(mesh, require_closed=True, max_pairs=2000000, max_tests=16000000)
        results[name] = dict(source=str(path), sourceGeometrySha256=source.geometry_sha,
                             emittedGeometrySha256=mesh.geometry_sha, audit=check,
                             encodingErrorM=float(np.max(np.linalg.norm(mesh.vertices-source.vertices, axis=1))))
        atomic_json(output / (name+'-check.json'), results[name])
        print(name, check['passed'], check['failures'], flush=True)
        if not check['passed']:
            raise ValueError('Invalid section classifier: ' + name)
        return dict(id=name, positions=append(positions), indices=append(source.faces.astype('<u4')),
                    vertexCount=len(positions), indexCount=int(source.faces.size), origin=origin.tolist(),
                    geometrySha256=mesh.geometry_sha)

    volumes = []
    for name in ['right', 'left', 'quadrate', 'caudate']:
        descriptor = volume(name, base / 'closed-liver-cells-v102' / (name+'.npz'))
        descriptor['partId'] = 'TA98PART:ZA:Liver:' + name
        volumes.append(descriptor)
    exclusions = []
    for role in ['artery', 'portal', 'biliary', 'hepatic_vein']:
        exclusions.append(volume(role, base/'local-course-clearance-v106'/(role+'-outer.npz')))
    # The source hila also meet retained parent conduits, which must not acquire
    # an artificial parenchymal end disk at their connection to these networks.
    exclusions.append(volume('retained-native-parents', base/'shared-liver-flows-v123/reserved-native-parents.npz'))
    atomic_json(output/'checks.json', results)
    bank = output / 'section.bin'
    bank.write_bytes(packed)
    bank_url = '/local-models/' + str(bank.relative_to(ROOT/'.local-models'))
    marker = dict(schema='liver-section-volumes/v1', bankUrl=bank_url,
                  bankSha256=hashlib.sha256(packed).hexdigest(),
                  sourceCatalogueSha256=hashlib.sha256(raw).hexdigest(),
                  volumes=volumes, exclusions=exclusions, renderOnly=True,
                  remainingLimits='Render classifiers preserve the installed source hila; this does not certify a new parenchymal tissue partition or microscopic anatomy.')
    chunk = next(p['chunk'] for p in atlas['parts'] if p['id']=='TA98PART:ZA:Liver:right')
    assert all(p['chunk']==chunk for p in atlas['parts'] if p['id'] in atlas['liverIntegration']['partIds'])
    appended = {**atlas['chunks'][chunk], 'url':bank_url, 'bytes':len(packed)}
    appended.pop('gzip', None)
    new_chunk = len(atlas['chunks'])
    atlas['chunks'].append(appended)
    for part in atlas['parts']:
        if part['id'] in atlas['liverIntegration']['partIds']:
            part['chunk'] = new_chunk
    atlas['liverIntegration']['bankUrl'] = bank_url
    atlas['liverIntegration']['bankSha256'] = marker['bankSha256']
    atlas['liverSectionIntegration'] = marker
    # sectionSolid remains false: the displayed lobes still have open source
    # hila. Only the separately audited classifiers may provide their cut faces.
    atomic_json(output/'candidate.json', atlas)
    atomic_json(output/'evidence.json', dict(passed=True, sourceCatalogueSha256=marker['sourceCatalogueSha256'],
                marker=marker, existingBankPrefixIdentical=packed[:len(old.read_bytes())]==old.read_bytes(),
                addedBytes=len(packed)-len(old.read_bytes()), totalBytes=len(packed),
                elapsedSeconds=time.monotonic()-start, installed=False, addedTa98Terms=0))


if __name__ == '__main__':
    main()
