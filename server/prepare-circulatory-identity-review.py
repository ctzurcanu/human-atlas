"""Render literal male hand/foot vascular sources in their native bone context.

This is a source-identity review, not a finite-wall or organ-interface acceptance.
No model or source bank is modified.
"""
from pathlib import Path
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader


def main():
    started = time.monotonic()
    source = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    expected = '1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39'
    assert digest_file(source) == expected
    atlas = json.loads(source.read_bytes())
    reader = AtlasReader()
    work = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-digital-source-identity-v1'
    work.mkdir(exist_ok=False)
    parts = {p['id']: p for p in atlas['parts']}
    results = []
    for region, bone in [('hand', 'metacarpal'), ('foot', 'metatarsal')]:
        artery = f'ZA:Dorsal digital arteries of {region}.r'
        ids = [artery] + sorted(p['id'] for p in parts.values()
            if not p.get('suppressed') and p['system'] == 'skeletal'
            and p['id'].startswith('ZA:') and p['id'].endswith('.r:Bone')
            and bone in p['name'].lower())
        assert len(ids) == 6
        entries = []
        for i, pid in enumerate(ids):
            p = parts[pid]
            mesh = reader.read_part(atlas, p, frame='model:local-ta98:assembled-metres', catalogue_path=source)
            output = work / f'{region}-{i}.npz'
            mesh.save(output)
            entries.append(dict(mesh=str(output), name=p['name'], part_id=pid,
                color=[.8,.14,.12,1] if pid == artery else [.72,.68,.55,.35]))
        spec = dict(frame='model:local-ta98:assembled-metres', meshes=entries,
            wall_seconds=45, resolution=800,
            views=[dict(name=region+'-native-dorsal', direction=[0,.5,-1]),
                   dict(name=region+'-native-oblique', direction=[1,.4,-1])])
        spec_path = work / (region+'-review.json')
        atomic_json(spec_path, spec)
        results.append(dict(region=region, spec=str(spec_path), part_ids=ids,
            artery_geometry_sha256=reader.read_part(atlas,parts[artery],frame=spec['frame'],catalogue_path=source).geometry_sha))
    assert digest_file(source) == expected
    receipt = dict(schema='circulation-source-identity-review/v1', source=str(source),
        source_sha256=expected, reviews=results, source_unchanged=True,
        wall_lumen_and_organ_interfaces_accepted=False, elapsed_seconds=time.monotonic()-started)
    atomic_json(work/'packet.json',receipt)
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
