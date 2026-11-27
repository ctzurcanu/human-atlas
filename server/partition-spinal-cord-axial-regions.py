"""Non-overlapping root-registered regional tissue, never duplicate envelopes.

Partition only at boundaries needed by named whole-cord controls. A segment
table is retained for later nerve-root completion, without fabricating 31
independently verified root interfaces from a relative-length template.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import manifold3d as md
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.boolean_cells import BooleanCells

ap=argparse.ArgumentParser();ap.add_argument('--regions',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
regions_path=Path(args.regions).resolve();registration=json.loads(regions_path.read_text())
segments={r['name']:r['y_range_m']for r in registration['segments']}
b,stored=BooleanCells.load(base/'female-regions-v2/boolean-source.npz');md.Manifold.reserve_ids(max(b.primitives)+1)
sources={name:b.derived_source(s,name+'-native-cell',retain_parents=True)for name,s in stored.items()}
# Seven real shared planes make eight axial bands. The same tissues are
# regrouped for cervical/thoracic/lumbar/sacral parts and both enlargements.
bands=[('coccygeal',1.13,segments['Co1'][1]),
       ('sacral-distal',segments['Co1'][1],segments['S2'][0]),
       ('sacral-proximal',segments['S2'][0],segments['S1'][1]),
       ('lumbar',segments['L5'][0],segments['L1'][1]),
       ('thoracic',segments['T12'][0],segments['T1'][0]),
       ('thoracic-t1',segments['T1'][0],segments['T1'][1]),
       ('cervical-enlargement',segments['C8'][0],segments['C3'][1]),
       ('cervical-upper',segments['C3'][1],1.58)]
for a,c in zip(bands,bands[1:]):assert abs(a[2]-c[1])<1e-12,(a,c)
solids={};attributes={}
for ordinal,(region,low,high)in enumerate(bands):
    print(json.dumps(dict(stage='partition',region=region,elapsed_seconds=time.monotonic()-start)),flush=True)
    # A single Boolean classifier per shared band, not independent offsets.
    cube=md.Manifold.cube((.08,high-low,.3)).translate((-.04,low,-.2))
    mask=b.derived_source(cube,'axial-'+region)
    for name,solid in sources.items():
        cell=solid^mask
        assert str(cell.status())=='Error.NoError'and not cell.is_empty(),(region,name)
        key=name+'@'+region;solids[key]=cell
        attributes[key]=dict(native_region=name,axial_region=region,y_range_m=[low,high])
before=sum(float(s.volume())for s in stored.values());after=sum(float(s.volume())for s in solids.values())
assert abs(after-before)<1e-12,(before,after)
b.save(solids,out/'boolean-source.npz')
print(json.dumps(dict(stage='shared-surface-assembly',cells=len(solids),elapsed_seconds=time.monotonic()-start)),flush=True)
result=b.assemble(solids,progress=lambda i,n,o:print(json.dumps(dict(stage='assemble',source_faces=i,total=n,overlays=o,elapsed_seconds=time.monotonic()-start)),flush=True))
atomic_npz(out/'cells.npz',positions=result['positions'],**result['surfaces'])
atomic_json(out/'cells.json',dict(cells=result['cells'],surfaces={k:len(fs)for k,fs in result['surfaces'].items()},
    attributes=attributes,construction=result['evidence'],bands=bands,registration_sha256=digest_file(regions_path),
    source_envelope_preserved=True,central_canal_is_void=True,whole_volume_difference_m3=after-before,
    boundary_status='Root-registered regional partitions; distal boundaries estimated from published relative segment lengths.',
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(cells=len(solids),surfaces=len(result['surfaces']),volume_difference_m3=after-before,elapsed_seconds=time.monotonic()-start)),flush=True)
