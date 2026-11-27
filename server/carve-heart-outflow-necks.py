"""Remove myocardial outflow plugs while keeping all apparatus cells intact."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.protected_assembly_carve import carve_with_protected_components
from geometry_tools.algorithms.coupled_mask_partition import partition_mask

ap=argparse.ArgumentParser()
ap.add_argument('--bank',required=True)
ap.add_argument('--fluids',required=True)
ap.add_argument('--stage',required=True)
ap.add_argument('--resume-carve')
ap.add_argument('--construction-tolerance-meters',type=float,default=2e-12)
ap.add_argument('--native-boolean-intersections',action='store_true')
ap.add_argument('--source-stage', default='outflow-protected-coupled-v289',
                help='Immutable original source positions for this assembly')
ap.add_argument('--paired-partition', action='store_true')
args=ap.parse_args(); start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart'; out=h/args.stage; out.mkdir(exist_ok=False)
emission=json.loads((h/args.bank/'emission.json').read_text())
assert digest_file(h/args.bank/'cells.bin')==emission['bank_sha256']
assert emission['passed_storage']
source=h/emission['source_stage']; meta=json.loads((source/'ownership.json').read_text())
z=np.load(source/'arrangement.npz'); p=z['positions'].astype('<f4').astype(float)
surfaces={r['id']:z[r['faces']] for r in meta['surfaces']}
normals={r['id']:z[r['normals']] for r in meta['surfaces']}
uvs={r['id']:z[r['uvs']] for r in meta['surfaces']}
field=json.loads((h/args.fluids/'fields.json').read_text())
assert field['passed'] and field['source_bank_sha256']==emission['bank_sha256']
fluids=[Mesh.load(h/args.fluids/f'fluid-{i}.npz') for i in range(2)]
solids=[make_solid(m.vertices,m.faces)[0] for m in fluids]
mutual=float((solids[0]^solids[1]).volume())
assert mutual<1e-16, ('Arterial computational lumens overlap',mutual)
atomic_json(out/'fluid-preflight.json',dict(passed=True,mutual_lumen_overlap_m3=mutual,
    caps_outside_tissue=[r['cap_outside_tissue_proof'] for r in field['records']]))
offset=np.r_[0,np.cumsum([len(m.vertices) for m in fluids])]
faceoffset=np.r_[0,np.cumsum([len(m.faces) for m in fluids])]
v=np.concatenate([m.vertices for m in fluids]); f=np.concatenate([m.faces+offset[i] for i,m in enumerate(fluids)])
caps=np.concatenate([np.load(h/args.fluids/f'fluid-{i}-caps.npz')['cap_face_ids']+faceoffset[i] for i in range(2)])
owners=np.repeat(['aortic-outflow-void','pulmonary-outflow-void'],[len(m.faces) for m in fluids])
print('Carving ten wall cells; keeping the other 28 tissue components intact',flush=True)
def save_carve(result):
    arrays=dict(positions=result['positions']); rows=[]
    for i,(sid,fs) in enumerate(result['surfaces'].items()):
        arrays[f'f{i}']=fs; arrays[f'p{i}']=result['charts'][sid]['parentFaces']; arrays[f'b{i}']=result['charts'][sid]['barycentrics']
        rows.append(dict(id=sid,f=f'f{i}',p=f'p{i}',b=f'b{i}'))
    atomic_npz(out/'raw-carve.npz',**arrays)
    atomic_json(out/'raw-carve.json',dict(rows=rows,cells=result['cells'],surfaceRows=result['surfaceRows'],
        evidence=result['evidence'],bank_sha256=emission['bank_sha256'],fluid_sha256=digest_file(h/args.fluids/'fields.json')))
    print('Raw Boolean result saved; checking protected source complements',flush=True)
prepared=None
if args.resume_carve:
    old=h/args.resume_carve; raw=json.loads((old/'raw-carve.json').read_text()); zz=np.load(old/'raw-carve.npz')
    assert raw['bank_sha256']==emission['bank_sha256'] and raw['fluid_sha256']==digest_file(h/args.fluids/'fields.json')
    prepared=dict(positions=zz['positions'],surfaces={r['id']:zz[r['f']] for r in raw['rows']},
        charts={r['id']:dict(parentFaces=zz[r['p']],barycentrics=zz[r['b']]) for r in raw['rows']},
        cells=raw['cells'],surfaceRows=raw['surfaceRows'],evidence=raw['evidence'])
def boundary_failure(q,actual,expected,record):
    atomic_npz(out/'failed-parent.npz',positions=q,actual=actual,expected=expected)
    atomic_json(out/'failed-parent.json',record)
if args.paired_partition:
    assert not args.resume_carve
    references=[]
    for i,m in enumerate(fluids):
        chart_file=h/args.fluids/f'fluid-{i}-source-charts.json'
        references.extend(json.loads(chart_file.read_text()) if chart_file.exists() else [None]*len(m.faces))
    def output_checkpoint(value):
        arrays=dict(positions=value['positions']); rows=[]
        for i,(sid,fs) in enumerate(value['surfaces'].items()):
            arrays[f'f{i}']=fs; arrays[f'n{i}']=value['normals'][sid]; arrays[f'u{i}']=value['uvs'][sid]
            rows.append(dict(id=sid,faces=f'f{i}',normals=f'n{i}',uvs=f'u{i}'))
        atomic_npz(out/'paired-checkpoint.npz',**arrays)
        atomic_json(out/'paired-checkpoint.json',dict(cells=value['cells'],surfaces=rows,
            selected_volumes_m3=value['selected_volumes_m3'],validation_pending=True,
            bank_sha256=emission['bank_sha256'],fluid_sha256=digest_file(h/args.fluids/'fields.json')))
    result=partition_mask(p,surfaces,meta['cells'],normals,uvs,field['wall_part_ids'],v,f,
        cap_faces=set(map(int,caps)),tolerance=args.construction_tolerance_meters,
        wall_seconds=105.,checkpoint=lambda pid,proof: print('Paired partition',pid,flush=True),
        output_checkpoint=output_checkpoint,mask_face_sources=references)
else:
    result=carve_with_protected_components(p,surfaces,meta['cells'],normals,uvs,field['wall_part_ids'],v,f,
        cutter_owners=owners,classifier_cap_faces=caps,prepared_carve=prepared,
        save_carve=save_carve,boundary_failure=boundary_failure,
        tolerance=args.construction_tolerance_meters,
        exact_source_replay=not args.native_boolean_intersections,
        require_protected_charts_literal=bool(field.get('protected_source_charts')),
        prove_cap_exclusion_by_output_ancestry=bool(field.get('prove_cap_exclusion_by_output_ancestry')))
arrays=dict(positions=result['positions']); definitions=[]
for i,(sid,fs) in enumerate(result['surfaces'].items()):
    arrays[f'surface_{i}_faces']=fs
    arrays[f'surface_{i}_normals']=result['normals'][sid]
    arrays[f'surface_{i}_uvs']=result['uvs'][sid]
    old_sid=sid.rsplit('/mask-',1)[0]
    old_role=next((r.get('role','')for r in meta['surfaces'] if r['id']==old_sid),'')
    definitions.append(dict(id=sid,faces=f'surface_{i}_faces',normals=f'surface_{i}_normals',uvs=f'surface_{i}_uvs',role=old_role))
atomic_npz(out/'arrangement.npz',**arrays)
# Only construction controls can receive a bounded Float32 correction later.
original=np.load(h/args.source_stage/'source-arrangement.npz')['positions']
assert np.array_equal(result['positions'][:len(original)],original)
atomic_npz(out/'source-arrangement.npz',positions=original)
meta.update(source_stage=args.stage,surfaces=definitions,cells=result['cells'],outflow_fluid_stage=args.fluids,
    outflow_carve_evidence=result['evidence'],source_bank_sha256=emission['bank_sha256'],
    numerical_acceptance_pending=True,installed=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'ownership.json',meta)
print(json.dumps(dict(constructed=True,cells=len(result['cells']),
    exposed_protected_patches=len(result['evidence'].get('protected_exposures',[])),
    elapsed_seconds=time.monotonic()-start,installed=False)),flush=True)
