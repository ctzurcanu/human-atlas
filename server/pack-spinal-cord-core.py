"""Encode and independently audit the shared grey/white compartments."""
from pathlib import Path
import argparse,json,time,sys
import numpy as np
from geometry_tools.common import digest_file,atomic_json,atomic_bytes
from geometry_tools.tissue import pack,normal_words,decode_packed,validate_cells
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output',required=True);ap.add_argument('--sex',choices=['female','male'],default='female');ap.add_argument('--max-pairs',type=int,default=250000);ap.add_argument('--max-tests',type=int,default=2000000);ap.add_argument('--precheck-start',type=int,default=0);ap.add_argument('--precheck-count',type=int,default=0);ap.add_argument('--checkpoints',help='Explicit saved audit stages; geometry, corner attributes and implementation must still match');ap.add_argument('--shared-only',action='store_true',help='Store the decodable bank once, without full-array NPZ duplicates per cell');ap.add_argument('--origin',type=float,nargs=3,default=[0.,1.35,-.05]);ap.add_argument('--emit-stage',action='store_true',help='Save a source-bound encoded bank, explicitly unvalidated and not installable');ap.add_argument('--from-stage',help='Audit an already encoded source-bound bank without rebuilding its attributes');args=ap.parse_args()
checkpoints=Path(args.checkpoints)if args.checkpoints else Path(args.output).with_name(Path(args.output).name+'-checks')
source=Path(args.source).resolve();data=json.loads((source/'cells.json').read_text())
frame='model:local-'+('female-'if args.sex=='female'else'')+'ta98:assembled-metres'
if args.from_stage:
    from types import SimpleNamespace
    import os
    stage=Path(args.from_stage).resolve();pending=json.loads((stage/'pending.json').read_text())
    assert pending['source_sha256']==digest_file(source/'cells.npz')and pending['bank_sha256']==digest_file(stage/'cells.bin')
    assert pending['source_metadata_sha256']==digest_file(source/'cells.json')
    assert pending['frame']==frame and pending['sex']==args.sex
    packed=SimpleNamespace(data=(stage/'cells.bin').read_bytes(),cells=pending['parts'])
    meshes=decode_packed(packed,frame,np.array(pending['origin']))
    if args.precheck_count:
        names=list(meshes)[args.precheck_start:args.precheck_start+args.precheck_count];assert names
        selected={n:meshes[n]for n in names}
    else:names=list(meshes);selected=meshes
    checks=validate_cells(selected,max_pairs=args.max_pairs,max_tests=args.max_tests,checkpoint_directory=checkpoints)
    if args.precheck_count:
        target=Path(args.output).with_name(Path(args.output).name+'-precheck-'+str(args.precheck_start)+'-'+str(len(names))+'.json')
        atomic_json(target,dict(bank_sha256=pending['bank_sha256'],source_sha256=pending['source_sha256'],names=names,checks=checks,total_cells=len(meshes),complete_assembly_check=False,installed=False))
        print(json.dumps(dict(passed=checks['passed'],audited_cells=len(names),failures=checks['failures'],total_bytes=len(packed.data))),flush=True)
        sys.exit(0 if checks['passed']else 1)
    assert checks['passed'],checks['failures']
    target=Path(args.output).resolve();target.mkdir(parents=True,exist_ok=False);os.link(stage/'cells.bin',target/'cells.bin')
    with np.load(source/'cells.npz')as source_arrays:
        positions=source_arrays['positions'];origin=np.array(pending['origin'])
        encoding_error=float(np.linalg.norm((positions-origin).astype('<f4').astype(float)+origin-positions,axis=1).max())
    manifest=dict(pending,schema='shared-finite-tissues/v1',bank=str(target/'cells.bin'),
        arrangement=str(source/'cells.npz'),surfaces={sid:sid for sid in data['surfaces']},cells=data['cells'],
        part_ids={n:p['id']for n,p in pending['parts'].items()},checks=checks,passed=True,
        atlas_modified=False,anatomical_acceptance=False,decoded_cells_available_from_bank=True,
        decoded_cell_files_omitted=True,installable=True,maximum_position_encoding_error_m=encoding_error)
    atomic_json(target/'tissue.json',manifest)
    print(json.dumps(dict(passed=True,cells=len(meshes),bytes=len(packed.data),occupied_checks=len(checks['pair_occupied_volume_checks']))),flush=True);sys.exit(0)
z=np.load(source/'cells.npz')
surfaces={k:z[k]for k in data['surfaces']}
origin=np.array(args.origin);assert np.isfinite(origin).all();local=(z['positions']-origin).astype('<f4')
# All tissues use the same encoded position bank. Exact coincidences created
# by Float32 rounding must be welded globally before normals are evaluated.
# Remove only triangles collapsed by that same weld; full closure, crossing
# and occupied-volume checks still reject any damaged candidate.
encoded,mapping=np.unique(local.astype(float),axis=0,return_inverse=True)
collapsed=0
for sid,faces in surfaces.items():
    faces=mapping[faces];q=encoded[faces];keep=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0
    collapsed+=int((~keep).sum());surfaces[sid]=faces[keep]
    assert len(surfaces[sid]),'A whole tissue interface collapsed at encoding'
words=normal_words(encoded,surfaces);corrected=0
for sid,faces in surfaces.items():
    q=encoded[faces];fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
    bad=np.einsum('tki,ti->tk',words[sid].astype(float),fn)<=0
    if bad.any():
        # Integer rounding of an almost-tangent smooth normal must not make
        # a stored corner point into its own tissue. Use the actual facet.
        facet=np.rint(fn/np.linalg.norm(fn,axis=1)[:,None]*32767).astype('<i2')
        words[sid][bad]=np.broadcast_to(facet[:,None],words[sid].shape)[bad]
        corrected+=int(bad.sum())
part_ids={n:'TA98CORD:'+args.sex+':'+n for n in data['cells']}
if args.emit_stage:
    target=Path(args.output).resolve();target.mkdir(parents=True,exist_ok=False)
    packed=pack_finite_cells(encoded+origin,surfaces,words,data['cells'],part_ids,origin)
    atomic_bytes(target/'cells.bin',packed.data)
    atomic_json(target/'pending.json',dict(schema='unvalidated-spinal-bank-stage/v1',passed=False,installable=False,
        source_sha256=digest_file(source/'cells.npz'),source_metadata_sha256=digest_file(source/'cells.json'),
        bank_sha256=digest_file(target/'cells.bin'),origin=origin.tolist(),frame=frame,sex=args.sex,
        parts=packed.cells,total_bytes=len(packed.data),position_bytes=int(packed.positions.nbytes),normal_bytes=int(packed.normal_words.nbytes),
        positions_encoded_once=True,reciprocal_normals_use_same_words_with_opposite_sign=True,
        construction=dict(data.get('construction',{}),quantized_outward_normal_corners_corrected=corrected,encoded_collapsed_faces_removed=collapsed)))
    print(json.dumps(dict(encoded=True,validated=False,installable=False,cells=len(packed.cells),bytes=len(packed.data))),flush=True);sys.exit(0)
if args.precheck_count:
    packed=pack_finite_cells(encoded+origin,surfaces,words,data['cells'],part_ids,origin)
    meshes=decode_packed(packed,frame,origin);names=list(meshes)[args.precheck_start:args.precheck_start+args.precheck_count];assert names
    checks=validate_cells({n:meshes[n]for n in names},max_pairs=args.max_pairs,max_tests=args.max_tests,checkpoint_directory=checkpoints)
    atomic_json(Path(args.output).with_name(Path(args.output).name+'-precheck-'+str(args.precheck_start)+'.json'),dict(source_sha256=digest_file(source/'cells.npz'),names=names,checks=checks,total_cells=len(meshes),complete_assembly_check=False,installed=False))
    print(json.dumps(dict(passed=checks['passed'],names=names,failures=checks['failures'],total_bytes=len(packed.data))),flush=True)
    sys.exit(0 if checks['passed']else 1)
if args.shared_only:
    from geometry_tools.finite_bank_artifacts import pack_shared_only
    pack=pack_shared_only
result=pack(encoded+origin,surfaces,data['cells'],part_ids,
    'model:local-'+('female-'if args.sex=='female'else'')+'ta98:assembled-metres',args.output,origin=origin.tolist(),max_pairs=args.max_pairs,max_tests=args.max_tests,
    corner_words=words,
    construction={**data.get('construction',{}),'source_cell_bank_sha256':digest_file(source/'cells.npz'),
        'quantized_outward_normal_corners_corrected':corrected,
        'encoded_coincident_position_welding':True,'encoded_collapsed_faces_removed':collapsed,
        'maximum_source_position_encoding_error_m':float(np.linalg.norm(local.astype(float)+origin-z['positions'],axis=1).max())},
    checkpoint_directory=checkpoints)
assert result['passed'],result['checks']['failures']
print(json.dumps(dict(passed=True,cells=len(result['cells']),bytes=result['total_bytes'],
    encoded_positions=result['position_bytes']//12,maximum_encoding_error_m=result['maximum_position_encoding_error_m'],
    total_tissue_volume_mL=sum(result['checks']['cell_volumes_m3'].values())*1e6,installed=False)),flush=True)
