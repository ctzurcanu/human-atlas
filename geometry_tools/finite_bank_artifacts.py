"""Write a validated shared bank without duplicating its arrays per tissue.

Per-cell meshes remain decodable from the bank and explicit draw records.
Numerical validation uses the same emitted meshes and checkpoint gates as
the regular tissue packer; this only changes redundant artifact storage.
"""
from pathlib import Path
import numpy as np
from .common import atomic_json,atomic_npz,atomic_bytes,digest_file
from .tissue import normal_words,decode_packed,validate_cells
from .algorithms.finite_cell_bank import pack_finite_cells


def pack_shared_only(positions,surfaces,cells,part_ids,frame,output,*,origin=(0.,0.,0.),corner_words=None,max_pairs=250000,max_tests=2000000,construction=None,checkpoint_directory=None):
    output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    origin=np.asarray(origin,float);positions=np.asarray(positions,float);encoded=(positions-origin).astype('<f4').astype(float)
    packed=pack_finite_cells(positions,surfaces,corner_words or normal_words(encoded,surfaces),cells,part_ids,origin)
    meshes=decode_packed(packed,frame,origin)
    checks=validate_cells(meshes,max_pairs=max_pairs,max_tests=max_tests,checkpoint_directory=checkpoint_directory)
    if construction and construction.get('finite_wall_cells_pass')is False:
        checks['passed']=False;checks['failures'].append('finite_wall_cells')
    if not checks['passed']:raise ValueError('Rejected finite bank: '+str(checks['failures']))
    output.mkdir(parents=True);bank=output/'cells.bin';atomic_bytes(bank,packed.data)
    arrays={'positions':positions};manifest_surfaces={}
    for i,(sid,faces)in enumerate(surfaces.items()):
        key='surface_'+str(i);arrays[key]=np.asarray(faces,int);manifest_surfaces[sid]=key
    atomic_npz(output/'arrangement.npz',**arrays)
    manifest=dict(schema='shared-finite-tissues/v1',frame=frame,bank=str(bank),bank_sha256=digest_file(bank),arrangement=str(output/'arrangement.npz'),surfaces=manifest_surfaces,
        cells=cells,part_ids=part_ids,parts=packed.cells,position_bytes=int(packed.positions.nbytes),normal_bytes=int(packed.normal_words.nbytes),total_bytes=len(packed.data),
        positions_encoded_once=True,reciprocal_normals_use_same_words_with_opposite_sign=True,
        maximum_position_encoding_error_m=float(np.linalg.norm(encoded+origin-positions,axis=1).max()),construction=construction,checks=checks,passed=True,
        decoded_cell_files_omitted=True,decoded_cells_available_from_bank=True,atlas_modified=False,anatomical_acceptance=False)
    atomic_json(output/'tissue.json',manifest);return manifest
