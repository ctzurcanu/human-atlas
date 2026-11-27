"""Checkpointed checks of actual emitted finite cells, with no installation.

Each audit/pair result is bound to the same packed bank and implementation.
An interrupted batch retains completed checks; incomplete output stays failed.
"""
from pathlib import Path
import itertools
import json
import time
from types import SimpleNamespace
import numpy as np

from .common import atomic_bytes, atomic_json, digest_file, identity
from .model import Mesh


def stage(positions, surfaces, words, cells, part_ids, frame, origin, output):
    from .algorithms.finite_cell_bank import pack_finite_cells
    from .tissue import decode_packed
    from .jobs import implementation_signature
    output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    packed=pack_finite_cells(positions,surfaces,words,cells,part_ids,origin)
    output.mkdir(parents=True)
    bank=output/'cells.bin';atomic_bytes(bank,packed.data)
    paths={}
    # Decode one chart at a time: each packed cell shares the same large bank.
    # Keeping a full unused copy for every cell defeats that representation.
    for name,row in packed.cells.items():
        mesh=decode_packed(SimpleNamespace(data=packed.data,cells={name:row}),frame,np.asarray(origin,float))[name]
        used,indices=np.unique(mesh.faces,return_inverse=True)
        compact=Mesh(mesh.vertices[used],indices.reshape(-1,3),mesh.frame,{k:v[used]for k,v in mesh.attributes.items()})
        path=output/(part_ids[name].replace(':','_').replace('/','_')+'.npz');compact.save(path)
        paths[name]=dict(path=str(path),sha256=digest_file(path))
    manifest=dict(schema='checkpointed-finite-tissues/v1',frame=frame,bank=str(bank),bank_sha256=digest_file(bank),
                  parts=packed.cells,cells=cells,part_ids=part_ids,geometry_origin=list(origin),mesh_files=paths,
                  implementation_sha256=implementation_signature(),total_bytes=len(packed.data),
                  position_bytes=int(packed.positions.nbytes),normal_bytes=int(packed.normal_words.nbytes),
                  positions_encoded_once=True,reciprocal_normals_use_same_words_with_opposite_sign=True,
                  passed=False,validation_state='pending',atlas_modified=False,anatomical_acceptance=False)
    atomic_json(output/'tissue.json',manifest)
    return manifest


def _read(folder):
    from .jobs import implementation_signature
    folder=Path(folder).resolve();manifest=json.loads((folder/'tissue.json').read_text())
    if digest_file(manifest['bank'])!=manifest['bank_sha256']or implementation_signature()!=manifest['implementation_sha256']:
        raise ValueError('Packed bank or check implementation changed')
    for row in manifest['mesh_files'].values():
        if digest_file(row['path'])!=row['sha256']:raise ValueError('An emitted tissue cell changed')
    return folder,manifest


def _cache_path(folder,manifest,kind,key,max_pairs):
    signature=identity(manifest['bank_sha256'],manifest['implementation_sha256'],kind,key,max_pairs)
    return folder/'validation'/(kind+'-'+signature+'.json')


def _finish(folder,manifest,max_pairs):
    checks={};volumes={};failures=[];missing=[];pairs=[]
    for name in manifest['mesh_files']:
        path=_cache_path(folder,manifest,'cell',name,max_pairs)
        if not path.exists():missing.append('cell:'+name);continue
        row=json.loads(path.read_text());checks[name]=row['audit']
        if not row['passed']:failures.append('cell:'+name)
        else:volumes[name]=row['volume_m3']
    for first,second in itertools.combinations(sorted(volumes),2):
        path=_cache_path(folder,manifest,'pair',[first,second],max_pairs)
        if not path.exists():missing.append('pair:'+first+':'+second);continue
        row=json.loads(path.read_text());pairs.append(row)
        if not row['passed']:failures.append('volume_overlap:'+first+':'+second)
    result=dict(passed=not failures and not missing,failures=failures,missing_checks=missing,cells=checks,
                cell_volumes_m3=volumes,pair_occupied_volume_checks=pairs,
                all_emitted_cells_checked=len(checks)==len(manifest['mesh_files']),anatomical_acceptance=False)
    manifest.update(checks=result,passed=result['passed'],validation_state='passed'if result['passed']else'failed'if failures else'pending')
    atomic_json(folder/'tissue.json',manifest,overwrite=True)
    return dict(passed=result['passed'],state=manifest['validation_state'],checked_cells=len(checks),
                total_cells=len(manifest['mesh_files']),checked_pairs=len(pairs),missing_checks=len(missing),failures=failures)


def check_cells(folder,*,names=None,max_pairs=1000000):
    from .audit import audit
    from .algorithms.finite_tissue_boolean import make_solid
    folder,manifest=_read(folder)
    for name in names or manifest['mesh_files']:
        if name not in manifest['mesh_files']:raise ValueError('Unknown staged cell: '+name)
        path=_cache_path(folder,manifest,'cell',name,max_pairs)
        if path.exists():continue
        began=time.monotonic();mesh=Mesh.load(manifest['mesh_files'][name]['path'])
        result=audit(mesh,require_closed=True,max_pairs=max_pairs);volume=None
        if result['passed']:
            solid,_=make_solid(mesh.vertices,mesh.faces);volume=float(solid.volume())
        row=dict(cell=name,audit=result,volume_m3=volume,passed=bool(result['passed']and volume is not None and volume>0),elapsed_seconds=time.monotonic()-began)
        atomic_json(path,row)
        print(json.dumps(dict(cell=name,passed=row['passed'],elapsed_seconds=row['elapsed_seconds'])),flush=True)
        _finish(folder,manifest,max_pairs)
    return _finish(folder,manifest,max_pairs)


def check_pairs(folder,*,max_pairs=1000000):
    from .algorithms.finite_tissue_boolean import make_solid
    folder,manifest=_read(folder);meshes={};solids={}
    for name in manifest['mesh_files']:
        path=_cache_path(folder,manifest,'cell',name,max_pairs)
        if not path.exists()or not json.loads(path.read_text())['passed']:
            raise ValueError('Every cell must pass before occupied-volume checks')
    for first,second in itertools.combinations(sorted(manifest['mesh_files']),2):
        path=_cache_path(folder,manifest,'pair',[first,second],max_pairs)
        if path.exists():continue
        began=time.monotonic()
        for name in (first,second):
            if name not in meshes:meshes[name]=Mesh.load(manifest['mesh_files'][name]['path'])
        a,b=meshes[first],meshes[second]
        if np.any(a.bounds[0]>b.bounds[1])or np.any(b.bounds[0]>a.bounds[1]):
            row=dict(first=first,second=second,bounds_disjoint=True,passed=True)
        else:
            for name in (first,second):
                if name not in solids:solids[name],_=make_solid(meshes[name].vertices,meshes[name].faces)
            volume=float((solids[first]^solids[second]).volume())
            extent=max(float(np.ptp(np.r_[a.vertices,b.vertices],axis=0).max()),np.finfo(float).tiny)
            bound=2048*np.finfo(float).eps*extent**3*(len(a.faces)+len(b.faces))
            row=dict(first=first,second=second,overlap_volume_m3=volume,numerical_volume_bound_m3=bound,passed=volume<=bound)
        row['elapsed_seconds']=time.monotonic()-began;atomic_json(path,row)
    return _finish(folder,manifest,max_pairs)
