"""Repair shading across urinary tissue selections without changing geometry.

Only normals/render aliases on the newly built physical wall are changed.
Every oriented geometric triangle and protected kidney corner word remains
literal. Numerical geometry proofs are reused only after those indexed
identities are replayed, and actual rendered review remains required.
"""
from pathlib import Path
import copy
import json

import numpy as np

from .algorithms.surface_normals import grouped_normal_words,source_normal_words
from .algorithms.finite_cell_bank import pack_finite_cells
from .algorithms.tube_rim_bridge import end_annulus
from .common import atomic_bytes, atomic_json, atomic_npz, digest_file
from .model import Mesh
from .tissue import decode_packed


def build(delivery_dir, arrangement_dir, coordinates_dir, output,*,physical_source_normals=False):
    delivery_dir, arrangement_dir, coordinates_dir, output=map(lambda p:Path(p).resolve(),
        (delivery_dir,arrangement_dir,coordinates_dir,output))
    if output.exists():
        raise FileExistsError(output)
    delivery=json.loads((delivery_dir/'delivery.json').read_text())
    if delivery.get('passed') is not True:
        raise ValueError('Use independently passed encoded urinary charts')
    spec=json.loads((arrangement_dir/'arrangement.json').read_text())
    coordinates=json.loads((coordinates_dir/'coordinates.json').read_text())
    with np.load(coordinates_dir/'coordinates.npz',allow_pickle=False)as z:
        inner_ids=z['source_inner_face_ids'].copy()
    wall=Mesh.load(coordinates['wall'])
    with np.load(Path(coordinates['wall']).parent/'carriers.npz',allow_pickle=False)as z:
        key='urethral_neck_ring' if 'urethral_neck_ring' in z else 'external_urethral_orifice_ring'
        neck_faces,_=end_annulus(wall.faces,wall.faces[inner_ids],z[key])
    with np.load(delivery['source'],allow_pickle=False)as z:
        arrays={k:z[k].copy()for k in z.files}
    with np.load(arrangement_dir/'arrangement.npz',allow_pickle=False)as z:
        original_tags={sid:z[key].copy()for sid,key in spec['source_face_arrays'].items()}
    sources=delivery['surfaces'];cells=delivery['cells'];surface_keys=delivery['corner_words']
    free={};groups={};source_tags={}
    protected_surfaces=set(delivery.get('source_exterior_normal_overrides',[]))
    for pid,uses in cells.items():
        for use in uses:
            sid=use['surface']
            if pid=='ZA:Prostate' or sid in protected_surfaces or \
                    use.get('neighbor')not in (None,'ZA:Prostate') or not sid.startswith('new:'):
                continue
            original=sid[4:]
            # Root facets are retained from their original UPJ tissue charts.
            root_ids=set(map(int,coordinates['source_root_faces']))
            if sid in delivery.get('source_face_arrays',{}):tags=arrays[delivery['source_face_arrays'][sid]]
            else:
                tags=original_tags[original];tags=tags[~np.isin(tags,list(root_ids))]
            faces=arrays[sources[sid]]
            if len(tags)!=len(faces):
                raise ValueError('The physical surface lost its face lineage')
            free[sid]=faces
            source_tags[sid]=tags
            groups[sid]=np.where(np.isin(tags,inner_ids),'lining',
                np.where(np.isin(tags,neck_faces),'neck-rim','outer'))
    if physical_source_normals:
        labels=np.where(np.isin(np.arange(len(wall.faces)),inner_ids),'lining',
            np.where(np.isin(np.arange(len(wall.faces)),neck_faces),'neck-rim','outer'))
        base,_=grouped_normal_words(wall.vertices,{'wall':wall.faces},{'wall':labels})
        fields,evidence=source_normal_words(arrays['positions'],free,source_tags,wall.vertices,wall.faces,base['wall'])
    else:
        fields,evidence=grouped_normal_words(arrays['positions'],free,groups)
    output.mkdir(parents=True)
    for sid,field in fields.items():
        arrays[surface_keys[sid]]=field
    atomic_npz(output/'source.npz',**arrays)
    result=copy.deepcopy(delivery);protected=0;changed=0
    for group,row in result['groups'].items():
        previous=delivery_dir/group
        if digest_file(previous/'cells.bin')!=row['bank_sha256']:
            raise ValueError('A source encoded bank changed')
        meta=json.loads((previous/'delivery.json').read_text())
        local_cells=meta['cells']
        used={u['surface']for uses in local_cells.values()for u in uses}
        uvs=None
        if group=='prostate':
            uvs={sid:arrays[delivery['corner_uvs'][sid]]for sid in used}
        packed=pack_finite_cells(arrays['positions'],{sid:arrays[sources[sid]]for sid in used},
            {sid:arrays[surface_keys[sid]]for sid in used},local_cells,{pid:pid for pid in local_cells},row['origin'],corner_uvs=uvs)
        meshes=decode_packed(packed,delivery['frame'],row['origin'])
        folder=output/group;folder.mkdir()
        atomic_bytes(folder/'cells.bin',packed.data)
        for pid in row['cells']:
            mesh=meshes[pid];fs=mesh.faces;used_controls,inv=np.unique(fs,return_inverse=True)
            previous_mesh=Mesh.load(previous/(pid.replace(':','_')+'.npz'))
            previous_words=previous_mesh.attributes['normals'][previous_mesh.faces].astype(np.int32)*previous_mesh.attributes['normalSigns'][previous_mesh.faces,None]
            actual_words=mesh.attributes['normals'][fs].astype(np.int32)*mesh.attributes['normalSigns'][fs,None]
            if ':ureter-wall:'not in pid and pid.startswith('RENALDETAIL:'):
                if not np.array_equal(actual_words,previous_words):
                    raise ValueError('A protected kidney tissue normal changed')
                protected+=1
            changed+=int(np.any(actual_words!=previous_words,axis=2).sum())
            mesh=Mesh(mesh.vertices[used_controls],inv.reshape(-1,3),delivery['frame'],
                {k:v[used_controls]for k,v in mesh.attributes.items()})
            if mesh.geometry_sha!=meta['checks']['cells'][pid]['geometry_sha256']:
                raise ValueError('A normal-only correction changed the audited geometry identity')
            meta['checks']['cells'][pid]['attributes_sha256']=mesh.attributes_sha
            mesh.save(folder/(pid.replace(':','_')+'.npz'))
        for pid,part in packed.cells.items():
            for draw,source in zip(part['tissueDrawGroups'],cells[pid]):
                if source.get('neighbor'):
                    draw['neighbor']=source['neighbor']
        meta['parts']=packed.cells
        meta['position_bytes']=int(packed.positions.nbytes);meta['normal_bytes']=int(packed.normal_words.nbytes);meta['total_bytes']=len(packed.data)
        row['geometry_bytes']=len(packed.data)
        row['folder']=str(folder);row['bank_sha256']=digest_file(folder/'cells.bin')
        meta['bank']=str(folder/'cells.bin');meta['bank_sha256']=row['bank_sha256']
        meta['source_geometry_audits']=dict(path=str(previous/'tissue.json'),sha256=digest_file(previous/'tissue.json'),
            every_actual_oriented_triangle_literal=True,actual_geometry_identities_replayed=True)
        meta['normal_field']=evidence
        atomic_json(folder/'delivery.json',meta)
        print('normal chart',group,'geometry literal',flush=True)
    result['source']=str(output/'source.npz')
    from .urinary_delivery import decoded_interface_ledger
    result['external_interfaces']=decoded_interface_ledger(output,result)
    result['normal_field']=dict(**evidence,changed_corner_words=changed,protected_original_tissue_normal_fields=protected,
        protected_source_normal_surfaces=sorted(protected_surfaces),source_uv_fields_retained=True,
        every_actual_oriented_triangle_literal=True,previous_delivery_sha256=digest_file(delivery_dir/'delivery.json'))
    result['total_geometry_bytes']=sum(r['geometry_bytes']for r in result['groups'].values())
    verified=delivery_dir/'whole-verification/verification.json'
    if verified.exists():
        proof=json.loads(verified.read_text())
        if not proof['passed']or proof['delivery_sha256']!=digest_file(delivery_dir/'delivery.json'):
            raise ValueError('The previous whole geometry verification does not bind this source')
        result['whole_geometry_verification']=dict(path=str(verified),sha256=digest_file(verified),
            actual_indexed_geometry_unchanged=True,only_render_normals_and_aliases_changed=True)
    result['anatomical_acceptance']=False;result['installed']=False
    atomic_json(output/'delivery.json',result)
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ['delivery-dir','arrangement-dir','coordinates-dir','output']:
        p.add_argument('--'+flag,required=True)
    r=build(**vars(p.parse_args()))
    print(json.dumps(r['normal_field'],indent=2))
