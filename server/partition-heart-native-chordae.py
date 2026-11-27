"""Finite source-preserving sail/strand cells with shared solid neck sections."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,atomic_bytes,digest_file
from geometry_tools.geometry_instance import instance_pose
from geometry_tools.distance import barycentric
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank,split_closed_cell
from geometry_tools.algorithms.shared_plane_refinement import safe_plane_coordinate
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.tissue import decode_packed,validate_cells
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--parts',nargs='+',type=int,default=[1,5,8,9,10]);args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';src=h/'apparatus-source-v203';out=h/args.stage;out.mkdir(exist_ok=False);reader=AtlasReader();rows=[];binary=bytearray();catalogue=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(catalogue.read_text());assert digest_file(catalogue)==json.loads((src/'source.json').read_text())['catalogue_sha256']
for index in args.parts:
    source=Mesh.load(src/f'part-{index}.npz');assert json.loads((src/f'part-{index}-full-audit.json').read_text())['passed'];record=source.source['record'];original=reader.reload(source.source,source.frame);assert original.geometry_sha==source.geometry_sha
    old_bank=reader._bank(Path(source.source['bank_path']));used=np.array(source.source['source_vertex_ids'],int);raw=reader.view(old_bank,record['positions'],record['vertexCount'],'<f4',3)[used].astype(float);p,inv=np.unique(raw,axis=0,return_inverse=True);f=inv[source.faces];old_words=reader.view(old_bank,record['normals'],record['vertexCount'],'<i2',3)[used][source.faces].copy()
    uv=source.attributes.get('uvs');old_uv=uv[source.faces]if uv is not None else None
    if 'normalSigns'in source.attributes:assert np.all(source.attributes['normalSigns']==1)
    tri=p[f];nn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);bad=np.einsum('tki,ti->tk',old_words.astype(float),nn)<=0;repaired=old_words.copy();face_words=np.rint(nn/np.linalg.norm(nn,axis=1)[:,None]*32767).astype('<i2');repaired[bad]=np.broadcast_to(face_words[:,None],repaired.shape)[bad]
    profile=json.loads((h/'chordal-sections-v207'/f'part-{index}.json').read_text());candidate=profile['measured_neck_candidate'];assert candidate
    world_n=np.array(profile['normal']);world_point=np.array(profile['origin_m'])+world_n*candidate['coordinate_m'];pose=instance_pose(record);rotation,scale=pose if pose is not None else(np.eye(3),1.);geoscale=float(record.get('geometryScale',1.));bias=np.array(record.get('geometryOrigin',[0.,0.,0.]))+(np.array(record.get('sourceOffset',[0.,0.,0.]))@rotation.T)*scale;n=world_n@rotation;desired=float((world_point-bias)@world_n/(scale*geoscale));value=p@n
    coordinate=safe_plane_coordinate(value,desired,search_radius_m=.00008/(scale*geoscale),minimum_separation_m=1e-9/(scale*geoscale));plane_evidence=dict(desired_coordinate_raw=desired,chosen_coordinate_raw=coordinate,minimum_source_vertex_distance_m=float(np.min(abs(value-coordinate))*scale*geoscale),search_radius_m=.00008)
    bank=SharedCutBank(p);split=split_closed_cell(bank,f,n,coordinate,f'cusp-{index}-native-neck',tolerance=1e-12);q=np.array(bank.positions)
    cap=split.interface_faces;count=len(cap);surfaces={'chord-exterior':split.negative_faces[:-count],'leaflet-exterior':split.positive_faces[:-count],'neck':cap};parents={'chord-exterior':split.negative_ancestry[:-count],'leaflet-exterior':split.positive_ancestry[:-count]};words={};uvs={}if old_uv is not None else None
    for sid,pa in parents.items():
        faces=surfaces[sid];controls=p[f[pa]];field=np.empty((len(faces),3,3),dtype='<i2');texture=np.empty((len(faces),3,2),dtype='<f4')if old_uv is not None else None
        for corner in range(3):
            weights,valid=barycentric(q[faces[:,corner]],controls);assert valid.all();v=np.einsum('tk,tkj->tj',weights,repaired[pa].astype(float));v/=np.linalg.norm(v,axis=1)[:,None];field[:,corner]=np.rint(v*32767).astype('<i2')
            if texture is not None:texture[:,corner]=np.einsum('tk,tkj->tj',weights,old_uv[pa]).astype('<f4')
            for k in range(3):
                literal=faces[:,corner]==f[pa,k];field[literal,corner]=repaired[pa[literal],k]
                if texture is not None:texture[literal,corner]=old_uv[pa[literal],k]
        words[sid]=field
        if uvs is not None:uvs[sid]=texture
    words['neck']=np.broadcast_to(np.rint(n*32767).astype('<i2'),(len(cap),3,3)).copy()
    if uvs is not None:
        axis=np.eye(3)[np.argmin(abs(n))];u=np.cross(n,axis);u/=np.linalg.norm(u);v=np.cross(n,u);uvs['neck']=((q[cap]-q[cap].mean(axis=(0,1)))@np.c_[u,v]*100).astype('<f4')
    source_id=record['id'];chord_id=f'TA98HEART:male:chordae:source-{index}';ids={'leaflet':source_id,'chordae':chord_id};cells={'leaflet':[dict(surface='leaflet-exterior',side=1),dict(surface='neck',side=-1,neighbor='chordae')],'chordae':[dict(surface='chord-exterior',side=1),dict(surface='neck',side=1,neighbor='leaflet')]};packed=pack_finite_cells(q,surfaces,words,cells,ids,corner_uvs=uvs)
    local=out/f'part-{index}';local.mkdir();atomic_bytes(local/'cells.bin',packed.data);atomic_npz(local/'arrangement.npz',positions=q,source_positions=p,source_faces=f,source_vertex_map=inv,chord_parent_faces=parents['chord-exterior'],leaflet_parent_faces=parents['leaflet-exterior'],**{k.replace('-','_'):v for k,v in surfaces.items()});decoded=decode_packed(packed,source.frame,np.zeros(3));checks=validate_cells(decoded,max_pairs=1000000,max_tests=7000000,checkpoint_directory=local/'checkpoints');atomic_json(local/'checked-geometry.json',checks)
    while len(binary)%4:binary.append(0)
    offset=len(binary);binary.extend(packed.data);layout=copy.deepcopy(packed.cells)
    for part in layout.values():
        for key in ['positions','normals','normalSigns','uvs','indices']:
            if key in part:part[key]+=offset
        part.update({k:copy.deepcopy(record[k])for k in ['geometryOrigin','geometryScale','geometryInstance','sourceOffset']if k in record})
    row=dict(source=index,source_part_id=source_id,source_geometry_sha256=source.geometry_sha,parts=layout,part_ids=ids,cells=cells,source_record=record,plane_normal_raw=n.tolist(),plane_coordinate_raw=coordinate,plane_evidence=plane_evidence,split=split.evidence,cut_loops=len(split.boundary_loops),source_normal_corners_repaired=int(bad.sum()),source_uvs_preserved=uv is not None,passed=checks['passed'],failures=checks['failures'],new_bytes=len(packed.data));rows.append(row);atomic_json(local/'construction.json',row);print(json.dumps({k:v for k,v in row.items()if k in ['source','cut_loops','source_normal_corners_repaired','passed','failures','new_bytes']}),flush=True)
atomic_bytes(out/'cells.bin',bytes(binary));atomic_json(out/'construction.json',dict(passed=all(r['passed']for r in rows),parts=rows,bank_sha256=digest_file(out/'cells.bin'),bytes=len(binary),new_images=0,source_catalogue_sha256=digest_file(catalogue),elapsed_seconds=time.monotonic()-start,installed=False))
