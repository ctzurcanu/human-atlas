"""Finite cells sharing one interface/attribute bank, checked after encoding."""
from pathlib import Path
from collections import Counter
import json,itertools,time
import numpy as np
from .common import atomic_json,atomic_bytes,atomic_npz,digest_file
from .model import Mesh
from .audit import audit,topology
from .contact import pair_contacts


def layer_arrangement(fields,faces,names):
    from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim,finite_wall_cells
    fields=np.asarray(fields,float);faces=np.asarray(faces)
    if fields.ndim!=3 or fields.shape[2]!=3 or len(fields)<2 or len(names)!=len(fields)-1 or len(set(names))!=len(names):raise ValueError('K corresponding interfaces require K-1 unique tissue names')
    if not np.isfinite(fields).all():raise ValueError('Finite explicitly constructed interfaces required')
    n=fields.shape[1];surfaces={};cells={name:[]for name in names}
    for k in range(len(fields)):
        sid='interface:'+str(k);surfaces[sid]=faces+k*n
        if k<len(names):cells[names[k]].append(dict(surface=sid,side=1,neighbor=names[k-1]if k else None))
        if k>0:cells[names[k-1]].append(dict(surface=sid,side=-1,neighbor=names[k]if k<len(names)else None))
    # Close wall thickness only at the actual open carrier boundaries.
    # No face spans a physiological lumen or a regional partition.
    mouths=boundary_loops(faces)
    for k,name in enumerate(names):
        for j,mouth in enumerate(mouths):
            sid=f'rim:{k}:{j}';surfaces[sid]=annular_rim(mouth+k*n,mouth+(k+1)*n);cells[name].append(dict(surface=sid,side=1,neighbor=None))
    return fields.reshape(-1,3),surfaces,cells,finite_wall_cells(fields,faces)


def normal_words(positions,surfaces):
    result={}
    for sid,faces in surfaces.items():
        # Only this chart contributes to its normal field. Avoid allocating
        # a whole-organ accumulator again for every small internal surface.
        used,local_faces=np.unique(faces,return_inverse=True);local_faces=local_faces.reshape(-1,3)
        q=positions[used][local_faces];normals=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);length=np.linalg.norm(normals,axis=1)
        if np.any(length==0):raise ValueError('Cannot encode a zero-area tissue surface')
        # Area-weighted within each explicit surface; hard tissue interfaces
        # keep distinct charts and reciprocal ownership uses opposite signs.
        sums=np.zeros((len(used),3),float)
        for k in range(3):np.add.at(sums,local_faces[:,k],normals)
        corner=sums[local_faces];n=np.linalg.norm(corner,axis=2);corner=np.divide(corner,n[:,:,None],out=np.zeros_like(corner),where=n[:,:,None]>0)
        bad=np.einsum('tki,ti->tk',corner,normals)<=0
        corner[bad]=np.broadcast_to((normals/length[:,None])[:,None],corner.shape)[bad]
        result[sid]=np.rint(np.clip(corner,-1,1)*32767).astype('<i2')
    return result


def decode_packed(packed,frame,origin):
    bank=packed.data;meshes={}
    for cell,row in packed.cells.items():
        n=row['vertexCount'];bank_positions=np.ndarray((n,3),'<f4',buffer=bank,offset=row['positions'])
        f=np.ndarray((row['indexCount']//3,3),'<u4',buffer=bank,offset=row['indices'])
        used,local_faces=np.unique(f,return_inverse=True);f=local_faces.reshape(-1,3)
        v=bank_positions[used].astype(float)+origin
        normal=np.ndarray((n,3),'<i2',buffer=bank,offset=row['normals'])[used];sign=np.ndarray((n,),'i1',buffer=bank,offset=row['normalSigns'])[used]
        attributes={'normals':normal,'normalSigns':sign}
        if 'uvs'in row:attributes['uvs']=np.ndarray((n,2),'<f4',buffer=bank,offset=row['uvs'])[used]
        meshes[cell]=Mesh(v,f,frame,attributes)
    return meshes


def validate_cells(meshes,*,max_pairs=250000,max_tests=2000000,checkpoint_directory=None):
    from geometry_tools.algorithms.finite_tissue_boolean import make_solid
    # A bounded job may stop between expensive cells or Boolean pairs. Reuse
    # only exact-geometry results from this implementation and native version.
    import importlib.metadata
    checkpoint=Path(checkpoint_directory).resolve() if checkpoint_directory else None
    signature=None
    if checkpoint is not None:
        checkpoint.mkdir(parents=True,exist_ok=True)
        from .common import identity
        signature=identity('finite-tissue-validation-checkpoint/v4',np.__version__,importlib.metadata.version('manifold3d'),importlib.metadata.version('libigl'),
            [(name,digest_file(Path(__file__).with_name(name))) for name in ['tissue.py','audit.py','contact.py','distance.py','spatial.py','model.py','geometry_instance.py']],
            digest_file(Path(__file__).with_name('algorithms')/'finite_tissue_boolean.py'),digest_file(Path(__file__).with_name('algorithms')/'shared_cell_disjoint.py'),max_pairs,max_tests)
    def checked(kind,geometry,run):
        if checkpoint is None:return run()
        from .common import identity
        key=identity(signature,kind,geometry);path=checkpoint/(key+'.json')
        if path.exists():
            record=json.loads(path.read_text())
            if record['signature']!=signature or record['geometry']!=geometry or record['kind']!=kind or record.get('result_sha256')!=identity(record['result']):raise ValueError('Mismatched tissue-validation checkpoint')
            return record['result']
        start=time.monotonic();result=run()
        atomic_json(path,dict(signature=signature,geometry=geometry,kind=kind,result=result,result_sha256=identity(result),elapsed_seconds=time.monotonic()-start))
        return result
    from .geometry_instance import instance_preimage
    instances={name:instance_preimage(mesh)for name,mesh in meshes.items()}
    poses=[value[1]for value in instances.values()if value is not None]
    common_instance=len(poses)==len(meshes)and len({json.dumps({k:p[k]for k in ['rotation','scale','origin']},sort_keys=True)for p in poses})==1
    boolean_meshes={name:value[0]for name,value in instances.items()}if common_instance else meshes
    volume_scale=poses[0]['scale']**3 if common_instance else 1.
    checks={};solids={};volumes={};failures=[]
    for i,(name,m) in enumerate(meshes.items()):
        checks[name]=checked('cell-audit',[m.geometry_sha,m.frame],lambda:audit(m,require_closed=True,max_pairs=max_pairs,max_tests=max_tests))
        if checkpoint is not None and ((i+1)%10==0 or i+1==len(meshes)):
            print(f'[tissue] cell audits {i+1}/{len(meshes)} saved',flush=True)
    for name,m in meshes.items():
        if not checks[name]['passed']:failures.append('cell:'+name);continue
        try:
            native=boolean_meshes[name];solid,_=make_solid(native.vertices,native.faces);solids[name]=solid;volumes[name]=float(solid.volume())*volume_scale
            if volumes[name]<=0:failures.append('nonpositive_cell:'+name)
        except ValueError as e:failures.append('solid:'+name+':'+str(e))
    overlap=[]
    for a,b in itertools.combinations(solids,2):
        ma,mb=meshes[a],meshes[b]
        if np.any(ma.bounds[0]>mb.bounds[1])or np.any(mb.bounds[0]>ma.bounds[1]):continue
        def pair_volume():
            ca,cb=boolean_meshes[a],boolean_meshes[b];amount=float((solids[a]^solids[b]).volume())*volume_scale;extent=max(float(np.ptp(np.r_[ca.vertices,cb.vertices],axis=0).max()),np.finfo(float).tiny)
            roundoff=2048*np.finfo(float).eps*extent**3*(len(ma.faces)+len(mb.faces))*volume_scale
            return dict(amount=amount,roundoff=roundoff)
        result=checked('occupied-volume',[ma.geometry_sha,ma.frame,mb.geometry_sha,mb.frame],pair_volume)
        amount,roundoff=result['amount'],result['roundoff'];ok=amount<=roundoff
        certificate=None
        if not ok:
            from .algorithms.shared_cell_disjoint import shared_cell_disjoint
            certificate=checked('shared-cell-disjoint',[ma.geometry_sha,ma.frame,mb.geometry_sha,mb.frame],lambda:shared_cell_disjoint(boolean_meshes[a],boolean_meshes[b],max_pairs=max_pairs,max_tests=max_tests))
            ok=certificate['passed']
        overlap.append(dict(first=a,second=b,overlap_volume_m3=amount,numerical_volume_bound_m3=roundoff,passed=ok,
                            **({'independent_disjointness_certificate':certificate}if certificate is not None else{})))
        if not ok:failures.append('volume_overlap:'+a+':'+b)
    return dict(passed=not failures,failures=failures,cells=checks,cell_volumes_m3=volumes,pair_occupied_volume_checks=overlap,
                all_emitted_cells_checked=True,anatomical_acceptance=False)


def pack(positions,surfaces,cells,part_ids,frame,output,*,origin=(0.,0.,0.),corner_words=None,corner_uvs=None,max_pairs=250000,max_tests=2000000,construction=None,checkpoint_directory=None):
    from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
    output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    origin=np.asarray(origin,float);local=np.asarray(positions,float)-origin;encoded=local.astype('<f4').astype(float)
    # Construct/check normals against the actual emitted Float32 faces.
    words=corner_words or normal_words(encoded,surfaces)
    packed=pack_finite_cells(positions,surfaces,words,cells,part_ids,origin,corner_uvs=corner_uvs)
    meshes=decode_packed(packed,frame,origin);checks=validate_cells(meshes,max_pairs=max_pairs,max_tests=max_tests,checkpoint_directory=checkpoint_directory)
    if construction and construction.get('finite_wall_cells_pass')is False:checks['passed']=False;checks['failures'].append('finite_wall_cells')
    output.mkdir(parents=True);bank=output/'cells.bin';atomic_bytes(bank,packed.data)
    arrays={'positions':np.asarray(positions,float)};manifest_surfaces={}
    for i,(sid,faces)in enumerate(surfaces.items()):key='surface_'+str(i);arrays[key]=np.asarray(faces,int);manifest_surfaces[sid]=key
    atomic_npz(output/'arrangement.npz',**arrays)
    for name,m in meshes.items():m.save(output/(part_ids[name].replace(':','_').replace('/','_')+'.npz'))
    manifest=dict(schema='shared-finite-tissues/v1',frame=frame,bank=str(bank),bank_sha256=digest_file(bank),arrangement=str(output/'arrangement.npz'),surfaces=manifest_surfaces,
                  cells=cells,part_ids=part_ids,parts=packed.cells,position_bytes=int(packed.positions.nbytes),normal_bytes=int(packed.normal_words.nbytes),
                  total_bytes=len(packed.data),positions_encoded_once=True,reciprocal_normals_use_same_words_with_opposite_sign=True,
                  maximum_position_encoding_error_m=float(np.linalg.norm(encoded+origin-np.asarray(positions),axis=1).max()),construction=construction,checks=checks,
                  passed=checks['passed'],atlas_modified=False,anatomical_acceptance=False)
    atomic_json(output/'tissue.json',manifest);return manifest


def build_layers(source,output,names,frame=None,origin=(0.,0.,0.),max_pairs=250000,max_tests=2000000):
    with np.load(source,allow_pickle=False)as z:
        actual_frame=frame or(str(z['frame'])if 'frame'in z else None)
        if not actual_frame:raise ValueError('An explicit coordinate frame is required')
        fields=z['fields'].astype(float)+(z['geometryOrigin']if 'geometryOrigin'in z else 0);faces=z['faces'].copy()
    positions,surfaces,cells,construction=layer_arrangement(fields,faces,names)
    return pack(positions,surfaces,cells,{n:n for n in names},actual_frame,output,origin=origin,max_pairs=max_pairs,max_tests=max_tests,construction=construction)


def add_parser(sub):
    p=sub.add_parser('tissue');s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('layers');q.add_argument('source');q.add_argument('--names',nargs='+',required=True);q.add_argument('--frame');q.add_argument('--origin',nargs=3,type=float,default=[0,0,0]);q.add_argument('--max-pairs',type=int,default=250000);q.add_argument('--max-tests',type=int,default=2000000);q.add_argument('--output',required=True)
    q=s.add_parser('pack');q.add_argument('spec');q.add_argument('--output',required=True)
    p.set_defaults(handler=run)


def run(args,db):
    if args.action=='layers':return build_layers(args.source,args.output,args.names,args.frame,args.origin,args.max_pairs,args.max_tests)
    spec=json.loads(Path(args.spec).read_text())
    with np.load(spec['source'],allow_pickle=False)as z:
        positions=z['positions'].copy();surfaces={sid:z[key].copy()for sid,key in spec['surfaces'].items()}
        words={sid:z[key].copy()for sid,key in spec['corner_words'].items()}if spec.get('corner_words')else None
        uvs={sid:z[key].copy()for sid,key in spec['corner_uvs'].items()}if spec.get('corner_uvs')else None
    return pack(positions,surfaces,spec['cells'],spec['part_ids'],spec['frame'],args.output,origin=spec.get('origin',[0,0,0]),corner_words=words,corner_uvs=uvs,max_pairs=spec.get('max_pairs',250000),max_tests=spec.get('max_tests',2000000))
