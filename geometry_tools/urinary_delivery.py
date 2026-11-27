"""Encode a continuous urinary arrangement in its existing renal charts.

Each kidney keeps its original Float32 origin and literal retained surfaces.
New cross-chart controls use a common representable lattice. Shared charts
are checked after decoding, including opposite corner-normal words.
"""
from collections import defaultdict
from pathlib import Path
import json

import numpy as np

from .algorithms.finite_tissue_boolean import CommonPositionBank
from .algorithms.constrained_boundary_disk import cancel_opposed_faces
from .algorithms.quantized_surface_topology import conform_stored_surfaces
from .audit import audit
from .bladder_carriers import recover
from .common import atomic_json, atomic_npz, digest_file, identity
from .model import AtlasReader, Mesh
from .tissue import normal_words, pack


def group_for(pid):
    for side in ('left', 'right'):
        if pid.startswith('RENALDETAIL:male:'+side+':'):
            return side
    if pid.startswith('BLADDER:male:'):
        return 'bladder'
    if pid.startswith('URETHRA:male:'):
        return 'urethra'
    if pid=='ZA:Prostate':
        return 'prostate'
    raise ValueError('Unassigned urinary chart: '+pid)


def shared_encoding(points, surfaces, cells, origins, *, chart_for=None):
    """Round new shared controls once, exactly representably in every chart."""
    membership = defaultdict(set)
    for pid, uses in cells.items():
        group = chart_for[pid] if chart_for is not None else group_for(pid)
        for use in uses:
            for i in np.unique(surfaces[use['surface']]):
                membership[int(i)].add(group)
    encoded = np.asarray(points, float).copy(); shared = 0
    for groups in set(tuple(sorted(v)) for v in membership.values()):
        ids = np.array([i for i,v in membership.items() if tuple(sorted(v)) == groups], int)
        os = np.asarray([origins[g] for g in groups]); q = points[ids]
        if len(groups) == 1:
            encoded[ids] = (q-os[0]).astype('<f4').astype(float)+os[0]
        else:
            shared += len(ids)
            extent = np.abs(q[:,None]-os).max(1)
            _, exponent = np.frexp(extent)
            step = np.exp2(np.maximum(exponent-24, -149))
            proposed = os[0]+np.rint((q-os[0])/step)*step
            for origin in os:
                if not np.array_equal((proposed-origin).astype('<f4').astype(float)+origin, proposed):
                    raise ValueError('A shared urinary control is not exactly representable in every chart')
            encoded[ids] = proposed
    return encoded, dict(shared_controls=shared, exact_cross_chart_positions=True,
        maximum_encoding_displacement_m=float(np.linalg.norm(encoded-points, axis=1).max()))


def decoded_interface_ledger(directory,delivery):
    """Replay actual cross-chart positions and words after a normal update."""
    directory=Path(directory);meta={};meshes={};result=[]
    with np.load(delivery['source'],allow_pickle=False)as z:
        positions=z['positions'].copy();surfaces={s:z[k].copy()for s,k in delivery['surfaces'].items()}
    for group,row in delivery['groups'].items():
        if digest_file(directory/group/'cells.bin')!=row['bank_sha256']:
            raise ValueError('A normal-updated chart bank changed')
        meta[group]=json.loads((directory/group/'delivery.json').read_text())
    def corners(pid,sid):
        group=group_for(pid)
        if pid not in meshes:
            m=Mesh.load(directory/group/(pid.replace(':','_')+'.npz'))
            check=meta[group]['checks']['cells'][pid]
            if m.geometry_sha!=check['geometry_sha256']or m.attributes_sha!=check['attributes_sha256']:
                raise ValueError('A decoded chart no longer matches its audited attributes')
            meshes[pid]=m
        m=meshes[pid];ordinal=next(i for i,u in enumerate(delivery['cells'][pid])if u['surface']==sid)
        draw=meta[group]['parts'][pid]['tissueDrawGroups'][ordinal]
        fs=m.faces[draw['start']//3:(draw['start']+draw['count'])//3]
        xyz=m.vertices[fs];words=m.attributes['normals'][fs].astype(int)*m.attributes['normalSigns'][fs,None]
        return {tuple(sorted(tuple(q)for q in tri)):{tuple(q):tuple(map(int,w))for q,w in zip(tri,n)}for tri,n in zip(xyz,words)}
    for pid,uses in delivery['cells'].items():
        for use in uses:
            peer=use.get('neighbor')
            if not peer or pid>=peer or group_for(pid)==group_for(peer):continue
            sid=use['surface'];a,b=corners(pid,sid),corners(peer,sid)
            expected={tuple(sorted(tuple(q)for q in tri))for tri in positions[surfaces[sid]]}
            if set(a)!=set(b)or set(a)!=expected or any(a[t][q]!=tuple(-x for x in b[t][q])for t in a for q in a[t]):
                raise ValueError('A current decoded interface is not literally reciprocal')
            result.append(dict(first=pid,second=peer,surface=sid,triangles=len(a),
                coordinate_normal_signature=identity([(k,sorted(v.items()))for k,v in sorted(a.items())]),
                actual_decoded_positions_and_opposed_words_match=True))
    return result


def stored_source_face_lineage(positions,faces,source_faces,source_tags):
    """Replay real stored triangles against their literal source fragments."""
    positions=np.asarray(positions);source_tags=np.asarray(source_tags,int)
    if np.all(source_tags<0):return np.full(len(faces),-1,int)
    lookup={}
    for f,parent in zip(source_faces,source_tags):
        q=positions[f]
        if any(np.array_equal(q[a],q[b])for a,b in ((0,1),(1,2),(2,0))):continue
        key=tuple(sorted(map(tuple,q)));lookup.setdefault(key,set()).add(int(parent))
    result=[]
    for f in faces:
        parents=lookup.get(tuple(sorted(map(tuple,positions[f]))),set())
        if len(parents)!=1:raise ValueError('A stored source surface lost unambiguous literal face lineage')
        result.append(next(iter(parents)))
    return np.array(result,int)


def build(arrangement_dir, catalogue, output, *, coordinates_dir):
    arrangement_dir, catalogue, output = Path(arrangement_dir).resolve(), Path(catalogue).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    spec = json.loads((arrangement_dir/'arrangement.json').read_text())
    atlas = json.loads(catalogue.read_text()); reader = AtlasReader(); byid = {p['id']:p for p in atlas['parts']}
    marker = atlas['finiteRenalIntegration']
    if marker['schema'] != 'renal-local-layers-v1':
        raise ValueError('Expected the installed finite three-coat renal assembly')
    with np.load(arrangement_dir/'arrangement.npz', allow_pickle=False) as z:
        positions = z['positions'].copy(); bank = CommonPositionBank(positions)
        new_surfaces = {sid:z[key].copy() for sid,key in spec['surfaces'].items()}
        tags = {sid:z[key].copy() for sid,key in spec['source_face_arrays'].items()}
        source_uvs={sid:z[key].copy()for sid,key in spec.get('corner_uvs',{}).items()}
        source_words={sid:z[key].copy()for sid,key in spec.get('corner_word_overrides',{}).items()}
    coordinates_dir = Path(coordinates_dir).resolve()
    if digest_file(coordinates_dir/'coordinates.npz') != spec['coordinates_sha256']:
        raise ValueError('The source-bound urinary coordinates changed')
    coordinates = json.loads((coordinates_dir/'coordinates.json').read_text())
    if digest_file(catalogue) != coordinates['catalogue_sha256']:
        raise ValueError('The installed renal source catalogue changed')
    root_ids = set(map(int, coordinates['source_root_faces']))
    pids = spec['part_ids']; cells = {pid:[] for pid in pids.values()}; surfaces = {}; old_words = {}
    origins = {'bladder': np.round(positions[np.unique(np.concatenate([new_surfaces[u['surface']]
        for name,uses in spec['cells'].items() if name.startswith('bladder-') for u in uses]))].mean(0)*1024)/1024}
    urethral_faces=[new_surfaces[u['surface']] for name,uses in spec['cells'].items()
                    if name.startswith('urethra-') for u in uses]
    if urethral_faces:
        origins['urethra']=np.round(positions[np.unique(np.concatenate(urethral_faces))].mean(0)*1024)/1024
    if 'prostate'in spec['cells']:
        # Native prostate exterior/UVs use absolute Float32 source positions.
        origins['prostate']=np.array(byid['ZA:Prostate'].get('geometryOrigin',[0.,0.,0.]),float)
    protected = {}; root_surfaces = [];uvs={}
    for group in marker['assemblies']:
        side = group['side']; origins[side] = np.array(byid[group['partIds'][0]]['geometryOrigin'], float)
        retained = [pid for pid in group['partIds'] if ':ureter-wall:' not in pid]
        if len(retained) != 30:
            raise ValueError('Expected thirty retained kidney/UPJ tissue cells')
        for pid in retained:
            cells[pid] = []
        for pid in retained:
            record = byid[pid]
            if not np.array_equal(record['geometryOrigin'], origins[side]):
                raise ValueError('Installed kidney cells do not share their declared chart')
            mesh = reader.read_part(atlas, record, frame=spec['frame'])
            mapping = np.array([bank.insert(q) for q in mesh.vertices], int)
            protected[pid] = mesh.surface_sha
            for number, use in enumerate(record['tissueDrawGroups']):
                peer = use.get('neighbor')
                if peer in retained and pid > peer:
                    continue
                if peer and peer not in cells:
                    raise ValueError('Unrepresented adjoining renal tissue: '+peer)
                fs = mesh.faces[use['start']//3:(use['start']+use['count'])//3]
                sid = 'retained:'+pid+':'+str(number)
                surfaces[sid] = mapping[fs]
                words = mesh.attributes['normals'][fs].astype(np.int32)*mesh.attributes['normalSigns'][fs,None]
                old_words[sid] = words.astype('<i2')
                cells[pid].append(dict(surface=sid, side=1, neighbor=peer))
                if peer:
                    cells[peer].append(dict(surface=sid, side=-1, neighbor=pid))
                    if ':ureter-wall:' in peer:
                        root_surfaces.append(sid)
    for name, uses in spec['cells'].items():
        pid = pids[name]
        for use in uses:
            peer_name = use.get('neighbor'); peer = pids[peer_name] if peer_name else None
            if peer and pid > peer:
                continue
            sid = use['surface']; keep = ~np.isin(tags[sid], list(root_ids))
            fs = new_surfaces[sid][keep]
            if not len(fs):
                continue
            if use['side'] < 0:
                fs = fs[:, [0,2,1]]
            new_sid = 'new:'+sid; surfaces[new_sid] = fs
            if sid in source_uvs:
                uv=source_uvs[sid][keep]
                uvs[new_sid]=uv if use['side']>0 else uv[:,[0,2,1]]
            if sid in source_words:
                field=source_words[sid][keep]
                old_words[new_sid]=field if use['side']>0 else -field[:,[0,2,1]]
            cells[pid].append(dict(surface=new_sid, side=1, neighbor=peer))
            if peer:
                cells[peer].append(dict(surface=new_sid, side=-1, neighbor=pid))
    positions = bank.positions
    encoded, encoding = shared_encoding(positions, surfaces, cells, origins)
    surfaces,stored_parents,storage_conformity=conform_stored_surfaces(encoded,surfaces,source_positions=positions)
    old_words={sid:field[stored_parents[sid]]for sid,field in old_words.items()}
    uvs={sid:field[stored_parents[sid]]for sid,field in uvs.items()}
    # Existing tissues are immutable here; only new collar/layer controls may
    # absorb quantization. A changed original chart is an error, not a fit.
    for pid in protected:
        fs = np.concatenate([surfaces[u['surface']] if u['side'] == 1 else surfaces[u['surface']][:,[0,2,1]] for u in cells[pid]])
        actual = Mesh(encoded, fs, spec['frame'])
        if actual.surface_sha != protected[pid]:
            raise ValueError('A retained kidney surface moved during encoding: '+pid)
    zero={}
    for sid,fs in surfaces.items():
        q=encoded[fs];area=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
        ids=np.flatnonzero(np.all(area==0,axis=1))
        if len(ids):zero[sid]=ids.tolist()
    if zero:
        output.mkdir(parents=True)
        witness={'positions':positions,'encoded_positions':encoded};keys={};uvkeys={};wordkeys={}
        for i,(sid,fs)in enumerate(surfaces.items()):
            key='surface_'+str(i);witness[key]=fs;keys[sid]=key
            if sid in uvs:key='uv_'+str(i);witness[key]=uvs[sid];uvkeys[sid]=key
            if sid in old_words:key='normal_'+str(i);witness[key]=old_words[sid];wordkeys[sid]=key
        atomic_npz(output/'encoding-failure.npz',**witness)
        atomic_json(output/'encoding-failure.json',dict(surfaces=keys,cells=cells,origins={k:v.tolist()for k,v in origins.items()},
            corner_uvs=uvkeys,corner_words=wordkeys,frame=spec['frame'],zero_area_triangles=zero,
            zero_area_count=sum(map(len,zero.values())),source_arrangement_sha256=digest_file(arrangement_dir/'arrangement.npz'),
            source_catalogue_sha256=digest_file(catalogue),coordinates_moved_after_encoding=False))
        raise ValueError('Actual encoded tissues contain '+str(sum(map(len,zero.values())))+' zero-area triangles; see encoding-failure.json')
    from .algorithms.axis_interface_retriangulation import retriangulate_axis_interfaces
    mutable={'new:'+sid for sid,field in tags.items()if np.all(field<0)}-set(old_words)-set(uvs)
    surfaces,axis_repairs=retriangulate_axis_interfaces(encoded,surfaces,mutable)
    words = normal_words(encoded, {sid:fs for sid,fs in surfaces.items() if sid not in old_words})
    words.update(old_words)
    output.mkdir(parents=True)
    arrays = {'positions': encoded}; surface_keys = {}; word_keys = {};uv_keys={};source_face_keys={}
    for i,(sid,fs) in enumerate(surfaces.items()):
        key = 'surface_'+str(i); arrays[key] = fs; surface_keys[sid] = key
        key = 'normal_'+str(i); arrays[key] = words[sid]; word_keys[sid] = key
        if sid in uvs:
            key='uv_'+str(i);arrays[key]=uvs[sid];uv_keys[sid]=key
        if sid.startswith('new:'):
            original=sid[4:];key='source_face_'+str(i)
            arrays[key]=stored_source_face_lineage(encoded,fs,new_surfaces[original],tags[original]);source_face_keys[sid]=key
    atomic_npz(output/'source.npz', **arrays)
    groups = {}; decoded = {}; external = []
    for group in origins:
        local_cells = {pid:[dict(u, neighbor=(u.get('neighbor') if u.get('neighbor') and group_for(u['neighbor']) == group else None))
                           for u in uses] for pid,uses in cells.items() if group_for(pid) == group}
        used = {u['surface'] for uses in local_cells.values() for u in uses}
        chart_uvs=None
        if group=='prostate':
            if not used<=set(uvs):
                raise ValueError('The prostate chart is missing its source-bound UV field')
            chart_uvs={sid:uvs[sid]for sid in used}
        result = pack(encoded, {sid:surfaces[sid] for sid in used}, local_cells,
                      {pid:pid for pid in local_cells}, spec['frame'], output/group,
                      origin=origins[group], corner_words={sid:words[sid] for sid in used},corner_uvs=chart_uvs, max_pairs=1000000,
                      construction=dict(source_catalogue_sha256=digest_file(catalogue), source_arrangement_sha256=digest_file(arrangement_dir/'arrangement.npz'),
                                        original_renal_surfaces_literal=True, encoding=encoding))
        if not result['passed']:
            raise ValueError('An encoded urinary tissue chart fails its independent geometry checks: '+group)
        for pid, row in result['parts'].items():
            decoded[pid] = Mesh.load(output/group/(pid.replace(':','_')+'.npz'))
            for actual_use, original_use in zip(row['tissueDrawGroups'], cells[pid]):
                peer = original_use.get('neighbor')
                if peer and group_for(peer) != group:
                    actual_use['neighbor'] = peer
        result['external_interfaces'] = {pid:[u for u in cells[pid] if u.get('neighbor') and group_for(u['neighbor']) != group] for pid in local_cells}
        atomic_json(output/group/'delivery.json', result)
        groups[group] = dict(folder=str(output/group), bank_sha256=result['bank_sha256'], geometry_bytes=result['total_bytes'],
                             cells=list(local_cells), origin=origins[group].tolist())
        print('chart',group,len(local_cells),'cells',result['total_bytes'],'bytes',flush=True)
    # Every shared chart has one coordinate/word signature globally. Verify
    # the actual decoded banks as well, rather than trusting source equality.
    for pid, uses in cells.items():
        record=json.loads((output/group_for(pid)/'delivery.json').read_text())['parts'][pid]
        mesh=decoded[pid]
        for use, draw in zip(uses,record['tissueDrawGroups']):
            peer=use.get('neighbor')
            if not peer or group_for(peer) == group_for(pid) or pid > peer:
                continue
            sid=use['surface']; expected=encoded[surfaces[sid]]
            def corners(who, surface):
                meta=json.loads((output/group_for(who)/'delivery.json').read_text())
                row=meta['parts'][who]; i=next(i for i,u in enumerate(cells[who]) if u['surface'] == surface)
                draw=row['tissueDrawGroups'][i]; m=decoded[who]
                fs=m.faces[draw['start']//3:(draw['start']+draw['count'])//3]
                xyz=m.vertices[fs]
                n=m.attributes['normals'][fs].astype(np.int32)*m.attributes['normalSigns'][fs,None]
                return {tuple(sorted(tuple(q) for q in tri)):{tuple(q):tuple(map(int,w)) for q,w in zip(tri,norm)} for tri,norm in zip(xyz,n)}
            a,b=corners(pid,sid),corners(peer,sid)
            if set(a) != set(b) or set(a) != {tuple(sorted(tuple(q) for q in t)) for t in expected}:
                raise ValueError('Encoded cross-chart interface coordinates differ')
            if any(a[t][q] != tuple(-x for x in b[t][q]) for t in a for q in a[t]):
                raise ValueError('Encoded cross-chart interface normals are not reciprocal')
            external.append(dict(first=pid,second=peer,surface=sid,triangles=len(a),
                coordinate_normal_signature=identity([(k,sorted(v.items())) for k,v in sorted(a.items())]),
                actual_decoded_positions_and_opposed_words_match=True))
    result=dict(schema='urinary-local-charts/v1', groups=groups, cells=cells, surfaces=surface_keys, corner_words=word_keys,
        part_ids={pid:pid for pid in cells}, frame=spec['frame'], source=str(output/'source.npz'),
        corner_uvs=uv_keys,source_exterior_normal_overrides=["new:"+sid for sid in source_words],
        source_face_arrays=source_face_keys,
        catalogue=str(catalogue), catalogue_sha256=digest_file(catalogue), external_interfaces=external,
        protected_original_cells=len(protected), original_renal_surfaces_literal=True,
        feature_plan=spec.get('feature_plan'),feature_plan_sha256=spec.get('feature_plan_sha256'),
        opening_curves=spec.get('opening_curves'),opening_surround_width_m=spec.get('opening_surround_width_m'),
        opening_support_method=spec.get('opening_support_method'),opening_windows=spec.get('opening_windows'),
        root_triangles_preserved=spec['root_faces_preserved'], encoding=encoding,
        exact_storage_conformity=storage_conformity,
        axis_interface_retriangulation=axis_repairs,
        total_geometry_bytes=sum(q['geometry_bytes'] for q in groups.values()),
        passed=True, installed=False, anatomical_acceptance=False)
    atomic_json(output/'delivery.json',result)
    return result


def verify_whole(delivery_dir, arrangement_dir, coordinates_dir, output):
    """Independently cancel decoded tissue interfaces and replay the fluid."""
    delivery_dir, arrangement_dir, coordinates_dir, output = map(lambda p:Path(p).resolve(),
        (delivery_dir,arrangement_dir,coordinates_dir,output))
    if output.exists():
        raise FileExistsError(output)
    delivery=json.loads((delivery_dir/'delivery.json').read_text())
    spec=json.loads((arrangement_dir/'arrangement.json').read_text())
    if delivery['passed'] is not True:
        raise ValueError('The encoded urinary charts did not pass')
    meshes=[]; volume=0.
    for group,row in delivery['groups'].items():
        if digest_file(delivery_dir/group/'cells.bin') != row['bank_sha256']:
            raise ValueError('An encoded urinary bank changed')
        meta=json.loads((delivery_dir/group/'delivery.json').read_text())
        for pid in row['cells']:
            if group not in ('bladder','urethra') and ':ureter-wall:' not in pid:
                continue
            mesh=Mesh.load(delivery_dir/group/(pid.replace(':','_')+'.npz'))
            used,inv=np.unique(mesh.faces,return_inverse=True)
            meshes.append(Mesh(mesh.vertices[used],inv.reshape(-1,3),mesh.frame))
            volume+=meta['checks']['cell_volumes_m3'][pid]
    points=np.concatenate([m.vertices for m in meshes]); faces=[]; offset=0
    for m in meshes:
        faces.append(m.faces+offset); offset+=len(m.vertices)
    v,inv=np.unique(points,axis=0,return_inverse=True)
    f,cancellation=cancel_opposed_faces(inv[np.concatenate(faces)])
    whole=Mesh(v,f,delivery['frame'])
    check=audit(whole,require_closed=True,expected_components=1,expected_genus=[2],max_pairs=1000000)
    output.mkdir(parents=True);whole.save(output/'whole-wall.npz')
    if not check['passed']:
        atomic_json(output/'whole-check-failure.json',check)
        raise ValueError('The complete decoded wall fails its independent geometry check')
    inner_keys=set();collapsed_inner_fragments=0
    with np.load(coordinates_dir/'coordinates.npz',allow_pickle=False)as z:
        original_inner=set(map(int,z['source_inner_face_ids']))
    with np.load(delivery['source'],allow_pickle=False)as z, np.load(arrangement_dir/'arrangement.npz',allow_pickle=False)as az:
        source=z['positions']
        for sid,key in spec['surfaces'].items():
            mask=np.isin(az[spec['source_face_arrays'][sid]],list(original_inner))
            if not mask.any():
                continue
            fs=az[key][mask]
            for t in fs:
                q=source[t]
                if any(np.array_equal(q[a],q[b])for a,b in ((0,1),(1,2),(2,0))):
                    collapsed_inner_fragments+=1;continue
                inner_keys.add(tuple(sorted(map(tuple,q))))
    indices=np.array([i for i,t in enumerate(f)if tuple(sorted(tuple(q)for q in v[t]))in inner_keys],int)
    if len(indices) != len(inner_keys):
        raise ValueError('Some actual lining faces disappeared from the decoded whole wall')
    outlet='external-urethral-orifice' if 'urethra' in delivery['groups'] else 'urethral-neck'
    carriers=recover(whole,output/'carriers',inner_face_ids=indices,
        mouth_labels=[outlet,'right-renal-root','left-renal-root'],cap_forbid_source_diagonals=True)
    integral=check['topology']['oriented_triangle_integral_m3'];residual=abs(volume-integral)
    extent=float(np.ptp(v,axis=0).max());bound=2048*np.finfo(float).eps*extent**3*len(f)
    if residual > bound:
        raise ValueError('Decoded tissue volumes do not sum to the one actual wall volume')
    result=dict(schema='decoded-continuous-urinary-replay/v1',passed=True,whole_check=check,
        carriers=carriers,internal_faces_canceled=cancellation,individual_tissue_volume_sum_m3=volume,
        whole_wall_volume_m3=integral,volume_replay_residual_m3=residual,numerical_volume_bound_m3=bound,
        delivery_sha256=digest_file(delivery_dir/'delivery.json'),
        all_actual_inner_faces_retained=True,one_connected_urinary_fluid_domain=True,
        storage_collapsed_inner_fragments=collapsed_inner_fragments,
        computational_disks_for_renderer=False,installed=False,anatomical_acceptance=False)
    atomic_json(output/'verification.json',result)
    return result


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('arrangement_dir');p.add_argument('--catalogue',required=True);p.add_argument('--output',required=True);p.add_argument('--coordinates-dir',required=True)
    args=p.parse_args();r=build(args.arrangement_dir,args.catalogue,args.output,coordinates_dir=args.coordinates_dir)
    print(json.dumps(dict(passed=r['passed'],total_geometry_bytes=r['total_geometry_bytes'],external_interfaces=len(r['external_interfaces'])),indent=2))
