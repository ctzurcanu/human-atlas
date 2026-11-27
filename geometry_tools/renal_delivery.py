"""Focused local delivery of finite collecting coats, using the fixed renal donor.

The original physical regional annuli constrain the construction tetrahedra.
Only the tissue is partitioned. No computational urinary disk is emitted.
"""
from collections import defaultdict
from pathlib import Path
import json,gzip,time
import numpy as np
import point_cloud_utils as pcu
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .renal_master import bound,OUT,FRAME
from .common import ROOT,atomic_json,atomic_npz,digest_file
from .harmonic_wall import tetrahedralize,layers
from .model import Mesh
from .audit import audit,topology
from .tissue import normal_words,decode_packed
from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank,make_solid
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces,contract_quantized_intersection_edges
from geometry_tools.algorithms.float32_diagonal_repair import repair_diagonals
from geometry_tools.algorithms.encoded_boundary_chain import conform_encoded_boundary

DEST=OUT/'layer-delivery'
ROLES=['adventitia','muscularis','mucosa']
REGIONS=['collecting-wall','upj','ureter-wall']

def combine(sf,uses):return np.concatenate([sf[u['surface']]if u['side']==1 else sf[u['surface']][:,[0,2,1]]for u in uses])
def point_key(tri):return tuple(sorted(map(tuple,tri)))

def donor(*,coordinate_map=None,destination=None,reflected=False):
    dest=Path(destination or DEST)
    if dest.exists():raise FileExistsError(dest)
    dest.mkdir();began=time.monotonic()
    old,z=bound('renal-delivery-v419');k=old['kidneys'][0]
    assert (k['sex'],k['side'])==('male','left')
    v=z['male_left_sourcePositions'];sf={s['id']:z[s['facesArray']]for s in k['surfaces']}
    _,q=bound('renal-ureter-fluid-v317');wv=q['positions'];wf=q['wholeWallFaces'];lookup={point_key(wv[f]):i for i,f in enumerate(wf)}
    face_regions={};face_neighbors={}
    for region in REGIONS:
        for u in k['cells'][region]:
            if u.get('neighbor')in REGIONS:continue
            for f in sf[u['surface']]:
                index=lookup[point_key(v[f])];face_regions[index]=region;face_neighbors[index]=u.get('neighbor')
    assert len(face_regions)==len(wf)
    bank=CommonPositionBank(wv);annuli=[]
    for sid in ['pelvic-ureteric-annulus','retained-ureter-root-annulus']:
        ids=np.array([bank.insert(p)for p in v],int);annuli.append(ids[sf[sid]])
    av=bank.positions;internal=np.concatenate(annuli)
    # Carrier identities belong to the source, not to distances measured after
    # a recipient's free-ureter course has been transported.
    distance=lambda pts,f:pcu.closest_points_on_mesh(np.ascontiguousarray(pts),wv,np.ascontiguousarray(f,dtype=np.int32))[0]
    tri=wv[wf];outer=(distance(tri.reshape(-1,3),q['outerClassifier']).reshape(-1,3).max(1)<2e-12)&(distance(tri.mean(1),q['outerClassifier'])<2e-12)
    inner=(distance(tri.reshape(-1,3),q['innerClassifier']).reshape(-1,3).max(1)<2e-12)&(distance(tri.mean(1),q['innerClassifier'])<2e-12)
    caps=wv[q['innerClassifier']];common=distance(caps.reshape(-1,3),q['outerClassifier']).reshape(-1,3).max(1)<2e-12
    for cap in caps[common]:
        c=cap.mean(0);nrm=np.cross(cap[1]-cap[0],cap[2]-cap[0]);nrm/=np.linalg.norm(nrm)
        end=(np.linalg.norm(tri.mean(1)-c,axis=1)<.006)&(np.max(abs((tri-c)@nrm),axis=1)<3e-9);outer&=~end;inner&=~end
    ambiguous=outer&inner;outer&=~ambiguous;inner&=~ambiguous
    shared=np.intersect1d(np.unique(wf[outer]),np.unique(wf[inner]))
    if coordinate_map is not None:
        v=coordinate_map(v);av=coordinate_map(av);wv=av[:len(wv)]
    if reflected:wf=wf[:,[0,2,1]]
    p,t,used,index,te=tetrahedralize(av,wf,internal)
    print('Literal collecting boundary and both regional annuli constrained',len(t),flush=True)
    pmap={tuple(x):i for i,x in enumerate(p)}
    barriers={tuple(sorted(pmap[tuple(av[i])]for i in f))for f in internal}
    source={tuple(sorted(pmap[tuple(wv[i])]for i in f)):j for j,f in enumerate(wf)}
    owners=defaultdict(list)
    for j,tet in enumerate(t):
        for face in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)]:owners[tuple(sorted(tet[list(face)]))].append(j)
    edges=[];seeds=defaultdict(set)
    for face,ts in owners.items():
        if len(ts)==2 and face not in barriers:edges.extend([ts,ts[::-1]])
        if len(ts)==1:seeds[ts[0]].add(face_regions[source[face]])
    edges=np.array(edges,int);graph=coo_matrix((np.ones(len(edges)),edges.T),shape=(len(t),len(t))).tocsr();n,component=connected_components(graph,directed=False)
    labels=defaultdict(set)
    for tid,regions in seeds.items():labels[int(component[tid])]|=regions
    if n!=3 or any(len(labels[i])!=1 for i in range(n)):raise ValueError('Physical regional annuli do not separate exactly the three original wall regions')
    regions=np.array([next(iter(labels[int(i)]))for i in component])
    r=layers(av,wf,wf[outer],wf[inner],ROLES,tetrahedra=(p,t),tetra_regions=regions,origin=k['geometryOrigin'],free_source_vertices=shared,method='boundary-harmonic')
    bank=CommonPositionBank(v);mapping=np.array([bank.insert(x)for x in r['positions']],int);pairs=defaultdict(list)
    def emit(fs,a,b=None):
        if b is not None and a>b:a,b=b,a;fs=fs[:,[0,2,1]]
        pairs[(a,b)].extend(fs)
    for name,uses in k['cells'].items():
        if name in REGIONS:continue
        for u in uses:
            if u.get('neighbor')in REGIONS:continue
            if u.get('neighbor')and name>u['neighbor']:continue
            fs=sf[u['surface']]if u['side']==1 else sf[u['surface']][:,[0,2,1]]
            if reflected:fs=fs[:,[0,2,1]]
            emit(fs,name,u.get('neighbor'))
    for name,uses in r['cells'].items():
        for u in uses:
            if u.get('neighbor')and name>u['neighbor']:continue
            fs=r['surfaces'][u['surface']];tags=r['source_faces'][u['surface']]
            if u['side']<0:fs=fs[:,[0,2,1]]
            fs=mapping[fs]
            if u.get('neighbor'):emit(fs,name,u['neighbor']);continue
            if np.any(tags<0):raise ValueError('An external coat face lost original wall ancestry')
            neighbors=np.array([face_neighbors[int(j)]for j in tags],object)
            for neighbor in set(neighbors):emit(fs[neighbors==neighbor],name,neighbor)
    surfaces={};cells={name:[]for name in [*[n for n in k['cells']if n not in REGIONS],*r['cells']]}
    for i,((a,b),fs)in enumerate(pairs.items()):
        sid='renal-layer-interface:'+str(i);surfaces[sid]=np.asarray(fs,int);cells[a].append(dict(surface=sid,side=1,neighbor=b))
        if b is not None:cells[b].append(dict(surface=sid,side=-1,neighbor=a))
    # The reviewed rounded hilum replaces four separate capsule perforations;
    # the native inner/outer capsule carriers and all other anatomy stay fixed.
    cv=Mesh.load(OUT/'capsule-v429/capsule.npz');cvp=coordinate_map(cv.vertices)if coordinate_map is not None else cv.vertices;cm=np.array([bank.insert(x)for x in cvp],int)
    for u in cells['capsule']:surfaces.pop(u['surface'])
    surfaces['capsule:rounded-hilum']=cm[cv.faces[:,[0,2,1]]if reflected else cv.faces];cells['capsule']=[dict(surface='capsule:rounded-hilum',side=1)]
    if coordinate_map is not None:
        cortex=Mesh.load(OUT/'cortex-storage-cleanup-common/cortex.npz');cm=np.array([bank.insert(x)for x in coordinate_map(cortex.vertices)],int)
        for u in cells['cortex']:surfaces.pop(u['surface'])
        surfaces['cortex:storage-cleaned']=cm[cortex.faces[:,[0,2,1]]if reflected else cortex.faces];cells['cortex']=[dict(surface='cortex:storage-cleaned',side=1)]
    positions=bank.positions
    # Shared chart refinement at a papillary edge also belongs to the adjoining
    # original medullary free surface, never to a second overlapping cap.
    from geometry_tools.algorithms.source_face_complement import conform_surface_boundary_edges
    surfaces,conformity=conform_surface_boundary_edges(positions,surfaces,cells,tolerance=2e-12)
    rows=[]
    for name,uses in cells.items():
        faces=combine(surfaces,uses);mesh=Mesh(positions,faces,FRAME)
        genus=11 if name.startswith('collecting-wall:')else 1 if name.startswith(('upj:','ureter-wall:'))else 0 if name=='capsule'or name.startswith(('pyramid:','papilla:'))else None
        topo=topology(mesh);ok=topo['closed_oriented_manifold']and(genus is None or[c['closed_orientable_genus']for c in topo['components']]==[genus])
        rows.append(dict(cell=name,passed=ok,topology=topo));print(name,len(faces),'topology',ok,flush=True)
        if not ok:mesh.save(dest/'failure.npz');atomic_json(dest/'failure.json',rows[-1]);raise ValueError('Donor tissue topology failed: '+name)
    arrays={'positions':positions};meta=[]
    for i,(sid,fs)in enumerate(surfaces.items()):key='surface_'+str(i);arrays[key]=fs;meta.append(dict(id=sid,array=key))
    atomic_npz(dest/'donor.npz',**arrays);atomic_json(dest/'donor.json',dict(surfaces=meta,cells=cells,topology=rows,construction=te,boundary_refinement=conformity,seconds=time.monotonic()-began,source_report_sha256=digest_file(ROOT/'reports/renal-delivery-v419.json'),donor_sha256=digest_file(dest/'donor.npz'),installed=False))
    return positions,surfaces,cells


def recipient(sex,side,*,replay_layers=False):
    import importlib.util
    spec=importlib.util.spec_from_file_location('renal_existing_delivery',ROOT/'scripts/stage-renal-delivery-v416.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    old,z=bound('renal-delivery-v419');stem='donor-clean'if(DEST/'donor-clean.json').exists()else'donor';r=json.loads((DEST/(stem+'.json')).read_text());data=np.load(DEST/(stem+'.npz'));pos=data['positions'];sf={s['id']:data[s['array']]for s in r['surfaces']};cells=r['cells']
    k=next(q for q in old['kidneys']if(q['sex'],q['side'])==(sex,side));label=sex+'-'+side;out=DEST/label;out.mkdir();a=json.loads((ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json').read_text())
    _,mz=bound('medullary-collecting-cells-v300');outer=mz['native_kidney_v'][np.unique(mz['native_kidney_f'])]
    _,uz=bound('renal-outlets-v257');jr,_=bound('renal-upj-native-v316');curve=helper.ordered_rings(uz,'left')[int(jr['targetStation']):]
    q=helper.L.load(a,k['sourceKidneyPartId']);native=q['v'][np.unique(q['f'])]
    if (sex,side)==('male','left'):proposed=pos.copy()
    else:
        if sex=='male':target=helper.ordered_rings(uz,side)
        else:
            uid='LOCAL:female:'+side[0]+'_ureter:'+('186'if side=='left'else'185');uq=helper.L.load(a,uid);uv,iv=np.unique(uq['v'],axis=0,return_inverse=True);rings=helper.F.trunk(uv,iv[uq['f']]);target=np.array([uv[g].mean(0)for g in rings])
        proposed,fit=helper.recipient_map(pos,outer,native,curve,target,side);assert abs(fit['uniformScale']-k['fit']['uniformScale'])<1e-14
        if any(k['fit'].get('wholeKidneyTranslationMeters',[0,0,0])):
            # The installed right-side placement was corrected from the
            # previous encoded bank; replay that recorded order literally.
            encoded_origin=np.array(k['geometryOrigin']);proposed=(proposed-encoded_origin).astype('<f4').astype(float)+encoded_origin
        proposed=helper.recipient_adjustment(proposed,k['fit'])
    if replay_layers:
        def mapped(points):
            result,fit=helper.recipient_map(points,outer,native,curve,target,side)
            if any(k['fit'].get('wholeKidneyTranslationMeters',[0,0,0])):
                o=np.array(k['geometryOrigin']);result=(result-o).astype('<f4').astype(float)+o
            return helper.recipient_adjustment(result,k['fit'])
        replay_path=DEST/(label+'-layer-replay')
        if(replay_path/'donor.json').exists():
            replay=json.loads((replay_path/'donor.json').read_text());rd=np.load(replay_path/'donor.npz');proposed=rd['positions'];sf={s['id']:rd[s['array']]for s in replay['surfaces']};cells=replay['cells']
        else:proposed,sf,cells=donor(coordinate_map=mapped,destination=replay_path,reflected=side=='right')
    prefix=z['male_left_sourcePositions'];expected=z[sex+'_'+side+'_sourcePositions'];err=float(np.linalg.norm(proposed[:len(prefix)]-expected,axis=1).max());assert err<1e-13,(label,err)
    origin=np.array(k['geometryOrigin'])
    if(sex,side)==('male','right'):
        root=next(s for s in k['surfaces']if s['id']=='retained-ureter-root-annulus');center=proposed[np.unique(z[root['facesArray']])].mean(0)
        origin=np.rint(center*1024)/1024
    sf={sid:(f[:,[0,2,1]]if side=='right'and not replay_layers else f.copy())for sid,f in sf.items()}
    sf,_,quant=conform_quantized_surfaces(proposed,sf,origin)
    sf,_,contract=contract_quantized_intersection_edges(proposed,sf,origin,len(prefix),len(proposed),maximum_edge_m=1e-8)
    diags={};encoding_error=float(np.linalg.norm((proposed-origin).astype('<f4').astype(float)+origin-proposed,axis=1).max())
    for sid,f in sf.items():
        try:sf[sid],diags[sid]=repair_diagonals(proposed,f,origin,plane_tolerance=max(2e-12,4*encoding_error))
        except ValueError:
            atomic_npz(out/'diagonal-failure.npz',positions=proposed,faces=f,origin=origin)
            raise
    stored=(proposed-origin).astype('<f4').astype(float)+origin
    conformity={}
    for sid,f in sf.items():sf[sid],_,conformity[sid]=conform_encoded_boundary(stored,f)
    rounding={}
    if(sex,side)==('male','right'):
        from .storage_interface_repair import repair
        stored,rounding=repair(stored,sf,cells,origin,len(prefix),['upj:mucosa','ureter-wall:adventitia','ureter-wall:muscularis','ureter-wall:mucosa'])
    words=normal_words(stored,sf);part_ids={name:'RENALDETAIL:'+sex+':'+side+':'+name for name in cells};packed=pack_finite_cells(stored,sf,words,cells,part_ids,origin)
    bank=out/'renal.bin';bank.write_bytes(packed.data);Path(str(bank)+'.gz').write_bytes(gzip.compress(packed.data,mtime=0))
    atomic_json(out/'packed.json',dict(parts=packed.cells,part_ids=part_ids,origin=origin,bank=str(bank),frame='model:'+('local-ta98'if sex=='male'else'local-female-ta98')+':assembled-metres'))
    meshes=decode_packed(packed,'model:'+('local-ta98'if sex=='male'else'local-female-ta98')+':assembled-metres',origin);checks=[]
    for name,mesh in meshes.items():
        genus=11 if name.startswith('collecting-wall:')else 1 if name.startswith(('upj:','ureter-wall:'))else 0 if name=='capsule'or name.startswith(('pyramid:','papilla:'))else None
        topo=topology(mesh)
        if not topo['closed_oriented_manifold']or genus is not None and[c['closed_orientable_genus']for c in topo['components']]!=[genus]:
            mesh.save(out/'failure.npz');atomic_json(out/'failure.json',dict(cell=name,topology=topo));raise ValueError(label+' encoded topology failed: '+name)
        result=audit(mesh,require_closed=True,expected_genus=[genus]if genus is not None else None,max_pairs=1000000)
        if not result['passed']:atomic_json(out/'failure.json',dict(cell=name,check=result));raise ValueError(label+' encoded audit failed: '+name)
        mesh.save(out/(name.replace(':','-')+'.npz'));checks.append(dict(cell=name,check=result));print(label,name,'passed',flush=True)
    arrays={'positions':stored};surface_rows=[]
    for i,(sid,f)in enumerate(sf.items()):key='surface_'+str(i);arrays[key]=f;surface_rows.append(dict(id=sid,array=key))
    atomic_npz(out/'arrangement.npz',**arrays)
    result=dict(sex=sex,side=side,bank=str(bank),bank_sha256=digest_file(bank),packedParts=packed.cells,partIds=part_ids,cells=cells,surfaces=surface_rows,geometryOrigin=origin.tolist(),fit=k['fit'],checks=checks,passed=True,source_control_transform_error_m=err,quantization=quant,derived_contractions=contract,diagonal_repairs=diags,encoded_boundary_refinement=conformity,internal_rounding_repair=rounding,new_images=0,new_uv_arrays=0,installed=False)
    atomic_json(out/'result.json',result);return result


def clean_cortex():
    """Use the bounded, independently audited source-vertex-only cleanup.

    Rebinding its old per-triangle provenance permits removal of numerical
    slivers. It is not a remodel or a signed-field resampling of the organ.
    """
    path=OUT/'cortex-storage-cleanup-common';a=json.loads((path/'audit.json').read_text());mesh=Mesh.load(path/'cortex.npz')
    assert a['passed']and a['geometry_sha256']==mesh.geometry_sha
    r=json.loads((DEST/'donor.json').read_text());z=np.load(DEST/'donor.npz');bank=CommonPositionBank(z['positions']);sf={s['id']:z[s['array']]for s in r['surfaces']}
    for u in r['cells']['cortex']:sf.pop(u['surface'])
    ids=np.array([bank.insert(x)for x in mesh.vertices],int);sf['cortex:storage-cleaned']=ids[mesh.faces];r['cells']['cortex']=[dict(surface='cortex:storage-cleaned',side=1)]
    arrays={'positions':bank.positions};rows=[]
    for i,(sid,f)in enumerate(sf.items()):key='surface_'+str(i);arrays[key]=f;rows.append(dict(id=sid,array=key))
    atomic_npz(DEST/'donor-clean.npz',**arrays);r.update(surfaces=rows,cortex_cleanup=dict(source_audit_sha256=digest_file(path/'audit.json'),source_mesh_sha256=digest_file(path/'cortex.npz'),source_vertex_simplification_tolerance_m=1e-7,no_signed_field_resampling=True),donor_sha256=digest_file(DEST/'donor-clean.npz'))
    atomic_json(DEST/'donor-clean.json',r)


if __name__=='__main__':donor()
