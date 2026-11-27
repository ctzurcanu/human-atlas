"""Independent coarse isosurfaces from private anatomical shape references.

Original faces, vertex sequences, UVs and pixels never enter the output bank.
The source is a local construction guide, not redistribution authority.
"""
from pathlib import Path
import argparse,json,sys,time,runpy
import numpy as np
import open3d as o3d
import manifold3d as md
from scipy.ndimage import gaussian_filter
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'scripts'),str(ROOT/'server')]
from ta98_glb_geometry import load_roots
SOURCE_CIDS=runpy.run_path(str(ROOT/'server/heart-reference-types.py'))['SOURCE_CIDS']
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.geometry_instance import instance_pose
from geometry_tools.algorithms.regular_grid_isosurface import extract

ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--material-end-caps',action='store_true');args=ap.parse_args()
started=time.monotonic();out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
reference=ROOT/'.local-models/models/Male_2024.glb'
metadata_path=ROOT.parent/'3d/ext/meshes/964db2dd4f98052f03baa9ca5f2dbcae/assets.json'
metadata=json.loads(metadata_path.read_text());wanted=set(SOURCE_CIDS.values())
rows={r['sni'][0]['id']:r for r in metadata['structures']if r.get('cid')in wanted}
sources=load_roots(reference,set(rows),allow_zero_weight_targets=True)
landmarks=json.loads((ROOT/'.local-models/geometry-workbench/heart/cardiac-label-landmarks-v604/report.json').read_text())
transform=np.asarray(landmarks['matrix']);by_cid={}
for m in sources:
    cid=rows[m['node']]['cid'];by_cid.setdefault(cid,[]).append(m)
atlas_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
atlas=json.loads(atlas_path.read_text());host=next(p for p in atlas['parts']if p['id']=='TA98HEART:male:epicardium')
R,scale=instance_pose(host);origin=np.asarray(host['geometryOrigin'])
records=[]
for term,cid in SOURCE_CIDS.items():
    if args.material_end_caps and term in {'A12.1.01.002','A12.1.03.002','A12.1.02.008','A12.1.08.002','A12.1.08.006'}:
        records.append(dict(term=term,cid=cid,status='native-region-or-paired-wall-required'));continue
    pieces=by_cid[cid];vv=[];ff=[];offset=0
    for m in pieces:
        v=m['unique']@transform[:3,:3].T+transform[:3,3];vv.append(v);ff.append(m['welded']+offset);offset+=len(v)
    v=np.concatenate(vv);f=np.concatenate(ff);centroid=v.mean(0)
    _,_,basis=np.linalg.svd(v-centroid,full_matrices=False)
    p=(v-centroid)@basis.T
    cap_records=[]
    if args.material_end_caps:
        directed=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);keys,inv,uses=np.unique(np.sort(directed,axis=1),axis=0,return_inverse=True,return_counts=True)
        boundary=directed[uses[inv]==1];adj={int(a):int(b)for a,b in boundary}
        if len(adj)!=len(boundary)or set(adj)!=set(adj.values()):
            records.append(dict(term=term,cid=cid,status='unresolved-source',reason='Source boundary is not a collection of oriented simple material ends'));continue
        unused=set(adj);caps=[]
        while unused:
            start=min(unused);ring=[start];here=adj[start]
            while here!=start:
                ring.append(here);here=adj[here]
                if len(ring)>len(adj):raise ValueError('Non-simple material end')
            unused.difference_update(ring);q=p[ring];_,sing,axes=np.linalg.svd(q-q.mean(0),full_matrices=False);uv=(q-q.mean(0))@axes[:2].T
            tri=np.asarray(md.triangulate([uv],allow_convex=False),dtype=int);cap=np.asarray(ring)[tri]
            if any((a==ring[0]and b==ring[1])for face in cap for a,b in zip(face,np.roll(face,-1))):cap=cap[:,::-1]
            caps.append(cap);cap_records.append(dict(edges=len(ring),planarity_ratio=float(sing[-1]/max(sing[0],1e-15))))
        if caps:f=np.concatenate([f,*caps])
    edges=np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1)
    _,counts=np.unique(edges,axis=0,return_counts=True);closed=bool(np.all(counts==2))
    widths=np.ptp(p,axis=0)
    # Limit sampling to a fixed small guide grid, with at least six samples
    # across the thinnest source extent where the source supplies a volume.
    pitch=max(min(float(widths.max()/72),float(np.min(widths[widths>1e-8])/8)),.00006)
    half_sheet=max(.00008,pitch*.7) if not closed else 0.
    if args.material_end_caps and not closed:
        records.append(dict(term=term,cid=cid,status='unresolved-source',reason='Material end closure did not establish a finite guide'));continue
    padding=max(pitch*3,half_sheet*2);lo=p.min(0)-padding;hi=p.max(0)+padding
    shape=np.ceil((hi-lo)/pitch).astype(int)+1
    if np.prod(shape)>900000:
        pitch*=float((np.prod(shape)/900000)**(1/3));shape=np.ceil((hi-lo)/pitch).astype(int)+1
    grid=np.indices(tuple(shape)).reshape(3,-1).T*pitch+lo
    scene=o3d.t.geometry.RaycastingScene(nthreads=1)
    scene.add_triangles(o3d.core.Tensor(p.astype('float32')),o3d.core.Tensor(f.astype('uint32')))
    query=o3d.core.Tensor(grid.astype('float32'))
    field=(scene.compute_signed_distance(query,nthreads=1,nsamples=3).numpy() if closed else
           scene.compute_distance(query,nthreads=1).numpy()-half_sheet).reshape(tuple(shape))
    # A light construction smoothing removes voxel stair steps. It does not
    # preserve any source control sequence or source triangulation.
    field=gaussian_filter(field,.35,mode='nearest')
    if not field.min()<0<field.max():
        records.append(dict(term=term,cid=cid,status='unresolved-source',closed_source=closed,reason='No finite independently sampled interval'));continue
    verts,faces,contour=extract(field,lo,pitch)
    # Bound the independent guide's triangle/memory cost without restoring
    # original source controls, source faces or image coordinates.
    solid=md.Manifold(md.Mesh64(verts,faces.astype('uint64')))
    if str(solid.status())=='Error.NoError':
        reduced=solid.set_tolerance(pitch*.08).simplify(pitch*.25).to_mesh64()
        verts=np.asarray(reduced.vert_properties[:,:3]);faces=np.asarray(reduced.tri_verts,dtype=int)
    world=verts@basis+centroid
    # Export into the literal shared heart-bank coordinate frame; the female
    # applies exactly its current proper similarity to these same controls.
    local=(world-origin)@R/scale
    volume=float(np.einsum('ij,ij->i',local[faces[:,0]],np.cross(local[faces[:,1]],local[faces[:,2]])).sum()/6)
    if volume<0:faces=faces[:,::-1]
    mesh=Mesh(local,faces.astype('int64'),'heart-native-local-metres');mesh.save(out/(term+'.npz'))
    records.append(dict(term=term,cid=cid,status='constructed-guide',source_triangles=len(f),triangles=len(faces),
                        source_vertices=len(v),vertices=len(verts),pitch_metres=pitch,closed_source=closed,
                        illustrative_sheet_half_thickness_m=half_sheet,grid_samples=len(grid),
                        explicit_material_end_caps=cap_records,
                        signed_volume_m3=abs(volume),mesh_sha256=digest_file(out/(term+'.npz')),
                        direct_source_faces_copied=0,direct_source_controls_copied=0,
                        independent_shape_not_anatomical_acceptance=True))
    print(term,cid,'reference',len(f),'->',len(faces),'triangles;',round(pitch*1000,3),'mm guide;',round(time.monotonic()-started,2),'s',flush=True)
atomic_json(out/'construction.json',dict(schema='cardiac-independent-reference-shapes/v1',records=records,
    source_sha256=digest_file(reference),registration_rms_m=landmarks['fit']['landmark_rms_m'],
    source_catalogue_sha256=digest_file(atlas_path),direct_source_faces_copied=0,source_uvs_copied=0,
    source_pixels_copied=0,elapsed_seconds=time.monotonic()-started,installed=False,anatomical_acceptance=False))
