"""Source-only floor-of-mouth attachment annotation and context review.

Does not reconstruct gingiva, move fixed anatomy or create a frenulum.
"""
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.interpolate import CubicSpline
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH,triangle_pairs
from tongue_whole_input import WORK,CATALOGUE

def bound(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=digest_file(p))

def main(output,annotate=False):
    began=time.monotonic()
    if output.exists() or not output.is_relative_to(WORK):raise ValueError('Use a new stage')
    gp=WORK/'floor-of-mouth-guides-v1/guides.json';guides=json.loads(gp.read_text())
    cp=Path(guides['contract']['path']);contract=json.loads(cp.read_text())
    for d in [guides['contract'],guides['gingiva']['mesh']]:
        if digest_file(d['path'])!=d['sha256']:raise ValueError('Bound source changed')
    gum=Mesh.load(guides['gingiva']['mesh']['path']);e=np.r_[gum.faces[:,[0,1]],gum.faces[:,[1,2]],gum.faces[:,[2,0]]]
    _,lab=connected_components(coo_matrix((np.ones(len(e)),(e[:,0],e[:,1])),shape=(len(gum.vertices),)*2).tocsr(),directed=False)
    valid=[k for k in np.unique(lab)if (lab[gum.faces[:,0]]==k).sum()>1000]
    lower=min(valid,key=lambda k:gum.vertices[lab==k,1].mean());face_ids=np.flatnonzero(lab[gum.faces[:,0]]==lower)
    used,inv=np.unique(gum.faces[face_ids],return_inverse=True)
    component=Mesh(gum.vertices[used],inv.reshape(-1,3),gum.frame,
        {k:v[used]for k,v in gum.attributes.items()},source=dict(parent=bound(guides['gingiva']['mesh']['path']),source_face_ids=face_ids.tolist(),source_vertex_ids=used.tolist()))
    output.mkdir();lower_path=output/'actual-lower-gingiva.npz';component.save(lower_path)
    atomic_npz(output/'lower-source-map.npz',source_faces=face_ids,source_vertices=used,frame=np.array(gum.frame))
    source_by_id={r['part_id']:r for r in contract['source_exports']}
    render=[dict(mesh=str(lower_path),name='Actual lower gingiva - unchanged native surface',color=[.74,.39,.42,1])]
    for pid,color in [('ZA:Sublingual gland.l',[.85,.70,.30,.8]),('ZA:Sublingual gland.r',[.85,.70,.30,.8]),
        ('ZA:Submandibular duct.l',[.15,.75,.85,1]),('ZA:Submandibular duct.r',[.15,.75,.85,1])]:
        row=source_by_id[pid]
        if digest_file(row['mesh']['path'])!=row['mesh']['sha256']:raise ValueError('Bound source changed')
        render.append(dict(mesh=row['mesh']['path'],name=row['name'],color=color))
    # The tongue is omitted from this source-inspection scene solely to expose
    # the lingual gingiva/duct endpoints; it remains unchanged in the atlas.
    annotation=None
    if annotate:
        index=BVH(component);controls=[];gland_rows={};source_crossings=[]
        for side in ['l','r']:
            row=source_by_id['ZA:Sublingual gland.'+side];gland=Mesh.load(row['mesh']['path'])
            centers=gland.triangles.mean(1);normals=np.cross(gland.triangles[:,1]-gland.triangles[:,0],gland.triangles[:,2]-gland.triangles[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
            roof=np.flatnonzero(normals[:,1]>.35);side_rows=[]
            for z in np.linspace(centers[roof,2].min(),centers[roof,2].max(),7):
                near=roof[np.argsort(np.abs(centers[roof,2]-z))[:20]];i=int(near[np.argmax(centers[near,1])]);q=centers[i]
                hit=index.closest_point(q);side_rows.append(dict(source_part_id=row['part_id'],source=bound(row['mesh']['path']),
                    source_face=i,source_barycentric=[1/3]*3,source_point_m=q.tolist(),
                    gingiva_face=int(hit['face_id']),gingiva_source_face=int(face_ids[hit['face_id']]),
                    gingiva_barycentric=hit['barycentric'].tolist(),point_m=hit['point'].tolist(),source_gap_m=float(hit['distance_m'])))
            gland_rows[side]=side_rows
            # Do not repair an existing native gland/gum clash while deriving
            # the FOM. Record literal triangle crossings as a remaining limit.
            count=0;examples=[];other=BVH(gland)
            for pairs in index.candidate_batches(other,0):
                r=triangle_pairs(index.triangles[pairs[:,0]],other.triangles[pairs[:,1]])
                mask=r['relation']=='crossing';count+=int(mask.sum());examples.extend(pairs[mask][:max(0,12-len(examples))].tolist())
            source_crossings.append(dict(first='ZA:Gingiva lower source component',second=row['part_id'],crossing_face_pairs=count,
                examples=examples,interpretation='Existing fixed native source crossing, not a new FOM defect or permitted tissue interface'))
        oral=[r for r in guides['duct_endpoint_loops']if r['classification'].startswith('anterior')]
        center=np.mean([r['center_m']for r in oral],axis=0);center[0]=0
        h=index.closest_point(center);middle=dict(source_part_id='paired anterior-medial duct endpoint mean',
            source_point_m=center.tolist(),gingiva_face=int(h['face_id']),gingiva_source_face=int(face_ids[h['face_id']]),
            gingiva_barycentric=h['barycentric'].tolist(),point_m=h['point'].tolist(),source_gap_m=float(h['distance_m']))
        controls=gland_rows['l']+[middle]+gland_rows['r'][::-1]
        p=np.array([r['point_m']for r in controls]);t=np.r_[0,np.linalg.norm(np.diff(p,axis=0),axis=1).cumsum()]
        if np.any(np.diff(t)<=0):raise ValueError('Repeated guide points need explicit source review')
        spline=CubicSpline(t,p,axis=0,bc_type='natural');queries=spline(np.linspace(0,t[-1],129));witnesses=[]
        for q in queries:
            hit=index.closest_point(q);i=int(hit['face_id']);tri=component.triangles[i];n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n/=np.linalg.norm(n)
            radial=np.array([hit['point'][0],0,hit['point'][2]-.04]);radial/=np.linalg.norm(radial)
            witnesses.append(dict(source_face=int(face_ids[i]),lower_component_face=i,barycentric=hit['barycentric'].tolist(),
                point_m=hit['point'].tolist(),source_surface_residual_m=float(np.linalg.norm(hit['barycentric']@tri-hit['point'])),
                projection_distance_m=float(hit['distance_m']),normal_m=n.tolist(),normal_radial_dot=float(n@radial)))
        if max(r['normal_radial_dot']for r in witnesses)>-.6:raise ValueError('Witness curve left the inward lingual gingiva')
        pts=np.array([r['point_m']for r in witnesses]);steps=np.linalg.norm(np.diff(pts,axis=0),axis=1)
        if steps.max()>.001:raise ValueError('Source-projected curve has an unreviewed >1mm jump')
        atomic_npz(output/'attachment-arc-witnesses.npz',positions=pts,
            source_face=np.array([r['source_face']for r in witnesses]),
            barycentric=np.array([r['barycentric']for r in witnesses]),frame=np.array(gum.frame))
        marked=np.unique([r['lower_component_face']for r in witnesses]);remaining=np.setdiff1d(np.arange(len(component.faces)),marked)
        # Literal complementary face subsets highlight the source witnesses;
        # no offset, duplicate overlapping face or new anatomical mesh is used.
        base_mesh=Mesh(component.vertices,component.faces[remaining],component.frame,component.attributes,source={'source':bound(lower_path),'inspection_only':True})
        stripe=Mesh(component.vertices,component.faces[marked],component.frame,component.attributes,source={'source':bound(lower_path),'inspection_only':True})
        base_mesh.save(output/'lower-gingiva-without-witness-faces.npz');stripe.save(output/'actual-witness-faces.npz')
        render[0]['mesh']=str(output/'lower-gingiva-without-witness-faces.npz')
        render.append(dict(mesh=str(output/'actual-witness-faces.npz'),name='Proposed source-bound gingiva-FOM joint witnesses (not an accepted anatomical margin)',color=[1,.55,.06,1]))
        duct_identity=[];atlas=json.loads(CATALOGUE.read_text())
        if digest_file(CATALOGUE)!=contract['source_catalogue']['sha256']:raise ValueError('Current catalogue changed; do not mix fixed source family')
        for side in ['l','r']:
            pid='ZA:Submandibular gland.'+side;part=next(p for p in atlas['parts']if p['id']==pid)
            native=AtlasReader().read_part(atlas,part,frame=gum.frame,catalogue_path=CATALOGUE,catalogue_sha=digest_file(CATALOGUE))
            path=output/('actual-submandibular-gland-'+side+'.npz');native.save(path);b=BVH(native)
            for r in guides['duct_endpoint_loops']:
                if r['part_id']!='ZA:Submandibular duct.'+side:continue
                hit=b.closest_point(r['center_m']);duct_identity.append(dict(part_id=r['part_id'],endpoint=r['classification'],source_loop=r['loop'],
                    submandibular_gland=bound(path),closest_gland_face=int(hit['face_id']),closest_gland_barycentric=hit['barycentric'].tolist(),
                    closest_gland_point_m=hit['point'].tolist(),gap_to_gland_m=float(hit['distance_m']),
                    accepted_endpoint_identity='oral'if r['classification'].startswith('anterior')else'gland',
                    actual_through_mucosa_ostium_built=False,finite_duct_wall_built=False))
        annotation=dict(schema='floor-mouth-attachment-witnesses/v1',source_gingiva=guides['gingiva']['mesh'],
            source_map=bound(output/'lower-source-map.npz'),curve=bound(output/'attachment-arc-witnesses.npz'),
            construction_controls=controls,witnesses=witnesses,
            marked_source_faces=face_ids[marked].tolist(),maximum_curve_step_m=float(steps.max()),
            maximum_projection_distance_m=max(r['projection_distance_m']for r in witnesses),
            anatomical_boundary_status='Proposed joint line inside demonstrated native lingual gingival surface. Source has no presegmented mucogingival/FOM border or fascia insertion; exact attachment height remains provisional.',
            arc_extent='Anterior arc and paired lateral fields beside the fixed sublingual gland roofs; not the full posterior FOM/piriform continuation.',
            source_gland_gum_crossings=source_crossings,duct_endpoint_identity=duct_identity,
            next_construction=dict(field='Curved finite mucosa plus underlying FOM fascia between this literal arc and existing 64-point tongue reflection',
                posterior='Retain paired lateral continuation seams, not a solid posterior disk',
                ducts='Finite annular duct collars through the field at actual anterior-medial source loops; lumen open',
                frenulum='Central continuous mucosa/fascia reflection in present resting pose; no free cord, triangle, or forced tongue elevation',
                unresolved=['Donor-independent fascia/mucosa thickness selection','Pre-existing fixed native gland/gum crossings','Reciprocal opening/replacement of current tongue closure charts','Posterolateral continuation boundary']),
            numerical_surface_witnesses_ready=True,exact_anatomical_attachment_accepted=False,
            installed=False,no_existing_mesh_modified=True,no_new_tissue_constructed=True)
        atomic_json(output/'attachment-witnesses.json',annotation)
    spec=dict(frame=gum.frame,resolution=768,wall_seconds=45,meshes=render,
        views=[dict(name='lower-gingiva-superior',direction=[0,1,0]),
            dict(name='lower-gingiva-lingual-oblique',direction=[0,.75,-.65]),
            dict(name='lower-gingiva-right',direction=[1,.15,0])],
        qualification='Literal source-context inspection only. Highlighted source facets indicate a proposed gingiva-FOM joint line, not a donor-verified anatomical boundary. No FOM roof, frenulum or through-mucosa opening constructed.')
    atomic_json(output/'review-spec.json',spec)
    summary=dict(schema='floor-mouth-source-inspection/v1',contract=bound(cp),guides=bound(gp),
        lower_gingiva=dict(mesh=bound(lower_path),source_map=bound(output/'lower-source-map.npz'),
            source_component=int(lower),source_faces=len(face_ids),vertices=len(used),bounds_m=component.bounds.tolist()),
        actual_duct_endpoint_loops=guides['duct_endpoint_loops'],ready_to_construct=False,
        attachment_annotation=None if annotation is None else bound(output/'attachment-witnesses.json'),
        no_existing_mesh_modified=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'inspection.json',summary)
    print(json.dumps(dict(output=str(output),lower_gingiva_faces=len(face_ids),elapsed_seconds=time.monotonic()-began)),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--annotate',action='store_true');a=ap.parse_args();main(a.output.resolve(),a.annotate)
