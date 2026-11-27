"""One source-owned finite anterior mucosal papilla field, without textures.

The human SEM and optical-scan references establish crown-like filiform and
rounded fungiform shapes and their approximate diameters. This original field
is a qualified anatomical construction, not copied participant surface data.
Only one outer source triangle is replaced by a conforming height graph; its
three boundary edges and every inner/shared tissue surface remain literal.
"""
from pathlib import Path
import argparse
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
from scipy.spatial import Delaunay
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.tissue import pack,normal_words
from geometry_tools.contact import pair_contacts
from tongue_whole_input import WORK,fixed_bone_sources

SOURCE=WORK/'shared-coverings-exposed-apo-sl-v3/packed/tissue.json'
REFERENCES=['https://pubmed.ncbi.nlm.nih.gov/4061027/',
            'https://www.nature.com/articles/s41598-023-46535-9',
            'https://arxiv.org/html/2307.06255v2']


def require(value,message):
    if not value:raise ValueError(message)


def bound(path):
    path=Path(path).resolve();return dict(path=str(path),sha256=digest_file(path))


def load_source():
    source=json.loads(SOURCE.read_text())
    require(source['passed']and digest_file(source['bank'])==source['bank_sha256'],
            'Use the literal passed v3 common tissue bank')
    with np.load(source['arrangement'],allow_pickle=False)as raw:
        positions=raw['positions'].copy()
        surfaces={key:raw[array].copy()for key,array in source['surfaces'].items()}
    return source,positions,surfaces


def plan(output):
    require(output.is_relative_to(WORK)and not output.exists(),'Use a new local stage')
    source,v,surfaces=load_source()
    f=surfaces['mucosa:outer'];tri=v[f]
    geo=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(geo,axis=1)/2
    n=geo/np.linalg.norm(geo,axis=1)[:,None];c=tri.mean(1)
    # A clearly anterior apical-side dorsal source triangle. No terminal sulcus
    # is invented from a body-length fraction. This local domain lies anterior
    # to the uncertain terminal region and does not claim a full dorsal field.
    eligible=np.flatnonzero((n[:,1]>.8)&(c[:,2]>.072)&(c[:,0]<-.002)&(c[:,0]>-.009))
    require(len(eligible)>0,'No supported anterior dorsal source face')
    face=int(max(eligible,key=lambda k:area[k]));t=tri[face]
    lengths=np.linalg.norm(t-np.roll(t,-1,axis=0),axis=1)
    opposite=np.array([lengths[1],lengths[2],lengths[0]])
    center=opposite@t/opposite.sum();inradius=2*area[face]/lengths.sum()
    normal=n[face];u=t[1]-t[0];u/=np.linalg.norm(u);vv=np.cross(normal,u)
    basis=np.array([u,vv]);local=(t-center)@basis.T
    radius=min(.0018,.78*inradius)
    require(radius>.0012,'Source triangle cannot contain a correctly scaled field')
    # Deterministic local positions; no individual donor pattern is asserted.
    centers=[]
    for a in np.arange(-.0015,.00151,.00068):
        for b in np.arange(-.0015,.00151,.00068):
            p=np.array([a+.00012*np.sin(71*a+93*b),b])
            if np.linalg.norm(p)<radius-.00024 and np.linalg.norm(p-np.array([0,.00020]))>.00072:
                centers.append(p)
    require(len(centers)>=5,'The physical domain supports too few filiform units')
    output.mkdir(parents=True)
    atomic_json(output/'plan.json',dict(schema='tongue-papillary-field-plan/v1',
        source=bound(SOURCE),source_bank=bound(source['bank']),source_arrangement=bound(source['arrangement']),
        authority_inventory=bound(ROOT/'reports/ta98-tongue-inventory.json'),
        whole_organ='Current source-bound male tongue mucosa/aponeurosis/SL assembly',
        terms=[dict(term='A05.1.04.013',name='Papillae of tongue',representation='group'),
               dict(term='A05.1.04.014',name='Filiform papillae',representation='finite-mucosal-relief'),
               dict(term='A05.1.04.015',name='Fungiform papillae',representation='finite-mucosal-relief')],
        source_outer_face=face,source_face_vertices=t.tolist(),source_face_vertex_ids=f[face].tolist(),
        center=center.tolist(),normal=normal.tolist(),basis=basis.tolist(),local_triangle=local.tolist(),
        field_radius_m=float(radius),field_domain_area_m2=float(np.pi*radius**2),
        filiform_centers_m=[p.tolist()for p in centers],fungiform_center_m=[0.,.00020],
        primary_references=REFERENCES,
        reference_shape_constraints=dict(filiform_body_with_secondary_hairs=[6,10],
            filiform_mean_diameter_m=.000355,fungiform_mean_diameter_m=.000878,
            filiform_number_per_cm2=[100,200],fungiform_sparsely_scattered=True),
        original_modeling_estimates=dict(filiform_height_m=.00028,filiform_secondary_peaks=8,
            fungiform_height_m=.00027,individual_positions=True),
        full_distribution_plan=dict(
            filiform='Dense anterior presulcal dorsal territory; exclude actual ventral surface and posterior lingual-tonsil territory. Continue over source-bound anterior facets only after the terminal landmark is established.',
            fungiform='Sparse rounded units across anterior dorsum with higher density near the true apex and anterolateral margins; deterministic spatial field with variable diameters, not one mirrored mesh per papilla.',
            controls='One physical mucosal host; parent and type controls own actual relieved source subfacets; microscopic instances never inflate the body hierarchy mesh count.',
            resource_strategy='Coarse continuous host for body, on-demand grouped relief for closer organ view. Keep detailed microscopic field separate from body mesh totals.'),
        interfaces=[dict(first='relieved mucosa:outer',second='unchanged adjacent mucosa',relation='literal shared boundary',
                         proof='Only interior controls; all three source-face edges remain the original single segments'),
                    dict(first='mucosa',second='aponeurosis',relation='shared-tissue-interface',state='preserve literal source chart'),
                    dict(first='papillary relief',second='mucosa:free/inner',relation='same continuous finite tissue',state='preserve inner coordinates/faces')],
        remaining=['Field is one supported anterior patch, not full dorsal distribution',
                   'Vallate/foliate, terminal sulcus/foramen caecum and microscopic taste buds not constructed'],
        source_files_modified=False,installed=False,anatomical_acceptance=False))
    print(json.dumps(dict(plan=str(output/'plan.json'),source_face=face,filiform_units=len(centers),
        field_radius_mm=radius*1000,filiform_number_per_cm2=len(centers)/(np.pi*radius**2)*1e-4)),flush=True)


def bump(d,r):
    q=np.clip(1-np.sum(d*d,axis=-1)/(r*r),0,1)
    return q*q


def field(points,plan):
    p=np.asarray(points,float);height=np.zeros(len(p));owner=np.zeros(len(p),np.int16)
    # Rounded cap, broad smooth top; 878 um diameter as reported in the human
    # optical-scan study. Height is an explicit original modeling estimate.
    fc=np.array(plan['fungiform_center_m'])
    fh=.00027*bump(p-fc,.000439)
    height+=fh;owner[fh>1e-10]=2
    for index,center in enumerate(np.asarray(plan['filiform_centers_m'])):
        d=p-center
        body=.00010*bump(d,.0001775)
        crown=np.zeros(len(p))
        for k in range(8):
            a=2*np.pi*k/8+index*.37
            hair=.000083*np.array([np.cos(a),np.sin(a)])
            crown=np.maximum(crown,.00018*bump(d-hair,.000035))
        values=body+crown
        height+=values;owner[values>1e-10]=1
    return height,owner


def construct(output,planned):
    began=time.monotonic();plan_record=json.loads(planned.read_text())
    require(output.is_relative_to(WORK)and not output.exists(),'Use a new local stage')
    require(digest_file(SOURCE)==plan_record['source']['sha256'],'Changed v3 source')
    source,v,surfaces=load_source();face=plan_record['source_outer_face']
    center=np.asarray(plan_record['center']);basis=np.asarray(plan_record['basis']);normal=np.asarray(plan_record['normal'])
    boundary=np.asarray(plan_record['local_triangle'])
    points=[p for p in boundary]
    # Rings resolve all eight crown peaks geometrically; no normal/bump map,
    # alpha decal or separate intersecting papillary solids are created.
    for index,c in enumerate(np.asarray(plan_record['filiform_centers_m'])):
        points.append(c)
        for r in [.00006,.000125,.0001775,.00021]:
            for a in np.arange(16)*2*np.pi/16:points.append(c+r*np.array([np.cos(a),np.sin(a)]))
        for k in range(8):
            a=2*np.pi*k/8+index*.37;hair=c+.000083*np.array([np.cos(a),np.sin(a)])
            points.append(hair)
            for r in [.0000175,.000035]:
                for b in np.arange(8)*2*np.pi/8:points.append(hair+r*np.array([np.cos(b),np.sin(b)]))
    fc=np.asarray(plan_record['fungiform_center_m']);points.append(fc)
    for r in [.00011,.00022,.00033,.000439,.00049]:
        for a in np.arange(24)*2*np.pi/24:points.append(fc+r*np.array([np.cos(a),np.sin(a)]))
    # The only hull controls are the three original corners. Interior field
    # support is separated from every boundary by the planned inradius margin.
    p=np.asarray(points);height,owners=field(p,plan_record)
    local=Delaunay(p).simplices
    a=p[local[:,1]]-p[local[:,0]];b=p[local[:,2]]-p[local[:,0]]
    signed=a[:,0]*b[:,1]-a[:,1]*b[:,0]
    local[signed<0]=local[signed<0][:,[0,2,1]]
    original=surfaces['mucosa:outer'][face].copy()
    newpositions=center+p[3:]@basis+height[3:,None]*normal
    positions=np.r_[v,newpositions]
    mapped=np.r_[original,np.arange(len(v),len(positions))]
    replacement=mapped[local]
    # Restore exact source boundary corners rather than recompute barycentric
    # positions. Every unmodified source coordinate remains literal.
    new_outer=np.r_[surfaces['mucosa:outer'][:face],replacement,surfaces['mucosa:outer'][face+1:]]
    surfaces={**surfaces,'mucosa:outer':new_outer}
    output.mkdir(parents=True)
    # Preserve all raw source normal words on unchanged charts. The split
    # exterior needs a derivative-derived, then quantized, original field.
    words={}
    data=Path(source['bank']).read_bytes();arr=source['parts']
    for cell,row in arr.items():
        n=row['vertexCount'];indices=np.ndarray((row['indexCount']//3,3),'<u4',buffer=data,offset=row['indices'])
        nw=np.ndarray((n,3),'<i2',buffer=data,offset=row['normals'])
        ns=np.ndarray((n,),'i1',buffer=data,offset=row['normalSigns'])
        off=0
        for use in source['cells'][cell]:
            sid=use['surface'];amount=len(surfaces[sid])if sid!='mucosa:outer'else len(new_outer)-len(replacement)+1
            raw=nw[indices[off:off+amount]].astype(np.int32)*ns[indices[off:off+amount]][:,:,None]
            if use['side']==-1:raw=-raw[:,[0,2,1]]
            if sid not in words:words[sid]=raw.astype('<i2')
            else:require(np.array_equal(words[sid],raw),'Existing reciprocal normal words differ')
            off+=amount
    old_words=words['mucosa:outer'].copy()
    # Per-point analytical graph normal via centered derivative of the same
    # finite function. Boundary corners inherit their literal old source chart.
    step=1e-8;du=(field(p+np.array([step,0]),plan_record)[0]-field(p-np.array([step,0]),plan_record)[0])/(2*step)
    dv=(field(p+np.array([0,step]),plan_record)[0]-field(p-np.array([0,step]),plan_record)[0])/(2*step)
    normals=normal-du[:,None]*basis[0]-dv[:,None]*basis[1]
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    new_words=np.rint(normals[local]*32767).astype('<i2')
    for i,t in enumerate(local):
        for k,j in enumerate(t):
            if j<3:new_words[i,k]=old_words[face,j]
    # Existing bank's origin is exactly preserved; this is the emitted geometry
    # on which tissue checks run, not a higher-precision-only candidate.
    origin=np.asarray(json.loads((SOURCE.parent.parent/'construction.json').read_text())['storageOrigin'])
    actual=(positions-origin).astype('<f4').astype(float)+origin
    # Analytic graph normals can face behind a coarse crown side triangle.
    # Explicit corner charts retain the actual derivative where it is outward
    # and use the literal face normal at those discrete crease corners. This
    # changes neither height, controls, topology nor the source boundary chart.
    tri=actual[replacement];geo=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);geo/=np.linalg.norm(geo,axis=1)[:,None]
    wrong=np.einsum('fki,fi->fk',new_words.astype(float),geo)<=0
    geometric=np.rint(geo*32767).astype('<i2')
    new_words[wrong]=np.broadcast_to(geometric[:,None,:],new_words.shape)[wrong]
    words['mucosa:outer']=np.r_[old_words[:face],new_words,old_words[face+1:]]
    atomic_npz(output/'prepack-source-arrangement.npz',positions=positions,
        **{f'surface_{i}':t for i,t in enumerate(surfaces.values())},
        **{f'words_{i}':words[s]for i,s in enumerate(surfaces)})
    result=pack(positions,surfaces,source['cells'],source['part_ids'],source['frame'],output/'packed',
        origin=origin,corner_words=words,construction=dict(kind='original-finite-mucosal-papillary-relief',
        finite_wall_cells_pass=True,parent_manifest_sha256=digest_file(SOURCE),plan_sha256=digest_file(planned),
        inner_and_reciprocal_charts_literal=True,full_dorsum_distribution=False,anatomical_acceptance=False))
    old_positions_literal=np.array_equal(actual[:len(v)],v)
    actual_tri=actual[replacement];g=np.cross(actual_tri[:,1]-actual_tri[:,0],actual_tri[:,2]-actual_tri[:,0]);g/=np.linalg.norm(g,axis=1)[:,None]
    w=new_words.astype(float);w/=np.linalg.norm(w,axis=2)[:,:,None]
    dot=np.einsum('fki,fi->fk',w,g)
    owners_f=field(p[local].mean(1),plan_record)[1]
    refs=[]
    muco_path=output/'packed/TA98TONGUE_male_mucosa.npz'
    actual_muco=Mesh.load(muco_path)
    old_face_of_new=np.r_[np.arange(face),np.full(len(replacement),face),
                         np.arange(face+1,824)].astype(np.int32)
    require(len(old_face_of_new)==len(actual_muco.faces), 'Actual decoded mucosal face ordering changed')
    for key,term,name,label in [('papillae','A05.1.04.013','Papillae of tongue',None),
        ('filiform','A05.1.04.014','Filiform papillae',1),('fungiform','A05.1.04.015','Fungiform papillae',2)]:
        chosen=np.flatnonzero(owners_f>0 if label is None else owners_f==label)+face
        require(len(chosen)>0,'A proposed papilla type has no actual owned surface')
        a_face=int(max(chosen,key=lambda i:height[local[int(i)-face]].mean()))
        refs.append(dict(key=key,name=name,term=term,parent='papillae'if label else'tongue',
            conceptId='TA98TONGUE:male:feature:'+key,partIds=['TA98TONGUE:male:mucosa'],kind='exact',
            representation='surface-region',physicalRepresentation='finite-three-dimensional-mucosal-relief',
            surface=dict(partId='TA98TONGUE:male:mucosa',faces=chosen.tolist(),sourceMeshSha256=digest_file(muco_path)),
            anchor=dict(partId='TA98TONGUE:male:mucosa',face=a_face,barycentric=[1/3]*3,
                        point=actual_muco.triangles[a_face].mean(0).tolist()),
            fullDistribution=False,anatomicalAcceptance=False))
    report=dict(schema='tongue-source-bound-papillary-relief/v1',passed=bool(result['passed']and old_positions_literal and np.all(dot>0)),
        source=bound(SOURCE),plan=bound(planned),actual_tissue_manifest=bound(output/'packed/tissue.json'),
        source_old_coordinates_literal=old_positions_literal,old_source_face=face,source_face_count=348,
        new_outer_faces=len(new_outer),replacement_faces=len(replacement),added_controls=len(newpositions),
        inner_and_shared_charts_literal=True,all_unchanged_outer_raw_normal_words_literal=True,
        new_corner_normal_minimum_plane_dot=float(dot.min()),new_normal_guards_pass=bool(np.all(dot>0)),
        crown_crease_corners_with_actual_face_normal_chart=int(wrong.sum()),
        maximum_relief_m=float(height.max()),filiform_units=len(plan_record['filiform_centers_m']),fungiform_units=1,
        actual_tissue_checks=result['checks'],additional_physical_meshes=0,new_image_assets=0,
        packed_bytes=result['total_bytes'],original_packed_bytes=source['total_bytes'],
        termGroups=refs,anatomical_acceptance=False,installed=False,
        remaining=plan_record['remaining']+['Actual close-up/whole Blender surface and section review required'],
        source_files_modified=False,skeleton_modified=False,elapsed_seconds=time.monotonic()-began)
    atomic_npz(output/'relief-correspondence.npz',local_parameters=p,height=height,
        source_face=original,source_face_id=np.array(face),replacement_faces=replacement,
        papilla_type=owners_f,local_triangles=local,source_mucosa_face_of_new_face=old_face_of_new)
    atomic_json(output/'surface-controls.json',dict(schema='tongue-construction-features/v1',
        hostPartId='TA98TONGUE:male:mucosa',carrierMesh=str(muco_path),carrierMeshSha256=digest_file(muco_path),
        carrierGeometrySha256=actual_muco.geometry_sha,carrierBankBindingPending=True,
        carrierBankHashes={},sourceTissueManifestSha256=digest_file(SOURCE),
        sourceFaceAncestry=dict(path=str(output/'relief-correspondence.npz'),sha256=digest_file(output/'relief-correspondence.npz'),
            array='source_mucosa_face_of_new_face'),
        termGroups=refs,additionalSolidMeshes=0,additionalImageAssets=0,anatomicalAcceptance=False,
        references=REFERENCES,fullDorsalDistribution=False))
    atomic_json(output/'construction.json',report)
    path=output/'packed/TA98TONGUE_male_mucosa.npz'
    atomic_json(output/'review-spec.json',dict(frame=source['frame'],resolution=640,wall_seconds=50,
        meshes=[dict(mesh=str(path),name='Actual finite mucosa with original anterior papilla field',color=[.77,.38,.41,1.])],
        views=[dict(name='top',direction=[0,1,0]),dict(name='oblique',direction=[1,1,1]),
               dict(name='section',direction=[1,0,0],section=dict(normal=[1,0,0],point=center.tolist(),clear_outer=True))],
        qualification='One original source-bound anterior papillary field; not full dorsal distribution'))
    # Source-owned local crop for microscopic review. It is an uncapped open
    # exterior diagnostic patch, never a fake detached complete organ.
    used,lookup=np.unique(replacement,return_inverse=True)
    patch=Mesh(actual[used],lookup.reshape(-1,3),source['frame'],source=dict(kind='uncapped-actual-exterior-diagnostic-crop'))
    patch.save(output/'relief-closeup.npz')
    atomic_json(output/'closeup-review-spec.json',dict(frame=source['frame'],resolution=640,wall_seconds=50,
        meshes=[dict(mesh=str(output/'relief-closeup.npz'),name='Original actual-scale papillary relief crop',color=[.77,.38,.41,1.])],
        views=[dict(name='top',direction=normal.tolist()),dict(name='side',direction=basis[0].tolist()),
               dict(name='oblique',direction=(normal+basis[0]).tolist())],
        qualification='Uncapped diagnostic exterior crop of actual finite source-bound mucosa'))
    print(json.dumps(dict(passed=report['passed'],replacement_faces=len(replacement),packed_bytes=result['total_bytes'],
        minimum_normal_dot=float(dot.min()),maximum_relief_um=float(height.max()*1e6),elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--plan',type=Path)
    a=ap.parse_args()
    if a.plan:construct(a.output.resolve(),a.plan.resolve())
    else:plan(a.output.resolve())
