"""Source-bound muscle extensions; adventitial integration is a separate gate."""
from collections import defaultdict
from pathlib import Path
import heapq,json
import numpy as np
from scipy.interpolate import CubicSpline
from .algorithms.solid_attachment import bridge_patches
from .algorithms.surface_normals import grouped_normal_words
from .audit import topology
from .bladder_apical import patch_centroid
from .common import atomic_json,atomic_npz,digest_file
from .model import Mesh
from .spatial import BVH
from .tissue import pack


def contact_patch(mesh, allowed, target, minimum_area_m2, *, maximum_radius_m=.008):
    """Grow a connected source-facet disk; never fabricate a solid opening."""
    allowed=np.asarray(allowed,int);target=np.asarray(target,float)
    sub=Mesh(mesh.vertices,mesh.faces[allowed],mesh.frame)
    query,_=BVH(sub).closest_points(target[None])
    seed=int(allowed[query['face_ids'][0]])
    triangles=mesh.triangles;centers=triangles.mean(1)
    areas=np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)/2
    edges=defaultdict(list)
    for i in allowed:
        for k in range(3):edges[tuple(sorted([tuple(triangles[i,k]),tuple(triangles[i,(k+1)%3])]))].append(int(i))
    neighbors=defaultdict(set)
    for faces in edges.values():
        for i in faces:neighbors[i].update(set(faces)-{i})
    chosen=[seed];front=set(neighbors[seed]);area=areas[seed]
    while area<minimum_area_m2:
        added=False
        for _,i in sorted((float(np.linalg.norm(centers[i]-query['points'][0])),i)for i in front):
            if np.linalg.norm(centers[i]-query['points'][0])>maximum_radius_m:continue
            top=topology(Mesh(mesh.vertices,mesh.faces[chosen+[i]],mesh.frame))
            if len(top['components'])!=1 or top['components'][0]['euler']!=1 or top['invalid_vertex_links']:continue
            chosen.append(i);front.remove(i);front.update(neighbors[i]-set(chosen));area+=areas[i];added=True;break
        if not added:raise ValueError('No connected source disk reaches the explicit attachment area')
    return sorted(chosen),dict(query_point=query['points'][0].tolist(),target_distance_m=float(query['distance_m'][0]),source_patch_area_m2=float(area))


def plan(source_directory, output):
    source_directory,output=Path(source_directory).resolve(),Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    manifest=source_directory/'sources.json';sources=json.loads(manifest.read_text())
    if digest_file(sources['source_catalogue'])!=sources['source_catalogue_sha256']:raise ValueError('Muscular source catalogue changed')
    rows={r['part_id']:r for r in sources['sources']};first_id='BLADDER:male:neck:neck-external-longitudinal'
    first=Mesh.load(rows[first_id]['path']);record=rows[first_id]['record']
    allowed=np.array([i for g in record['tissueDrawGroups']if g.get('neighbor')=='BLADDER:male:neck:adventitia'
        for i in range(g['start']//3,(g['start']+g['count'])//3)],int)
    normals=np.cross(first.triangles[:,1]-first.triangles[:,0],first.triangles[:,2]-first.triangles[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
    specs=[]
    for name,target,second_id,second_target,area,radii,origin,term in [
        ('pubovesical-left',[.006,.894,-.003],'ZA:Hip bone.l:Bone-1',[.012,.898,.033],1e-12,[.0035,.0015],[0,.9140625,0],'A08.3.01.020'),
        ('pubovesical-right',[-.006,.894,-.003],'ZA:Hip bone.r:Bone-1',[-.012,.898,.033],1e-12,[.0035,.0015],[0,.9140625,0],'A08.3.01.020'),
        ('vesicoprostatic',[0,.892,-.023],'ZA:Prostate',[0,.890,-.027],8e-6,[.003,.001],[0,.9140625,-.0244140625],'A08.3.01.021M')]:
        region=allowed[normals[allowed,2]<-.2]if name=='vesicoprostatic'else allowed[normals[allowed,2]>.2]
        fs,info=contact_patch(first,region,target,area)
        second=Mesh.load(rows[second_id]['path']);sn=np.cross(second.triangles[:,1]-second.triangles[:,0],second.triangles[:,2]-second.triangles[:,0]);sn/=np.linalg.norm(sn,axis=1)[:,None]
        region=np.flatnonzero(sn[:,1]>.3)if name=='vesicoprostatic'else np.flatnonzero(sn[:,2]<-.2)
        ss,second_info=contact_patch(second,region,second_target,1e-12)
        specs.append(dict(name=name,ta98_term=term,first_part_id=first_id,first_faces=fs,second_part_id=second_id,second_faces=ss,
            section_radii_m=radii,geometry_origin=origin,first_patch=info,second_patch=second_info,
            source_adventitia_must_be_repartitioned=True,source_bone_is_fixed=name!='vesicoprostatic'))
    output.mkdir();result=dict(schema='bladder-muscular-extension-plan/v1',source_manifest=str(manifest),source_manifest_sha256=digest_file(manifest),
        source_catalogue_sha256=sources['source_catalogue_sha256'],extensions=specs,maximum_step_m=.001,section_controls=16,
        anatomical_acceptance=False,installed=False)
    atomic_json(output/'plan.json',result);return result


def construct(plan_path, output):
    plan_path,output=Path(plan_path).resolve(),Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    plan=json.loads(plan_path.read_text());manifest=Path(plan['source_manifest'])
    if digest_file(manifest)!=plan['source_manifest_sha256']:raise ValueError('Source manifest changed')
    source=json.loads(manifest.read_text());rows={r['part_id']:r for r in source['sources']}
    if digest_file(source['source_catalogue'])!=plan['source_catalogue_sha256']:raise ValueError('Current bladder changed')
    output.mkdir();reports=[]
    for row in plan['extensions']:
        parents=[];patches=[]
        for label in ['first','second']:
            r=rows[row[label+'_part_id']]
            if digest_file(r['path'])!=r['sha256']:raise ValueError('A fixed source changed')
            m=Mesh.load(r['path']);parents.append(m);patches.append(Mesh(m.vertices,m.faces[row[label+'_faces']],m.frame,m.attributes))
        first,second=patches;c0,c1=map(patch_centroid,patches)
        def normal(m):
            n=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0]).sum(0);return n/np.linalg.norm(n)
        length=np.linalg.norm(c1-c0);spline=CubicSpline([0.,1.],np.array([c0,c1]),bc_type=((1,normal(first)*length),(1,-normal(second)*length)))
        centers=spline(np.linspace(0,1,max(8,int(np.ceil(length/plan['maximum_step_m']))+1)))
        radii=np.tile(row['section_radii_m'],(len(centers),1))
        p,surfaces,evidence=bridge_patches(first,second,centers,radii,first_side=-1,second_side=-1,section_controls=plan['section_controls'],projected_phase_alignment=False)
        words,_=grouped_normal_words(p,surfaces,{k:np.full(len(v),k,dtype=object)for k,v in surfaces.items()})
        for label,sid in [('first','first_contact'),('second','second_contact')]:
            m=parents[0 if label=='first'else 1];f=m.faces[row[label+'_faces']]
            w=m.attributes['normals'][f].astype(int)*m.attributes.get('normalSigns',np.ones(len(m.vertices),int))[f,None]
            words[sid]=(-w[:,[0,2,1]]).astype('<i2')
        folder=output/row['name'];pid='BLADDER:male:'+row['name'];cell=[dict(surface=s,side=1)for s in surfaces]
        result=pack(p,surfaces,{pid:cell},{pid:pid},source['frame'],folder,origin=row['geometry_origin'],corner_words=words,max_pairs=1000000)
        actual=Mesh.load(folder/(pid.replace(':','_')+'.npz'))
        first_literal=np.array_equal(actual.triangles[:len(first.faces)],first.triangles[:,[0,2,1]])
        second_literal=np.array_equal(actual.triangles[-len(second.faces):],second.triangles[:,[0,2,1]])
        report=dict(**row,part_id=pid,packed=result,bridge=evidence,first_source_cap_literal=first_literal,second_source_cap_literal=second_literal,
            passed=bool(result['passed']and first_literal and second_literal),adventitial_repartition_required=True,
            complete_body_packing_verified=False,anatomical_acceptance=False,installed=False)
        reports.append(report);atomic_json(folder/'candidate.json',report)
    result=dict(plan=str(plan_path),plan_sha256=digest_file(plan_path),extensions=reports,passed=all(r['passed']for r in reports),
        anatomical_acceptance=False,installed=False,adventitial_integration_complete=False)
    atomic_json(output/'candidate.json',result);return result
