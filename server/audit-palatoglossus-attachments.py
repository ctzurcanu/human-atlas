"""Independent saved-array palatoglossus attachment audit; no geometry edits."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from tongue_whole_input import WORK,new_output,whole_input


def cap_proof(mesh,face_ids,host,output,name):
    cap=Mesh(mesh.vertices,mesh.faces[face_ids],mesh.frame)
    result,contacts=pair_contacts(cap,host,0,max_pairs=100000,max_tests=1000000)
    atomic_npz(output/(name+'-boundary-witnesses.npz'),**result)
    tri=cap.triangles
    samples=np.r_[tri.reshape(-1,3),tri.mean(1),(tri[:,0]+tri[:,1])/2,(tri[:,1]+tri[:,2])/2,(tri[:,2]+tri[:,0])/2]
    winding=np.abs(igl.winding_number(host.vertices,host.faces,samples))
    distances=np.sqrt(igl.point_mesh_squared_distance(samples,host.vertices,host.faces)[0])
    contained=bool(contacts['complete_within_pair_budget']and contacts['crossing_face_pairs']==0 and
        contacts['coplanar_overlap_face_pairs']==0 and contacts['touching_face_pairs']==0 and np.all(winding>.99))
    return dict(passed=contained,scope='Complete finite cap triangle/host-surface contact test plus connected-component interior witnesses',
        cap_faces=face_ids.tolist(),contacts=contacts,minimum_absolute_winding=float(winding.min()),
        sampled_outside=int(np.sum(winding<.5)),minimum_sampled_host_surface_distance_m=float(distances.min()),
        samples_count=len(samples),witnesses=str(output/(name+'-boundary-witnesses.npz')))


def contact_components(face_ids,mesh):
    todo=set(map(int,face_ids));parts=[];faces=mesh.faces
    # Normal chart vertices can be duplicated. Connectivity must use literal
    # physical controls rather than split shading indices.
    _,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True);physical=inverse[faces]
    vertices={i:set(map(int,physical[i]))for i in todo}
    while todo:
        seed=min(todo);todo.remove(seed);group={seed};controls=set(vertices[seed]);changed=True
        while changed:
            connected={i for i in todo if vertices[i]&controls};changed=bool(connected)
            for i in connected:controls.update(vertices[i])
            group.update(connected);todo.difference_update(connected)
        parts.append(sorted(group))
    return parts


def main(stage,output):
    began=time.monotonic();output=new_output(output);output.mkdir(parents=True)
    candidate=json.loads((stage/'candidate.json').read_text());plan=json.loads((stage/'plan.json').read_text())
    if candidate['plan']['sha256']!=digest_file(stage/'plan.json'):raise ValueError('Candidate plan hash changed')
    binding=plan['whole_input'];tongue,_,actual=whole_input(binding['whole_mesh'],binding['fit_result_path'])
    if binding!=actual:raise ValueError('Whole source evidence changed')
    palate=Mesh.load(plan['soft_palate']['path'])
    if digest_file(plan['soft_palate']['path'])!=plan['soft_palate']['sha256']:raise ValueError('Palate source changed')
    source_bindings=[plan['source_catalogue'],plan['source_audit'],plan['whole'],plan['soft_palate'],
        *plan['fixed_bones'],*plan['fixed_adjoining']]
    if not all(digest_file(r['path'])==r['sha256']for r in source_bindings):raise ValueError('Fixed source file changed')
    sides={};sections=plan['sections'];radial=plan['ring_controls'];body_faces=(sections-1)*radial*2
    for label,row in candidate['sides'].items():
        if digest_file(row['mesh'])!=row['sha256']:raise ValueError('Candidate NPZ changed')
        mesh=Mesh.load(row['mesh'])
        if mesh.geometry_sha!=row['geometry_sha256']:raise ValueError('Candidate geometry changed')
        proximal=np.arange(body_faces,len(mesh.faces),2,dtype=int)
        distal=np.arange(body_faces+1,len(mesh.faces),2,dtype=int)
        tri=mesh.triangles
        rings=[]
        for station in range(sections-1):
            rings.append(np.array([tri[2*(station*radial+k),0]for k in range(radial)]))
        rings.append(tri[distal,1]);rings=np.asarray(rings)
        if not np.array_equal(rings[0],tri[proximal,2]):raise ValueError('Actual first body ring differs from cap ring')
        whole_cap=cap_proof(mesh,proximal,palate,output,label+'-palatal-cap')
        tongue_cap=cap_proof(mesh,distal,tongue,output,label+'-lingual-cap')
        ring_winding={};full_contacts={}
        for key,host in [('palate',palate),('tongue',tongue)]:
            winding=np.abs(igl.winding_number(host.vertices,host.faces,rings.reshape(-1,3))).reshape(sections,radial)
            states=winding>.99;counts=states.sum(1)
            entries=[];exits=[]
            for k in range(radial):
                state=states[:,k];entries.append(int(np.sum(~state[:-1]&state[1:])));exits.append(int(np.sum(state[:-1]&~state[1:])))
            ring_winding[key]=dict(inside_counts_by_station=counts.tolist(),per_radial_path_entries=entries,
                per_radial_path_exits=exits,minimum_first_ring_absolute_winding=float(winding[0].min()),
                minimum_last_ring_absolute_winding=float(winding[-1].min()))
            result,contacts=pair_contacts(mesh,host,0,max_pairs=100000,max_tests=1000000)
            ids=np.flatnonzero((result['relation']=='crossing')|(result['relation']=='coplanar-overlap'))
            muscle_faces=np.unique(result['pairs'][ids,0]);host_faces=np.unique(result['pairs'][ids,1])
            components=contact_components(muscle_faces,mesh)
            atomic_npz(output/(label+'-'+key+'-surface-contact-witnesses.npz'),**result)
            full_contacts[key]=dict(summary=contacts,muscle_faces=muscle_faces.tolist(),host_faces=host_faces.tolist(),
                muscle_contact_components=components,contact_components_count=len(components))
        widths=np.linalg.norm(rings-np.roll(rings,-1,axis=1),axis=2)
        corners=mesh.attributes['normals'][mesh.faces].astype(float);corners/=np.linalg.norm(corners,axis=2)[:,:,None]
        geo=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);geo/=np.linalg.norm(geo,axis=1)[:,None]
        minimum_dot=float(np.einsum('fki,fi->fk',corners,geo).min())
        resources=dict(actual_normal_chart_vertices=len(mesh.vertices),physical_unique_positions=len(np.unique(mesh.vertices,axis=0)),
            normal_chart_position_bytes=12*len(mesh.vertices),normal_chart_normal_bytes=6*len(mesh.vertices),
            actual_index_bytes=mesh.faces.size*4,total_runtime_geometry_bytes=len(mesh.vertices)*18+mesh.faces.size*4,
            source_report_estimated_packed_bytes=row['estimated_packed_bytes'],additional_images=0,UV_bytes=0)
        resources['estimate_matches_actual_expanded_charts']=resources['total_runtime_geometry_bytes']==row['estimated_packed_bytes']
        findings=[]
        if not whole_cap['passed']:findings.append('Complete palatal cap is not wholly contained in the finite palate despite fan control witnesses')
        if not tongue_cap['passed']:findings.append('Complete distal cap is not wholly contained in the actual tongue')
        for key in ['palate','tongue']:
            if full_contacts[key]['contact_components_count']!=1:findings.append(f'{key} crossing footprint has{full_contacts[key]["contact_components_count"]} disconnected candidate-surface components')
        sides[label]=dict(mesh=row['mesh'],mesh_sha256=digest_file(row['mesh']),geometry_sha256=mesh.geometry_sha,
            palatal_cap=whole_cap,lingual_cap=tongue_cap,ring_occupancy=ring_winding,
            full_host_surface_contacts=full_contacts,resources=resources,minimum_actual_corner_normal_face_dot=minimum_dot,
            actual_topology=topology(mesh),substantive_findings=findings,
            endpoint_cap_center_position_matches_plan=dict(
                palate=bool(np.array_equal(tri[proximal[0],0],(np.asarray(row['centerline_m'][0])-np.array([0,1.625,0])).astype('<f4').astype(float)+np.array([0,1.625,0]))),
                tongue=bool(np.array_equal(tri[distal[0],0],(np.asarray(row['centerline_m'][-1])-np.array([0,1.625,0])).astype('<f4').astype(float)+np.array([0,1.625,0])))))
    report=dict(schema='palatoglossus-independent-full-attachment-audit/v1',candidate_sha256=digest_file(stage/'candidate.json'),
        plan_sha256=digest_file(stage/'plan.json'),whole_input=binding,sides=sides,
        source_bindings_unchanged=True,geometry_modified=False,blender_retried=False,anatomical_acceptance=False,
        primary_reference_url='https://pubmed.ncbi.nlm.nih.gov/281278/',
        reference_scope='Adult dissection/radiography/histology supports fan/belly/taper; source-specific dimensions and fibre boundary/attachments remain modeling estimates',
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'audit.json',report)
    print(json.dumps(dict(sides={k:dict(findings=r['substantive_findings'],palatal_cap=r['palatal_cap']['passed'],
        lingual_cap=r['lingual_cap']['passed'],palate_contact_components=r['full_host_surface_contacts']['palate']['contact_components_count'],
        tongue_contact_components=r['full_host_surface_contacts']['tongue']['contact_components_count'],
        runtime_bytes=r['resources']['total_runtime_geometry_bytes'])for k,r in sides.items()},elapsed_seconds=report['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=Path,default=WORK/'palatoglossus-v3');p.add_argument('--output',required=True,type=Path)
    args=p.parse_args();main(args.stage.resolve(),args.output)
