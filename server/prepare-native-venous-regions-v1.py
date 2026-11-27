"""Literal native venous-region retrieval, without moving or constructing tissue."""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-native-venous-regions-v1'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'


def main():
    began=time.monotonic();assert digest_file(SOURCE)==EXPECTED;assert not BASE.exists();BASE.mkdir()
    atlas=json.loads(SOURCE.read_bytes());parts={p['id']:p for p in atlas['parts']}
    census_path=ROOT/'reports/ta98-blood-circulation-inventory-installed-v11.json'
    census=json.loads(census_path.read_bytes());terms={r['id']:r for r in census['terms']}
    body_path=ROOT/'.local-models/organ-models/spinal-cord/male-body.json';body=json.loads(body_path.read_bytes());bp={p['id']:p for p in body['parts']}
    reader=AtlasReader();frame='model:local-ta98:assembled-metres';meshes={};sources={}
    def load(pid):
        if pid in meshes:return meshes[pid]
        p=parts[pid];assert not p.get('suppressed')
        assert dict(bp[pid],chunk=body['sourceChunkIndices'][bp[pid]['chunk']])==p
        m=reader.read_part(atlas,p,frame=frame,catalogue_path=SOURCE,catalogue_sha=EXPECTED)
        meshes[pid]=m
        sources[pid]=dict(record=p,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,
                          bank_sha256=m.source['bank_sha256'],source_vertices=m.source['source_vertex_ids'])
        return m
    rows=[]
    for side,suffix in [('left','l'),('right','r')]:
        radius=load(f'ZA:Radius.{suffix}:Bone');annular=load(f'ZA:Annular ligament of radius.{suffix}')
        # The annular ligament identifies the radial head. The distal radius
        # endpoint fixes the wrist plane; anatomical boundary remains qualified.
        distal=radius.vertices[radius.vertices[:,1]<np.quantile(radius.vertices[:,1],.025)].mean(0)
        proximal=annular.vertices.mean(0);axis=proximal-distal;length=np.linalg.norm(axis);axis/=length
        for base,term,parent in [('Cephalic','A12.3.08.021','A12.3.08.015'),('Basilic','A12.3.08.022','A12.3.08.018')]:
            pid=f'ZA:{base} vein.{suffix}';p=parts[pid];m=load(pid)
            assert p['license']=='CC BY-SA 4.0' and p['provenance']['url']=='https://github.com/Z-Anatomy/Models-of-human-anatomy'
            assert not p.get('sourceOffset') and not p.get('geometryInstance') and p.get('geometryScale',1)==1
            assert pid in terms[parent]['active_exact_source_part_ids'] and terms[term]['representation_status']=='missing-direct-representation'
            components=recover_station_rings(m.vertices,m.faces);assert len(components)==1
            rings,replay=align_station_correspondence(m.vertices,m.faces,components[0]['rings'])
            rings=np.asarray(rings);centers=m.vertices[rings].mean(1);projection=(centers-distal)@axis
            selected=np.flatnonzero((projection>=0)&(projection<=length));assert len(selected)>2
            first,last=int(selected[0]),int(selected[-1]);assert np.array_equal(selected,np.arange(first,last+1));assert first>0 and last<len(rings)-1
            verts=np.unique(rings[first:last+1]);mask=np.isin(m.faces,verts).all(1);faces=np.flatnonzero(mask)
            local=Mesh(m.vertices,m.faces[faces],frame);loops=boundary_loops(local.faces);assert len(loops)==2
            assert set(np.concatenate(loops))==set(rings[[first,last]].reshape(-1))
            row=dict(term=term,name=terms[term]['name_en'],authority_parent=terms[term]['parent_id'],side=side,source_part_id=pid,source_exact_term=parent,
                     source_station_count=len(rings),selected_station_interval=[first,last],source_face_ids=faces.tolist(),source_region_faces=len(faces),
                     region_bounds_m=local.bounds.tolist(),source_region_length_m=float(np.linalg.norm(np.diff(centers[first:last+1],axis=0),axis=1).sum()),
                     anatomical_landmarks=dict(proximal_part=annular.source['record']['id'],distal_part=radius.source['record']['id'],proximal_m=proximal.tolist(),distal_m=distal.tolist(),axis=axis.tolist()),
                     boundary_status='Provisional local radial-head and distal-radius planes; literal complete source rings, no vertex move or new cap.',
                     shared_interfaces=[dict(source_ring_station=i,source_vertex_ids=rings[i].tolist(),center_m=centers[i].tolist(),adjacent_source='remaining original '+pid) for i in [first,last]],
                     source_outer_wall_only=True,native_registration='Same original baked Z-Anatomy catalogue and actual body record; no offset, instance or scale.',lumen_patency_unverified=True)
            atomic_npz(BASE/(term+'-'+side+'.npz'),source_face_ids=faces,rings=rings,source_vertices=m.vertices,source_faces=m.faces)
            rows.append(row)
    pid='ZA:Azygos vein';m=load(pid);svc=load('ZA:Superior vena cava');bronchus=load('LOCAL:male:pulmonary:bronchus:right_main_bronchus')
    assert pid in terms['A12.3.07.001']['active_exact_source_part_ids'] and terms['A12.3.07.002']['representation_status']=='missing-direct-representation'
    comps=recover_station_rings(m.vertices,m.faces);assert len(comps)==1
    rings,replay=align_station_correspondence(m.vertices,m.faces,comps[0]['rings']);rings=np.asarray(rings);centers=m.vertices[rings].mean(1)
    tree=BVH(svc);enddist=[tree.closest_point(p)['distance_m'] for p in centers[[0,-1]]]
    if enddist[0]<enddist[1]:rings=rings[::-1];centers=centers[::-1]
    tangent=np.gradient(centers,axis=0);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    upper=np.flatnonzero((centers[:,1]>np.mean(bronchus.bounds[:,1]))&(tangent[:,2]>.25))
    assert len(upper)>2
    first=int(upper[0]);assert first>0
    ids=np.unique(rings[first:]);mask=np.isin(m.faces,ids).all(1);faces=np.flatnonzero(mask)
    local=Mesh(m.vertices,m.faces[faces],frame);assert len(boundary_loops(local.faces))==2
    arrays,crossings=pair_contacts(local,bronchus,0.,max_pairs=30000,max_tests=700000)
    atomic_npz(BASE/'azygos-bronchus-contact.npz',**arrays)
    endpoint,_=tree.closest_points(m.vertices[rings[-1]])
    rows.append(dict(term='A12.3.07.002',name='arch of azygos vein',authority_parent='A12.3.07.001',source_part_id=pid,source_exact_term='A12.3.07.001',side='right',
                     source_station_count=len(rings),selected_station_interval=[first,len(rings)-1],source_face_ids=faces.tolist(),source_region_faces=len(faces),region_bounds_m=local.bounds.tolist(),
                     source_region_length_m=float(np.linalg.norm(np.diff(centers[first:],axis=0),axis=1).sum()),
                     boundary_status='Provisional start of anterior-directed terminal bend, referenced to current right main bronchus; retains literal source rings.',
                     shared_interfaces=[dict(source_ring_station=first,source_vertex_ids=rings[first].tolist(),center_m=centers[first].tolist(),adjacent_source='remaining original azygos vein'),
                                        dict(source_ring_station=len(rings)-1,source_vertex_ids=rings[-1].tolist(),center_m=centers[-1].tolist(),adjacent_source='ZA:Superior vena cava',nearest_surface_distance_m=endpoint['distance_m'].tolist(),blood_communication_unverified=True)],
                     current_right_bronchus_contact=crossings,source_outer_wall_only=True,native_registration='Same baked native Z-Anatomy catalogue; actual current right bronchus checked.',lumen_patency_unverified=True))
    atomic_npz(BASE/'A12.3.07.002-right.npz',source_face_ids=faces,rings=rings,source_vertices=m.vertices,source_faces=m.faces)
    # Save the actual source arrays once for the subsequent literal bank pack.
    for pid,mesh in meshes.items():
        filename=hashlib.sha256(pid.encode()).hexdigest()[:12]+'.npz';mesh.save(BASE/filename);sources[pid]['mesh']=str((BASE/filename).relative_to(ROOT))
        sources[pid]['source_vertices']=np.asarray(sources[pid]['source_vertices']).tolist()
    report=dict(schema='native-venous-region-source-plan/v1',source_catalogue=str(SOURCE.relative_to(ROOT)),source_catalogue_sha256=EXPECTED,
                actual_body=str(body_path.relative_to(ROOT)),actual_body_sha256=digest_file(body_path),inventory=str(census_path.relative_to(ROOT)),inventory_sha256=digest_file(census_path),
                regions=rows,sources=sources,unique_new_terms=sorted({r['term'] for r in rows}),
                source_reference_only_not_completed_finite_vessels=True,original_source_face_partition_required=True,
                prospective_controls='Whole source concepts must select all disjoint physical parts; each new region selects its own literal faces. Systems whole/component reconciliation must be verified before installation.',
                primary_references=['https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/Entity%20TA98%20EN/12.3.08.021%20Entity%20TA98%20EN.htm','https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/Entity%20TA98%20EN/12.3.08.022%20Entity%20TA98%20EN.htm','https://pubmed.ncbi.nlm.nih.gov/9358971/'],
                unverified_private_cerebral_network_not_used=True,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(BASE/'source-plan.json',report)
    assert digest_file(SOURCE)==EXPECTED
    print(json.dumps(dict(output=str(BASE),terms=report['unique_new_terms'],regions=[{k:r[k] for k in ['term','side','source_part_id','selected_station_interval','source_region_faces','source_region_length_m']} for r in rows],azygos_bronchus_crossing_pairs=rows[-1]['current_right_bronchus_contact']['crossing_face_pairs'],elapsed_seconds=report['elapsed_seconds'])))


if __name__=='__main__':main()
