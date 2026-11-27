"""Qualify unsigned wall distances with source-chart and volume incidence.

Temporary disks close mathematical outer/inner classifiers only. They are not
anatomical caps, branch openings, render geometry or an installable repair.
"""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import topology,winding_contains
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.algorithms.annular_pipe_envelopes import close_with_shared_disks

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/posterior-intercostal-native-v1'
OUT=B/'root-incidence-v1';C=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME='model:local-ta98:assembled-metres';ROUNDING_BOUND_M=2e-7

def literal_submesh(mesh,face_ids):
    faces=mesh.faces[face_ids];ids=np.unique(faces);mapping=np.full(len(mesh.vertices),-1,int);mapping[ids]=np.arange(len(ids))
    attrs={k:v[ids]for k,v in mesh.attributes.items()}
    return Mesh(mesh.vertices[ids],mapping[faces],mesh.frame,attrs,dict(source_geometry_sha256=mesh.geometry_sha,source_attributes_sha256=mesh.attributes_sha,source_face_ids=np.asarray(face_ids).tolist(),coordinates_and_attributes_unmodified=True))

def compact_exact(mesh):
    v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True)
    return Mesh(v,inv[mesh.faces],mesh.frame,source=dict(exact_coordinate_topological_chart=True,source_surface_sha256=mesh.surface_sha))

def nearest_polyline(points,centres):
    points=np.asarray(points);a,b=centres[:-1],centres[1:];d=b-a;den=np.einsum('ij,ij->i',d,d)
    t=np.clip(np.einsum('nsi,si->ns',points[:,None]-a[None],d)/den[None],0,1)
    q=a[None]+t[:,:,None]*d[None];length=np.linalg.norm(points[:,None]-q,axis=2);segments=np.argmin(length,axis=1)
    return dict(distance_m=length[np.arange(len(points)),segments].tolist(),segment_ids=segments.tolist(),segment_fraction=t[np.arange(len(points)),segments].tolist(),closest_points_m=q[np.arange(len(points)),segments].tolist(),maximum_station_center_residual_m=float(np.linalg.norm(points[:,None]-centres[None],axis=2).min(1).max()))

def run():
    started=time.monotonic();OUT.mkdir(exist_ok=False,parents=True);a=json.loads(C.read_bytes());assert digest_file(C)=='46262d8ab07682fcae3ff0ee634f652abad766c1e010f6e2ec311755aea9abf5'
    parts={p['id']:p for p in a['parts']};reader=AtlasReader();inventory=json.loads((B/'inventory.json').read_bytes());bindings={p['id']:p for p in inventory['native_sources']}
    parent_specs=[('supreme-l','ZA:Supreme intercostal artery.l',None,'Outer source envelope; finite inner blood boundary unavailable'),('supreme-r','ZA:Supreme intercostal artery.r',None,'Outer source envelope; finite inner blood boundary unavailable'),('thoracic-outer','TA98WALL:male:thoracic-aorta:externa',0,'Actual installed outer chart; not the blood-lumen boundary'),('thoracic-inner','TA98WALL:male:thoracic-aorta:intima',1,'Actual installed intimal inner chart; declared blood-space boundary'),('arch-outer','ZA:Aortic arch',None,'Outer source envelope; finite inner blood boundary unavailable'),('abdominal-outer','ZA:Abdominal aorta',None,'Outer source envelope; finite inner blood boundary unavailable')]
    classifiers={};parent_report=[]
    for key,pid,group,scope in parent_specs:
        full=reader.read_part(a,parts[pid],frame=FRAME)
        if group is None:chart=full;face_ids=np.arange(len(full.faces))
        else:
            g=parts[pid]['tissueDrawGroups'][group];face_ids=np.arange(g['start']//3,(g['start']+g['count'])//3);chart=literal_submesh(full,face_ids)
        logical=compact_exact(chart);recovered=recover_station_rings(logical.vertices,logical.faces);assert len(recovered)==1
        centres=logical.vertices[np.array(recovered[0]['rings'])].mean(1)
        v,f,disks=close_with_shared_disks(logical.vertices,logical.faces,expected_mouths=2);closed=Mesh(v,f,FRAME,source=dict(temporary_classifier_only=True,anatomical_disks=False,installed_parent_geometry_sha256=full.geometry_sha,source_chart_geometry_sha256=chart.geometry_sha,computational_cap_faces_start=len(logical.faces)))
        topo=topology(closed);assert topo['closed_oriented_manifold']and len(topo['components'])==1
        closed.save(OUT/(key+'-TEMPORARY-classifier.npz'));atomic_npz(OUT/(key+'-source-chart-centres.npz'),centres=centres,frame=np.array(FRAME),source_face_ids=face_ids)
        classifiers[key]=(chart,closed,centres,BVH(chart),BVH(closed))
        parent_report.append(dict(key=key,id=pid,installed_geometry_sha256=full.geometry_sha,installed_attributes_sha256=full.attributes_sha,chart_geometry_sha256=chart.geometry_sha,chart_surface_sha256=chart.surface_sha,source_face_ids=face_ids.tolist(),source_components=1,source_stations=len(centres),source_ring_size=len(recovered[0]['rings'][0]),temporary_classifier_geometry_sha256=closed.geometry_sha,computational_disk_faces=len(disks),computational_disks_not_rendered=True,scope=scope))
    rows=[]
    for side in ['l','r']:
        old=json.loads((B/f'{side}-physical-ports-v2.json').read_bytes())
        for course in old['courses']:
            source=reader.read_part(a,parts[course['source_id']],frame=FRAME);binding=bindings[course['source_id']]
            assert source.geometry_sha==binding['geometry_sha256']and source.attributes_sha==binding['attributes_sha256']
            original=json.loads((B/f'{side}-courses.json').read_bytes())['courses'][course['course_index']]
            ends=np.array([q['center_m']for q in original['source_mouths']])
            keys=['supreme-'+side]if course['level']<=2 else['thoracic-outer','thoracic-inner','arch-outer','abdominal-outer']
            declared=keys[0];line=nearest_polyline(ends,classifiers[declared][2]);root_end=int(np.argmin(line['distance_m']));root=original['source_mouths'][root_end]
            points=np.vstack([np.array(root['center_m']),np.array(root['world_points_m'])]);incidence=[]
            for key in keys:
                chart,closed,centres,index,closed_index=classifiers[key]
                nearest,near_evidence=index.closest_points(points,backend='bvh');mathematical,math_evidence=closed_index.closest_points(points,backend='bvh')
                winding=winding_contains(closed,points,boundary_tolerance_m=ROUNDING_BOUND_M,distance_backend='bvh');inside=np.array(winding['inside']);ambiguous=np.array(winding['boundary_ambiguous']);on_disk=(mathematical['face_ids']>=len(chart.faces))&(mathematical['distance_m']<=ROUNDING_BOUND_M)
                labels=['at-temporary-source-mouth-closure'if cap else'boundary-ambiguous'if amb else'inside-parent-envelope'if ins else'outside-parent-envelope'for cap,amb,ins in zip(on_disk,ambiguous,inside)]
                signed=[None if amb or cap else float(-d if ins else d)for d,ins,amb,cap in zip(nearest['distance_m'],inside,ambiguous,on_disk)]
                source_line=nearest_polyline(points[:1],centres)
                incidence.append(dict(parent_chart=key,point_order=['source-root-centre']+[f'source-rim-{k:02}'for k in range(len(points)-1)],points_m=points.tolist(),classification=labels,signed_distance_to_actual_envelope_surface_m=signed,unsigned_distance_to_actual_envelope_surface_m=nearest['distance_m'].tolist(),nearest_actual_chart_face_ids=nearest['face_ids'].tolist(),winding=winding,temporary_cap_ambiguity=on_disk.tolist(),source_centreline_incidence=source_line,centreline_incident_within_float32_rounding_bound=source_line['distance_m'][0]<=ROUNDING_BOUND_M,physiological_lumen_join_inferred=False))
            rows.append(dict(side=side,level=course['level'],source_id=course['source_id'],source_component=course['source_component'],declared_parent_id=course['expected_parent'],source_root_endpoint=root_end,previous_unsigned_distance_root_endpoint=course['root_endpoint'],declared_parent_source_centreline_endpoint_queries=line,previous_unsigned_rim_surface_distances_m=course['parent_endpoint_queries'][course['root_endpoint']],incidence=incidence,origin_not_relabelled=True,geometry_unchanged=True))
            atomic_json(OUT/'checkpoint.json',dict(complete=False,parents=parent_report,courses=rows,elapsed_seconds=time.monotonic()-started),overwrite=(OUT/'checkpoint.json').exists())
    report=dict(schema='posterior-intercostal-parent-envelope-incidence/v1',source_catalogue_sha256=digest_file(C),complete=True,parent_charts=parent_report,courses=rows,rounding_boundary_m=ROUNDING_BOUND_M,scope='Unsigned surface distance is not a gap. Inside-parent or source-mouth incidence can reflect a parent radius; actual finite-wall root ostia/shared interfaces remain unverified.',no_origin_relabel=True,no_aortic_or_skeletal_movement=True,prior_gap_interpretations_superseded=[str(B/f'{s}-physical-ports-v2.json')for s in ['l','r']],temporary_cap_disks_are_not_anatomy=True,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-started)
    atomic_json(OUT/'proof.json',report)
    print(json.dumps(dict(courses=len(rows),declared_parent_centreline_incident=sum(q['incidence'][0]['centreline_incident_within_float32_rounding_bound']for q in rows),classified_centres={key:sum(q['classification'][0]in['inside-parent-envelope','at-temporary-source-mouth-closure']for r in rows for q in r['incidence']if q['parent_chart']==key)for key in classifiers},elapsed_seconds=report['elapsed_seconds'])))

if __name__=='__main__':run()
