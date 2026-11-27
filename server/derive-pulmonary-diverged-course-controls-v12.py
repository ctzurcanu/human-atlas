"""Named course controls after donor-common takeoffs, using literal surface charts."""
import json
import time
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra,connected_components

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.tube_station_topology import recover_station_rings

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def intrinsic_profile(vertices,faces,seed_points,level):
    lookup={tuple(p):i for i,p in enumerate(vertices)}
    seeds=np.array([lookup[tuple(p)]for p in seed_points],int)
    edges=np.unique(np.sort(np.vstack([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
    length=np.linalg.norm(vertices[edges[:,0]]-vertices[edges[:,1]],axis=1)
    graph=coo_matrix((np.r_[length,length],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(vertices),len(vertices))).tocsr()
    distance=dijkstra(graph,directed=False,indices=seeds,min_only=True)
    segments=[];nodes={};points=[]
    for face in faces:
        vals=distance[face]
        if not vals.min()<level<vals.max():continue
        crossing=[]
        for a,b in zip(face,np.roll(face,-1)):
            if (distance[a]-level)*(distance[b]-level)<0:
                key=tuple(sorted((int(a),int(b))))
                if key not in nodes:
                    lo,hi=key;t=(level-distance[lo])/(distance[hi]-distance[lo]);nodes[key]=len(points)
                    points.append(vertices[lo]+t*(vertices[hi]-vertices[lo]))
                crossing.append(nodes[key])
        if len(crossing)!=2:return None,dict(reason='vertex-valued or non-simple scalar contour')
        segments.append(crossing)
    if not segments:return None,dict(reason='no contour')
    segments=np.asarray(segments,int);degree=np.bincount(segments.ravel(),minlength=len(points))
    if not np.all(degree==2):return None,dict(reason='contour not closed degree-two',degree_histogram=np.bincount(degree).tolist())
    count=connected_components(coo_matrix((np.ones(len(segments)),(segments[:,0],segments[:,1])),shape=(len(points),len(points))),directed=False,return_labels=False)
    if count!=1:return None,dict(reason='more than one course contour',components=int(count))
    return np.asarray(points),dict(controls=len(points),components=1,every_degree_two=True,
                                  intrinsic_coordinate='shortest weighted literal triangle-edge graph distance from actual entry boundary; not exact continuous geodesic')


def main():
    began=time.monotonic();prior=json.loads((BASE/'takeoff-and-root-controls-v10/evidence.json').read_text());family=json.loads((BASE/'family-contract-v2.json').read_text())
    path=Path(family['source_catalogue']);assert digest_file(path)==family['source_catalogue_sha256'];atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    inputs=np.load(BASE/'takeoff-and-root-controls-v10/literal-correspondence-and-root-profiles.npz',allow_pickle=False)
    out=BASE/'diverged-course-controls-v12';out.mkdir();arrays={};controls=[];unresolved=[]
    for ci,row in enumerate(prior['named_control_pairs']):
        # Separate target HRA B I/B II and tied source S4 origins need a real
        # branch topology correspondence, not an averaged artificial point.
        if len(row['target_ids'])!=1 or not row['source_unique_largest_entry_profile']:
            unresolved.append(dict(source_id=row['source_id'],reason='Combined or tied source/target branch topology; no averaged course control emitted'));continue
        source_mesh=reader.read_part(atlas,parts[row['source_id']],frame='model:male:assembled-metres')
        sv,si=np.unique(source_mesh.vertices,axis=0,return_inverse=True);sf=si[source_mesh.faces]
        tube=recover_station_rings(sv,sf)[row['source_entry_component']];sf=sf[tube['face_ids']]
        target_id=row['target_ids'][0];target_mesh=reader.read_part(atlas,parts[target_id],frame=source_mesh.frame)
        tv,ti=np.unique(target_mesh.vertices,axis=0,return_inverse=True);tf=ti[target_mesh.faces]
        source_seed=inputs[row['source_entry_array']];target_seed=inputs[row['target_array']]
        trials=[];chosen=None
        for level in [.003,.005,.008,.012]:
            sp,sm=intrinsic_profile(sv,sf,source_seed,level);tp,tm=intrinsic_profile(tv,tf,target_seed,level)
            trials.append(dict(edge_graph_distance_m=level,source=sm,target=tm))
            if sp is not None and tp is not None:
                chosen=(level,sp,tp,sm,tm);break
        if chosen is None:
            unresolved.append(dict(source_id=row['source_id'],target_id=target_id,reason='No homologous unique closed source/target intrinsic course chart',trials=trials));continue
        level,sp,tp,sm,tm=chosen;sk=f'case{ci}_source_course_profile';tk=f'case{ci}_target_course_profile';arrays[sk]=sp;arrays[tk]=tp
        controls.append(dict(source_id=row['source_id'],target_id=target_id,kind='named-after-divergence-course-chart',
                             source_array=sk,target_array=tk,edge_graph_distance_m=level,
                             source_point_m=sp.mean(0).tolist(),target_point_m=tp.mean(0).tolist(),
                             source_geometry_sha256=source_mesh.geometry_sha,target_geometry_sha256=target_mesh.geometry_sha,
                             exact_branch_registration_claimed=False,source_profile=sm,target_profile=tm,
                             initialization_only=True,trials=trials))
    # A common donor entry maps to the actual named HRA common-parent entry,
    # not to several contradictory downstream child ports.
    interface=json.loads((BASE/'source-correspondence-contract-v7/contract.json').read_text())['airway_interface_landmarks']
    charts=np.load(BASE/'source-correspondence-contract-v7/literal-fixed-airway-interface-controls.npz',allow_pickle=False)
    common=[]
    for group,parent,grand in [
        ([0,1,2],'right_superior_lobar_bronchus','right_main_bronchus'),
        ([3,4],'right_middle_lobar_bronchus','right_inferior_lobar_bronchus'),
        ([10,11],'left_superior_lobar_bronchus','left_main_bronchus'),
    ]:
        source_rows=[prior['named_control_pairs'][i]for i in group];source_points=np.array([r['source_point_m']for r in source_rows])
        pid='LOCAL:male:pulmonary:bronchus:'+parent;gpid='LOCAL:male:pulmonary:bronchus:'+grand
        possible=[r for r in interface if r['part_id']==pid and gpid in r['other_branch_ids'] and r['all_edges_have_literal_opposed_other_branch']]
        if len(possible)!=1:
            unresolved.append(dict(source_ids=[r['source_id']for r in source_rows],target_id=pid,reason='No unique complete actual common-parent entry chart'));continue
        target_points=charts[possible[0]['array_key']];key='common_'+parent;arrays[key+'_target_entry_chart']=target_points;arrays[key+'_source_node_representations']=source_points
        common.append(dict(kind='named-common-parent-entry',source_ids=[r['source_id']for r in source_rows],target_id=pid,
                           source_point_m=source_points.mean(0).tolist(),target_point_m=target_points.mean(0).tolist(),
                           maximum_source_common_node_representation_difference_m=float(np.linalg.norm(source_points[:,None]-source_points,axis=2).max()),
                           target_chart_array=key+'_target_entry_chart',source_points_array=key+'_source_node_representations',
                           target_not_an_average_of_separate_child_entries=True,initialization_only=True))
    np.savez_compressed(out/'literal-intrinsic-course-charts.npz',**arrays)
    result=dict(schema='pulmonary-named-diverged-course-initialization-controls/v12',source_catalogue_sha256=family['source_catalogue_sha256'],
                arrays_sha256=digest_file(out/'literal-intrinsic-course-charts.npz'),course_controls=controls,common_parent_controls=common,
                unresolved_correspondences=unresolved,
                graph_incidence_preserved_by_shared_common_parent_nodes=True,
                duplicate_native_entries_not_sent_to_different_target_child_ports=True,
                no_named_midcourse_point_created_from_whole_mesh_centroid=True,
                all_points_are_chart_centres_not_anatomical_walls_or_holes=True,
                new_source_initialization_not_historical_HRA_field_replay=True,
                source_field_not_yet_constructed=True,geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,
                elapsed_seconds=time.monotonic()-began)
    (out/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),course_controls=len(controls),common_parent_controls=len(common),unresolved=len(unresolved),
                          evidence_sha256=digest_file(out/'evidence.json'),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
