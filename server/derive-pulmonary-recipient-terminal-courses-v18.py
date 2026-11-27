"""Source-bound recipient terminal course graph from actual HRA airway surfaces.

Scalar-level contours and their far-side component lineage supply course
initialization charts. They are not vascular lumina, points drilled through an
organ, or a replacement anatomical hierarchy. No geometry is moved.
"""
import json
import time
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components,dijkstra
from geometry_tools.common import digest_file,identity
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def contour_components(v,f,distance,level):
    tri=distance[f];ids=np.flatnonzero((tri.min(1)<level)&(tri.max(1)>level));nodes={};points=[];source_edges=[];fractions=[];segments=[];face_ids=[]
    unresolved=[]
    for fi in ids:
        face=f[fi];cross=[]
        for a,b in zip(face,np.roll(face,-1)):
            if (distance[a]-level)*(distance[b]-level)<0.:
                key=tuple(sorted((int(a),int(b))))
                if key not in nodes:
                    a,b=key;t=(level-distance[a])/(distance[b]-distance[a]);nodes[key]=len(points)
                    points.append(v[a]+t*(v[b]-v[a]));source_edges.append(key);fractions.append(t)
                cross.append(nodes[key])
        if len(cross)==2:segments.append(cross);face_ids.append(int(fi))
        else:unresolved.append(dict(face=int(fi),reason='level passes a literal source vertex'))
    if not segments:return [],unresolved
    points=np.asarray(points);segments=np.asarray(segments);source_edges=np.asarray(source_edges);fractions=np.asarray(fractions)
    graph=coo_matrix((np.ones(len(segments)),(segments[:,0],segments[:,1])),shape=(len(points),len(points))).tocsr()
    count,labels=connected_components(graph,directed=False)
    degree=np.bincount(segments.ravel(),minlength=len(points));rows=[]
    for ci in range(count):
        pp=np.flatnonzero(labels==ci);ss=np.flatnonzero(labels[segments[:,0]]==ci)
        if not len(ss):continue
        if not np.all(degree[pp]==2):
            unresolved.append(dict(component=int(ci),source_faces=np.asarray(face_ids)[ss].tolist(),reason='contour is not closed degree two',degree_counts=np.bincount(degree[pp]).tolist()));continue
        local=np.full(len(points),-1,int);local[pp]=np.arange(len(pp))
        rows.append(dict(points=points[pp],source_edges=source_edges[pp],fractions=fractions[pp],
          segments=local[segments[ss]],faces=np.asarray(face_ids)[ss]))
    return rows,unresolved


def main():
    began=time.monotonic();family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue'])
    assert digest_file(path)==family['source_catalogue_sha256'];atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    interfaces=json.loads((BASE/'source-correspondence-contract-v7/contract.json').read_text())['airway_interface_landmarks']
    old_charts=np.load(BASE/'source-correspondence-contract-v7/literal-fixed-airway-interface-controls.npz',allow_pickle=False)
    out=BASE/'recipient-terminal-course-graph-v18';out.mkdir(exist_ok=True);arrays={};vertices=[];faces=[];owners=[];sources=[]
    for i,pid in enumerate(family['fixed_current_airway_ids']):
        mesh=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres');offset=len(vertices)
        vertices.extend(mesh.vertices);faces.extend(mesh.faces+offset);owners.extend([i]*len(mesh.faces))
        sources.append(dict(id=pid,geometry_sha256=mesh.geometry_sha,source_part_index=i,faces=len(mesh.faces)))
    v,inv=np.unique(np.asarray(vertices),axis=0,return_inverse=True);f=inv[np.asarray(faces)];owners=np.asarray(owners)
    arrays.update(literal_world_vertices=v,literal_world_faces=f,face_owner_source_index=owners)
    edges=np.unique(np.sort(np.vstack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0)
    length=np.linalg.norm(v[edges[:,0]]-v[edges[:,1]],axis=1);keep=length>0.;edges,length=edges[keep],length[keep]
    graph=coo_matrix((np.r_[length,length],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(v),len(v))).tocsr()
    # The separately stored main-bronchus charts have 32 shared edges each;
    # derive the free boundary of the actual combined surface, rather than
    # assuming these two overlapping charts are two separate root mouths.
    loops=boundary_loops(f)
    maincharts=[r for r in interfaces if 'main_bronchus'in r['part_id']and not r['all_edges_have_literal_opposed_other_branch']]
    represented=set(tuple(p)for r in maincharts for p in old_charts[r['array_key']])
    root_loops=[loop for loop in loops if all(tuple(p)in represented for p in v[loop])]
    assert len(root_loops)==1,'Need one actual proximal boundary identified by both named main-entry source charts'
    distal_loops=[loop for loop in loops if not all(tuple(p)in represented for p in v[loop])]
    rootcharts=[];seeds=[]
    for i,loop in enumerate(root_loops):
        pp=v[loop];key=f'literal_combined_main_entry_{i}';arrays[key]=pp;seeds.extend(map(int,loop))
        assert all(tuple(p)in represented for p in pp)
        rootcharts.append(dict(array_key=key,points=len(loop),source_main_charts=maincharts,
          actual_combined_surface_free_boundary=True,independent_left_right_root_mouths_not_assumed=True))
    distance=dijkstra(graph,directed=False,indices=np.asarray(seeds),min_only=True)
    assert np.isfinite(distance).all(),'Not every current airway surface vertex is reachable from the two actual main-entry charts'
    arrays['literal_edge_graph_distance_from_main_entries_m']=distance
    span=.001371;levels=np.arange(.000713,float(distance.max())+span,span);nodes=[];links=[];unresolved=[];previous_labels=None;previous_nodes={};vertex_label=np.full(len(v),-1,int)
    members_at_node={};node_by_key={}
    for li,level in enumerate(levels):
        active=np.flatnonzero(distance>level)
        if not len(active):break
        count,labels=connected_components(graph[active][:,active],directed=False)
        vertex_label.fill(-1);vertex_label[active]=labels
        contours,bad=contour_components(v,f,distance,float(level))
        if bad:unresolved.append(dict(level_m=float(level),contours=bad))
        candidate=defaultdict(list)
        for row in contours:
            high=np.where(distance[row['source_edges'][:,0]]>distance[row['source_edges'][:,1]],row['source_edges'][:,0],row['source_edges'][:,1])
            components=np.unique(vertex_label[high])
            if len(components)!=1 or components[0]<0:
                unresolved.append(dict(level_m=float(level),reason='closed contour has more than one far-side component'));continue
            candidate[int(components[0])].append(row)
        current_nodes={}
        for component,rows in sorted(candidate.items()):
            if len(rows)!=1:
                unresolved.append(dict(level_m=float(level),component=component,reason='more than one contour bounds one far-side component; no artificial merge',contours=len(rows)));continue
            row=rows[0];key=f'level{li}_component{component}';member=active[labels==component];members_at_node[key]=member
            arrays[key+'_profile_points']=row['points'];arrays[key+'_literal_source_edges']=row['source_edges'];arrays[key+'_edge_fractions']=row['fractions'];arrays[key+'_profile_segments']=row['segments'];arrays[key+'_source_face_ids']=row['faces']
            unique,nn=np.unique(owners[row['faces']],return_counts=True)
            owner_counts=[dict(source_part_id=sources[int(o)]['id'],crossed_faces=int(n))for o,n in zip(unique,nn)]
            owner_counts.sort(key=lambda x:-x['crossed_faces'])
            node=dict(key=key,level_m=float(level),points_array=key+'_profile_points',literal_source_edges_array=key+'_literal_source_edges',edge_fractions_array=key+'_edge_fractions',
              source_faces_array=key+'_source_face_ids',profile_segments_array=key+'_profile_segments',profile_points=len(row['points']),
              chart_centre_m=row['points'].mean(0).tolist(),actual_surface_owner_contributions=owner_counts,
              far_side_source_vertex_count=len(member),closed_degree_two_contour=True,initialization_chart_not_vascular_wall=True)
            nodes.append(node);node_by_key[key]=node;current_nodes[component]=key
            if previous_labels is not None:
                parent_ids=np.unique(previous_labels[member]);parent_ids=parent_ids[parent_ids>=0]
                if len(parent_ids)==1 and int(parent_ids[0])in previous_nodes:
                    links.append(dict(parent=previous_nodes[int(parent_ids[0])],child=key,nested_actual_far_side_vertex_components=True))
                elif len(parent_ids)>1:
                    unresolved.append(dict(level_m=float(level),reason='far-side component does not descend from one previous component',component=component))
        previous_labels=vertex_label.copy();previous_nodes=current_nodes
    children=defaultdict(list)
    for row in links:children[row['parent']].append(row['child'])
    terminals=[];membership={k:set(map(int,a))for k,a in members_at_node.items()}
    edge_faces=defaultdict(list)
    for fi,face in enumerate(f):
        for a,b in zip(face,np.roll(face,-1)):edge_faces[tuple(sorted((int(a),int(b))))].append(fi)
    for ti,loop in enumerate(distal_loops):
        key=f'literal_terminal_boundary_{ti}';arrays[key+'_profile']=v[loop];arrays[key+'_source_vertex_ids']=np.asarray(loop,int)
        source_faces=np.array([edge_faces[tuple(sorted((int(a),int(b))))][0]for a,b in zip(loop,np.roll(loop,-1))],int)
        assert all(len(edge_faces[tuple(sorted((int(a),int(b))))])==1 for a,b in zip(loop,np.roll(loop,-1)))
        arrays[key+'_source_face_ids']=source_faces
        unique,nn=np.unique(owners[source_faces],return_counts=True)
        counts=[dict(source_part_id=sources[int(o)]['id'],boundary_edges=int(n))for o,n in zip(unique,nn)]
        required=set(map(int,loop));minimum=float(distance[loop].min());maximum=float(distance[loop].max())
        possible=[row for row in nodes if row['level_m']<minimum and required<=membership[row['key']]]
        parent=max(possible,key=lambda row:row['level_m'])if possible else None
        terminals.append(dict(key=key,actual_terminal_boundary_profile_array=key+'_profile',source_face_ids_array=key+'_source_face_ids',
          source_vertex_ids_array=key+'_source_vertex_ids',profile_points=len(loop),
          minimum_root_edge_distance_m=minimum,maximum_root_edge_distance_m=maximum,
          literal_profile_centre_m=v[loop].mean(0).tolist(),actual_surface_owner_contributions=counts,
          final_closed_course_chart=parent['key']if parent else None,
          final_chart_to_terminal_root_distance_gap_m=minimum-parent['level_m']if parent else None,
          actual_boundary_edges_each_have_one_literal_owner_face=True,
          anatomical_vessel_identity_not_assigned_by_mouth_centre=True))
    forks=[dict(node=k,children=vv)for k,vv in children.items()if len(vv)>1]
    np.savez_compressed(out/'literal-recipient-contour-lineage-and-terminal-witnesses.npz',**arrays)
    result=dict(schema='current-fixed-HRA-airway-recipient-terminal-course-initialization-graph/v18',
      source_catalogue_sha256=family['source_catalogue_sha256'],ordinary_runtime_evidence_sha256=digest_file(BASE/'actual-runtime-context-v9/runtime-context.json'),
      literal_surface_geometry_sha256=identity('fixed-HRA-airway-surface-composite/v1',v,f),
      source_records=sources,root_actual_main_entry_charts=rootcharts,
      actual_world_vertices=len(v),actual_faces=len(f),actual_free_boundary_loops=len(loops),actual_distal_boundary_loops=len(distal_loops),distance_level_spacing_m=span,
      maximum_actual_edge_graph_distance_m=float(distance.max()),nodes=nodes,links=links,forks=forks,terminal_witnesses=terminals,
      unresolved_contours=unresolved,arrays_sha256=digest_file(out/'literal-recipient-contour-lineage-and-terminal-witnesses.npz'),
      method='Weighted literal triangle-edge shortest distance, actual closed level contours, and nested far-side vertex-component lineage. Not an exact continuous medial axis or exact geodesic.',
      source_vessel_to_recipient_branch_assignment_pending=True,
      no_averaged_branch_forks=True,terminal_courses_not_only_near_entry_controls=True,
      no_anatomical_or_vascular_patency_acceptance=True,geometry_changed=False,catalogue_changed=False,
      elapsed_seconds=time.monotonic()-began)
    (out/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),course_charts=len(nodes),lineage_links=len(links),forks=len(forks),
      actual_terminal_boundary_witnesses=len(terminals),unresolved_levels=len(unresolved),
      evidence_sha256=digest_file(out/'evidence.json'),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
