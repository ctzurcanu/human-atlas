"""One terminal-bound, rooted right S2 arterial course construction candidate.

This uses the source-qualified recurrent and ascending source fragments. It
maps the original bronchial reference graph to source-bound current HRA course
lineage and transports the source bronchovascular offsets, then rotates every
literal arterial section properly. Existing arterial origins remain fixed.
No finite communicating wall or anatomical installation is claimed here.
"""
import argparse
import json
import time
import itertools
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.optimize import linear_sum_assignment
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader,Mesh

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OUT=BASE/'right-S2-terminal-bound-construction-v19'
SOURCE_GUIDE='ZA:Posterior segmental bronchus of right lung (BII)'
TARGET_GUIDE='LOCAL:male:pulmonary:bronchus:right_posterior_bronchus'
RECURRENT='ZA:Apical segmental artery of right lung'
ASCENDING='ZA:Posterior segmental artery of right lung'


def graph_of_curves(rows,arrays,sections_key='sections_array'):
    points=[];curve_points=[];radii=[];offsets=[]
    for row in rows:
        ss=arrays[row[sections_key]];cc=ss.mean(1);offsets.append(len(points));points.extend(cc)
        radii.extend(np.linalg.norm(ss-cc[:,None],axis=2).max(1));curve_points.append(cc)
    points=np.asarray(points);parent=np.arange(len(points))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return int(i)
    for a,b in sorted(cKDTree(points).query_pairs(1e-6)):
        aa,bb=find(a),find(b)
        if aa!=bb:parent[bb]=aa
    labels=np.array([find(i)for i in range(len(points))]);unique=np.unique(labels);remap={int(x):i for i,x in enumerate(unique)}
    node=np.array([remap[int(x)]for x in labels]);coords=np.array([points[labels==x].mean(0)for x in unique])
    radius=np.array([np.asarray(radii)[labels==x].max()for x in unique]);edges=[];curve_nodes=[]
    for start,cc in zip(offsets,curve_points):
        ids=node[start:start+len(cc)];curve_nodes.append(ids);edges.extend(tuple(sorted((int(a),int(b))))for a,b in zip(ids,ids[1:])if a!=b)
    edges=np.array(sorted(set(edges)),int);length=np.linalg.norm(coords[edges[:,0]]-coords[edges[:,1]],axis=1)
    graph=coo_matrix((np.r_[length,length],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(coords),len(coords))).tocsr()
    return coords,graph,curve_nodes,radius


def read_family():
    family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue'])
    assert digest_file(path)==family['source_catalogue_sha256'];atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']}
    return family,atlas,parts,AtlasReader()


def ancestry(key,parent):
    result=[]
    while key is not None:result.append(key);key=parent.get(key)
    return result


def contract():
    began=time.monotonic();family,atlas,parts,reader=read_family();OUT.mkdir()
    source=json.loads((BASE/'rooted-source-recipient-graph-v17/evidence.json').read_text())
    src=np.load(BASE/'rooted-source-recipient-graph-v17/literal-source-nodes-and-recipient-course-charts.npz',allow_pickle=False)
    guides=[r for r in source['original_hidden_bronchial_reference_graph']['source_curves']if r['part_id']==SOURCE_GUIDE]
    coords,graph,curve_nodes,radius=graph_of_curves(guides,src)
    assert graph.nnz//2==len(coords)-1,'The named source reference must be an actual single centreline tree'
    old_roots=json.loads((BASE/'takeoff-and-root-controls-v10/evidence.json').read_text());entry=next(r for r in old_roots['named_control_pairs']if r['source_id']==SOURCE_GUIDE)
    root=int(np.linalg.norm(coords-np.array(entry['source_point_m']),axis=1).argmin());assert np.linalg.norm(coords[root]-entry['source_point_m'])<1e-6
    distance,pred=dijkstra(graph,directed=False,indices=root,return_predecessors=True);assert np.isfinite(distance).all()
    degree=np.diff(graph.indptr);leaves=np.flatnonzero((degree==1)&(np.arange(len(coords))!=root));critical=set(np.flatnonzero(degree!=2))|{root}
    children=defaultdict(list)
    for i,p in enumerate(pred):
        if p>=0:children[int(p)].append(i)
    desc={}
    def descend(i):
        dd=([i]if i in leaves else[])
        for child in children[i]:dd+=descend(child)
        desc[i]=dd;return dd
    descend(root)
    recipient=json.loads((BASE/'recipient-terminal-course-graph-v18/evidence.json').read_text())
    rr=np.load(BASE/'recipient-terminal-course-graph-v18/literal-recipient-contour-lineage-and-terminal-witnesses.npz',allow_pickle=False)
    nodes={r['key']:r for r in recipient['nodes']};parents={r['child']:r['parent']for r in recipient['links']}
    terminals=[r for r in recipient['terminal_witnesses']if all(x['source_part_id']==TARGET_GUIDE for x in r['actual_surface_owner_contributions'])]
    assert len(terminals)==6
    for row in terminals:assert ancestry(row['final_closed_course_chart'],parents)[-1]=='level0_component0'
    all_lineages=[]
    for row in recipient['terminal_witnesses']:
        trace=ancestry(row['final_closed_course_chart'],parents)if row['final_closed_course_chart']else[]
        all_lineages.append(dict(terminal=row['key'],trace=trace,continuous_to_actual_proximal_prefix=bool(trace and trace[-1]=='level0_component0'),
          disconnected_component_root=trace[-1]if trace else None))
    common=json.loads((BASE/'diverged-course-controls-v12/evidence.json').read_text())
    common=next(r for r in common['common_parent_controls']if r['target_id'].endswith('right_superior_lobar_bronchus'))
    assert SOURCE_GUIDE in common['source_ids']
    common_arrays=np.load(BASE/'diverged-course-controls-v12/literal-intrinsic-course-charts.npz',allow_pickle=False)
    common_chart=common_arrays[common['target_chart_array']];target_root=np.array(common['target_point_m'])
    target_world=rr['literal_world_vertices'];lookup={tuple(p):i for i,p in enumerate(target_world)}
    root_depth=float(np.median(rr['literal_edge_graph_distance_from_main_entries_m'][[lookup[tuple(p)]for p in common_chart]]))
    target_paths={};target_coords={'shared-superior-lobar-entry':target_root};target_parent={}
    for row in terminals:
        reverse=ancestry(row['final_closed_course_chart'],parents)[::-1]
        selected=[k for k in reverse if nodes[k]['level_m']>root_depth]
        assert selected
        last='shared-superior-lobar-entry';path=[last]
        for k in selected:
            point=np.array(nodes[k]['chart_centre_m']);target_coords[k]=point
            if k in target_parent:assert target_parent[k]==last
            else:target_parent[k]=last
            last=k;path.append(k)
        key=row['key'];target_coords[key]=np.array(row['literal_profile_centre_m']);target_parent[key]=last;path.append(key);target_paths[key]=path
    # Branch-end matching is an explicit geometric initialization within the
    # already identified S2 territory. It does not name subsegmental arteries.
    a=coords[leaves]-coords[root];a/=np.linalg.norm(a,axis=1)[:,None]
    b=np.array([target_coords[r['key']]-target_root for r in terminals]);b/=np.linalg.norm(b,axis=1)[:,None]
    cost=1.-a@b.T;ii,jj=linear_sum_assignment(cost)
    assert len(ii)==len(leaves),'Not enough actual S2 terminal witnesses for the source reference leaves'
    leaf_target={int(leaves[i]):terminals[j]['key']for i,j in zip(ii,jj)}
    target_node={root:'shared-superior-lobar-entry'}
    for i in sorted(critical,key=lambda x:distance[x]):
        if i==root:continue
        if i in leaf_target:target_node[i]=leaf_target[i];continue
        traces=[ancestry(leaf_target[leaf],target_parent)[::-1]for leaf in desc[i]]
        shared=[k for k in traces[0]if all(k in trace for trace in traces[1:])]
        assert shared;target_node[i]=shared[-1]
    arrays=dict(source_guide_nodes=coords,source_guide_distance_m=distance,source_guide_root=np.array(root),target_actual_common_entry=common_chart)
    edges=[];collapsed=[]
    for child in sorted(critical,key=lambda x:distance[x]):
        if child==root:continue
        path=[child];p=int(pred[child])
        while p not in critical:path.append(p);p=int(pred[p])
        path.append(p);path=path[::-1];left,right=target_node[p],target_node[child]
        reverse=ancestry(right,target_parent)
        if left not in reverse:raise AssertionError('Target branch mapping lost parent ancestry')
        route=reverse[:reverse.index(left)+1][::-1]
        points=np.array([target_coords[k]for k in route]);source_points=coords[path]
        key=f'guide_edge{len(edges)}';arrays[key+'_source_polyline']=source_points;arrays[key+'_target_polyline']=points
        length=float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum())
        if length<1e-6:collapsed.append(dict(source_parent=int(p),source_child=int(child),target_parent=left,target_child=right))
        edges.append(dict(key=key,source_parent=int(p),source_child=int(child),source_node_path=list(map(int,path)),
          target_node_path=route,source_polyline_array=key+'_source_polyline',target_polyline_array=key+'_target_polyline',
          source_length_m=float(distance[child]-distance[p]),target_course_length_m=length,
          named_parent_entry_precedes_child_divergence=True))
    arterial=json.loads((BASE/'native-curves-v3/evidence.json').read_text());art=np.load(BASE/'native-curves-v3/literal-section-centres.npz',allow_pickle=False)
    selected=[r for r in arterial['curves']if (r['part_id']==RECURRENT and r['component']in [0,1,2,3])or r['part_id']==ASCENDING]
    assert len(selected)==7
    source_audit=json.loads((BASE/'evidence.json').read_text())
    actual_components=np.load(BASE/'literal-components.npz',allow_pickle=False)
    roots=[];selected_rows=[]
    for row in selected:
        key=row['key'];arrays[key+'_source_sections']=art[key+'_literal_sections'];arrays[key+'_source_face_ids']=art[key+'_source_face_ids']
        selected_rows.append(dict(row,source_sections_array=key+'_source_sections',source_faces_array=key+'_source_face_ids'))
        if row['component']==1:
            end=0 if row['larger_radius_endpoint']==0 else-1;arrays[key+'_fixed_origin_profile']=art[key+'_literal_sections'][end]
            roots.append(dict(part_id=row['part_id'],component=1,endpoint=end,profile_array=key+'_fixed_origin_profile',
              TA98_construction_identity='A12.2.01.109'if row['part_id']==RECURRENT else'A12.2.01.108',
              parent_part=RECURRENT if row['part_id']==RECURRENT else'ZA:Right pulmonary artery',
              parent_source_component=9 if row['part_id']==RECURRENT else 0,
              parent_interface='Actual upper-trunk takeoff station47; side-port surface ownership still requires finite junction construction'if row['part_id']==RECURRENT else'Actual source-distal RPA open 12-point rim; phase/coordinate pairing still needs replay',
              fixed_origin_prefix_m=.004,transition_length_m=.016,communicating_wall_unestablished=True))
    np.savez_compressed(OUT/'literal-guides-target-routes-and-protected-origins.npz',**arrays)
    result=dict(schema='right-S2-two-origin-terminal-bound-course-construction-contract/v19',
      source_catalogue_sha256=family['source_catalogue_sha256'],source_catalogue=family['source_catalogue'],
      source_semantic_root_review_sha256=digest_file(BASE/'right-posterior-takeoff-review-v6/root-source-correspondence-review-v1.json'),
      current_recipient_terminal_evidence_sha256=digest_file(BASE/'recipient-terminal-course-graph-v18/evidence.json'),
      arrays_sha256=digest_file(OUT/'literal-guides-target-routes-and-protected-origins.npz'),
      selected_source_curves=selected_rows,qualified_arterial_source_identities=[dict(term='A12.2.01.109',source=RECURRENT,components=[0,1,2,3]),dict(term='A12.2.01.108',source=ASCENDING,components=[0,1,2])],
      original_bronchial_guide=SOURCE_GUIDE,fixed_current_bronchial_guide=TARGET_GUIDE,
      source_reference_tree=dict(nodes=len(coords),edges=graph.nnz//2,root=root,critical_nodes=sorted(map(int,critical)),leaf_nodes=leaves.tolist()),
      source_reference_leaf_to_actual_recipient_terminal=[dict(source_leaf=int(i),target_terminal=k,geometric_initialization_only=True)for i,k in leaf_target.items()],
      actual_target_common_parent=common,source_and_target_critical_node_pairs=[dict(source_node=int(i),target_node=k)for i,k in target_node.items()],
      mapped_source_guide_edges=edges,collapsed_guide_edge_mappings=collapsed,
      all_108_terminal_root_lineage_qualification=all_lineages,
      selected_terminal_mouths=terminals,
      native_arterial_origins_and_parent_ownership=roots,
      moving_family='Only seven actual right-S2 arterial source components; both full source parts retain every other component literally.',
      fixed_dependencies='All other52 qualified vascular family parts, all current heart/skeleton/lungs/airways/context, both origin profiles, upper-trunk component9, RPA root/recipient, and nine coupled lingular pieces.',
      bronchovascular_offset_transport='Original arterial centre minus nearest original BII reference course point is transported by one proper local rotation; artery centres are never put on the bronchial guide.',
      construction_rules=['No collapsed source reference edge may be used for construction.',
        'Every actual shared source vascular centre node uses one mapping and one rooted protection weight.',
        'Retain the actual literal arterial ring shape, area and local radii by proper rotations.',
        'Current HRA bronchial envelope is a fixed blocker; actual triangle checks must follow construction.',
        'Current root interface face ownership and finite communicating walls remain unaccepted until reconstructed and verified.'],
      construction_ready=not collapsed,geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,
      elapsed_seconds=time.monotonic()-began)
    (OUT/'contract.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),source_reference_leaves=len(leaves),actual_S2_terminals=len(terminals),
      source_guide_critical_edges=len(edges),collapsed_edges=len(collapsed),construction_ready=result['construction_ready'],
      contract_sha256=digest_file(OUT/'contract.json'),elapsed_seconds=result['elapsed_seconds'])))


def topology():
    began=time.monotonic();c=json.loads((OUT/'contract.json').read_text());a=np.load(OUT/'literal-guides-target-routes-and-protected-origins.npz',allow_pickle=False)
    arrays={k:a[k]for k in a.files};source=a['source_guide_nodes'];root=c['source_reference_tree']['root'];leaves=c['source_reference_tree']['leaf_nodes']
    original=json.loads((BASE/'recipient-terminal-course-graph-v18/evidence.json').read_text());nodes={r['key']:r for r in original['nodes']}
    parent={r['child']:r['parent']for r in original['links']};terminals=c['selected_terminal_mouths'];target_root=np.array(c['actual_target_common_parent']['target_point_m'])
    chart=a['target_actual_common_entry'];rr=np.load(BASE/'recipient-terminal-course-graph-v18/literal-recipient-contour-lineage-and-terminal-witnesses.npz',allow_pickle=False)
    lookup={tuple(p):i for i,p in enumerate(rr['literal_world_vertices'])};depth=float(np.median(rr['literal_edge_graph_distance_from_main_entries_m'][[lookup[tuple(p)]for p in chart]]))
    target_coords={'shared-superior-lobar-entry':target_root};target_parent={};paths={}
    for row in terminals:
        trace=[k for k in ancestry(row['final_closed_course_chart'],parent)[::-1]if nodes[k]['level_m']>depth]
        last='shared-superior-lobar-entry';path=[last]
        for k in trace:
            target_coords[k]=np.array(nodes[k]['chart_centre_m']);target_parent[k]=last;last=k;path.append(k)
        key=row['key'];target_coords[key]=np.array(row['literal_profile_centre_m']);target_parent[key]=last;path.append(key);paths[key]=path
    source_children=defaultdict(list)
    for edge in c['mapped_source_guide_edges']:source_children[edge['source_parent']].append(edge['source_child'])
    desc={}
    def descend(i):
        result=[i]if i in leaves else[]
        for child in source_children[i]:result+=descend(child)
        desc[i]=result;return result
    descend(root)
    aa=source[leaves]-source[root];aa/=np.linalg.norm(aa,axis=1)[:,None]
    bb=np.array([target_coords[r['key']]-target_root for r in terminals]);bb/=np.linalg.norm(bb,axis=1)[:,None];cost=1.-aa@bb.T
    lca={}
    for mask in range(1,1<<len(terminals)):
        traces=[paths[terminals[i]['key']]for i in range(len(terminals))if mask&(1<<i)]
        common=[k for k in traces[0]if all(k in p for p in traces[1:])];lca[mask]=common[-1]
    best=None;valid=0;checked=0
    for order in itertools.permutations(range(len(terminals))):
        checked+=1;assignment={leaf:order[i]for i,leaf in enumerate(leaves)};mapping={root:'shared-superior-lobar-entry'}
        for source_node in c['source_reference_tree']['critical_nodes']:
            if source_node==root:continue
            mapping[source_node]=lca[sum(1<<assignment[leaf]for leaf in desc[source_node])]
        okay=all(mapping[e['source_parent']]!=mapping[e['source_child']]and mapping[e['source_parent']]in ancestry(mapping[e['source_child']],target_parent)for e in c['mapped_source_guide_edges'])
        if not okay:continue
        valid+=1;score=float(sum(cost[i,order[i]]for i in range(len(leaves))))
        if best is None or score<best[0]:best=(score,order,mapping)
    out=BASE/'right-S2-topology-qualified-contract-v20';out.mkdir()
    result=dict(schema='right-S2-finite-terminal-ancestry-assignment-qualification/v20',source_contract_sha256=digest_file(OUT/'contract.json'),
      checked_terminal_bijections=checked,noncollapsed_ancestry_preserving_bijections=valid,
      source_catalogue_sha256=c['source_catalogue_sha256'],geometry_changed=False,catalogue_changed=False,
      construction_ready=best is not None,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    if best is not None:
        score,order,mapping=best;edges=[]
        for i,edge in enumerate(c['mapped_source_guide_edges']):
            left,right=mapping[edge['source_parent']],mapping[edge['source_child']];trace=ancestry(right,target_parent);route=trace[:trace.index(left)+1][::-1]
            points=np.array([target_coords[k]for k in route]);key=f'qualified_guide_edge{i}';arrays[key+'_target_polyline']=points
            edges.append(dict(edge,target_node_path=route,target_polyline_array=key+'_target_polyline',target_course_length_m=float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum())))
        result.update(qualified_source_guide_edges=edges,critical_node_target_mapping=[dict(source_node=i,target_node=k)for i,k in mapping.items()],
          selected_terminal_assignment=[dict(source_leaf=leaf,target_terminal=terminals[order[i]]['key'])for i,leaf in enumerate(leaves)],
          angular_initialization_cost=score,terminal_branch_identity_acceptance=False)
        np.savez_compressed(out/'literal-qualified-guides-and-roots.npz',**arrays);result['arrays_sha256']=digest_file(out/'literal-qualified-guides-and-roots.npz')
    (out/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),checked=checked,valid=valid,construction_ready=result['construction_ready'],
      evidence_sha256=digest_file(out/'evidence.json'),elapsed_seconds=result['elapsed_seconds'])))


def ordered():
    """A guide is allowed to branch along a recipient course before its fork.

    The new arterial geometry will retain its OWN source incidence. This
    guidance correspondence does not claim a bronchial graph isomorphism.
    """
    began=time.monotonic();c=json.loads((OUT/'contract.json').read_text());a=np.load(OUT/'literal-guides-target-routes-and-protected-origins.npz',allow_pickle=False)
    arrays={k:a[k]for k in a.files};recipient=json.loads((BASE/'recipient-terminal-course-graph-v18/evidence.json').read_text())
    nodes={r['key']:r for r in recipient['nodes']};parent={r['child']:r['parent']for r in recipient['links']}
    rr=np.load(BASE/'recipient-terminal-course-graph-v18/literal-recipient-contour-lineage-and-terminal-witnesses.npz',allow_pickle=False)
    lookup={tuple(p):i for i,p in enumerate(rr['literal_world_vertices'])};chart=a['target_actual_common_entry']
    depth=float(np.median(rr['literal_edge_graph_distance_from_main_entries_m'][[lookup[tuple(p)]for p in chart]]))
    coords={'shared-superior-lobar-entry':np.array(c['actual_target_common_parent']['target_point_m'])};target_parent={};paths={};target_depth={'shared-superior-lobar-entry':0.}
    for row in c['selected_terminal_mouths']:
        trace=[k for k in ancestry(row['final_closed_course_chart'],parent)[::-1]if nodes[k]['level_m']>depth]
        last='shared-superior-lobar-entry';path=[last]
        for k in trace:
            coords[k]=np.array(nodes[k]['chart_centre_m']);target_depth[k]=nodes[k]['level_m']-depth;target_parent[k]=last;last=k;path.append(k)
        k=row['key'];coords[k]=np.array(row['literal_profile_centre_m']);target_depth[k]=row['minimum_root_edge_distance_m']-depth;target_parent[k]=last;path.append(k);paths[k]=path
    root=c['source_reference_tree']['root'];leaves=c['source_reference_tree']['leaf_nodes'];sd=a['source_guide_distance_m']
    source_children=defaultdict(list)
    for edge in c['mapped_source_guide_edges']:source_children[edge['source_parent']].append(edge['source_child'])
    desc={}
    def descend(i):
        result=[i]if i in leaves else[]
        for child in source_children[i]:result+=descend(child)
        desc[i]=result;return result
    descend(root)
    assigned={r['source_leaf']:r['target_terminal']for r in c['source_reference_leaf_to_actual_recipient_terminal']}
    mapped=dict(assigned);allocation=[];unresolved=[]
    for i in sorted(c['source_reference_tree']['critical_nodes'],key=lambda j:-sd[j]):
        if i in assigned or i==root:continue
        traces=[paths[assigned[leaf]]for leaf in desc[i]]
        common=[k for k in traces[0]if all(k in trace for trace in traces[1:])]
        candidates=[k for k in common if k!='shared-superior-lobar-entry'and all(k!=mapped[child]and k in ancestry(mapped[child],target_parent)for child in source_children[i])]
        if not candidates:unresolved.append(dict(source_node=i,reason='No strictly ordered actual common-course chart before all mapped children'));continue
        desired=float(sd[i]/min(sd[leaf]for leaf in desc[i])*min(target_depth[assigned[leaf]]for leaf in desc[i]))
        chosen=min(candidates,key=lambda k:(abs(target_depth[k]-desired),target_depth[k]));mapped[i]=chosen
        allocation.append(dict(source_node=i,target_course_chart=chosen,desired_root_distance_m=desired,
          actual_root_distance_m=target_depth[chosen],actual_recipient_fork_homology_claimed=False,
          source_arterial_incidence_is_preserved_separately=True))
    mapped[root]='shared-superior-lobar-entry';edges=[]
    if not unresolved:
        for ei,edge in enumerate(c['mapped_source_guide_edges']):
            left,right=mapped[edge['source_parent']],mapped[edge['source_child']]
            trace=ancestry(right,target_parent)
            assert left in trace and left!=right
            route=trace[:trace.index(left)+1][::-1];points=np.array([coords[k]for k in route]);key=f'ordered_guide_edge{ei}'
            arrays[key+'_target_polyline']=points
            edges.append(dict(edge,target_node_path=route,target_polyline_array=key+'_target_polyline',
              target_course_length_m=float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum()),
              target_path_is_actual_course_chart_lineage=True))
    out=BASE/'right-S2-ordered-chart-contract-v21';out.mkdir();np.savez_compressed(out/'literal-ordered-guides-and-roots.npz',**arrays)
    result=dict(c,schema='right-S2-ordered-actual-course-chart-construction-contract/v21',
      predecessor_contract_sha256=digest_file(OUT/'contract.json'),terminal_fork_bijection_audit_sha256=digest_file(BASE/'right-S2-topology-qualified-contract-v20/evidence.json'),
      guide_fork_allocations=allocation,unresolved_ordered_charts=unresolved,
      ordered_source_guide_edges=edges,critical_node_target_mapping=[dict(source_node=int(i),target_node=k)for i,k in mapped.items()],
      arrays_sha256=digest_file(out/'literal-ordered-guides-and-roots.npz'),
      guide_correspondence_is_a_declared_nonisomorphic_construction=True,
      shared_sibling_guide_course_intervals_do_not_imply_shared_vascular_tissue=True,
      no_claim_that_arterial_and_bronchial_topologies_are_equal=True,
      collapsed_guide_edge_mappings=[],construction_ready=not unresolved,
      geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    (out/'contract.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),edges=len(edges),fork_allocations=len(allocation),unresolved=len(unresolved),
      construction_ready=result['construction_ready'],contract_sha256=digest_file(out/'contract.json'),elapsed_seconds=result['elapsed_seconds'])))


def proper_between(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float);a/=np.maximum(np.linalg.norm(a,axis=1,keepdims=True),1e-100);b/=np.maximum(np.linalg.norm(b,axis=1,keepdims=True),1e-100)
    axis=np.cross(a,b);cos=np.sum(a*b,axis=1);k=np.zeros((len(a),3,3));k[:,0,1]=-axis[:,2];k[:,0,2]=axis[:,1];k[:,1,0]=axis[:,2];k[:,1,2]=-axis[:,0];k[:,2,0]=-axis[:,1];k[:,2,1]=axis[:,0]
    result=np.eye(3)+k+(k@k)/np.maximum(1.+cos,1e-12)[:,None,None]
    for i in np.flatnonzero(cos<-.999999):
        basis=np.eye(3)[int(np.argmin(np.abs(a[i])))];normal=np.cross(a[i],basis);normal/=np.linalg.norm(normal);result[i]=Rotation.from_rotvec(np.pi*normal).as_matrix()
    return result


def blend_rotations(left,right,weight):
    q0=Rotation.from_matrix(left).as_quat();q1=Rotation.from_matrix(right).as_quat();dot=np.sum(q0*q1,axis=1)
    q1=np.where((dot<0)[:,None],-q1,q1);dot=np.clip(np.abs(dot),0.,1.);theta=np.arccos(dot);sine=np.sin(theta)
    w=np.asarray(weight,float);near=sine<1e-8
    a=np.sin((1.-w)*theta)/np.maximum(sine,1e-100);b=np.sin(w*theta)/np.maximum(sine,1e-100)
    a[near]=1.-w[near];b[near]=w[near];q=a[:,None]*q0+b[:,None]*q1;q/=np.linalg.norm(q,axis=1)[:,None]
    return Rotation.from_quat(q).as_matrix()


def smoothstep(x):
    t=np.clip(x,0.,1.);return t*t*t*(10.+t*(-15.+6.*t))


class CourseMap:
    def __init__(self,c,a):
        self.edges=c['ordered_source_guide_edges'];self.splines=[];starts=[];ends=[];edge_index=[];us=[];ue=[];self.endpoint_rotation={}
        for ei,row in enumerate(self.edges):
            old=a[row['source_polyline_array']];new=a[row['target_polyline_array']]
            oldlen=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(old,axis=0),axis=1))];newlen=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(new,axis=0),axis=1))]
            assert oldlen[-1]>1e-6 and newlen[-1]>1e-6
            keep=np.r_[True,np.diff(newlen)>1e-9];new=new[keep];newlen=newlen[keep]
            spline=CubicSpline(newlen/newlen[-1],new,bc_type='natural');self.splines.append(spline)
            starts.extend(old[:-1]);ends.extend(old[1:]);edge_index.extend([ei]*(len(old)-1));us.extend(oldlen[:-1]/oldlen[-1]);ue.extend(oldlen[1:]/oldlen[-1])
            # A shared critical source node has one proper rotation, taken
            # from its actual incoming guide course. The root uses its one
            # outgoing common-parent course and is reused by all children.
            self.endpoint_rotation[row['source_child']]=proper_between((old[-1]-old[-2])[None],spline(1.,1)[None])[0]
        root=c['source_reference_tree']['root'];first=next((i,row)for i,row in enumerate(self.edges)if row['source_parent']==root)
        ei,row=first;old=a[row['source_polyline_array']]
        self.endpoint_rotation[root]=proper_between((old[1]-old[0])[None],self.splines[ei](0.,1)[None])[0]
        self.starts=np.asarray(starts);self.ends=np.asarray(ends);self.direction=self.ends-self.starts;self.length2=np.sum(self.direction**2,axis=1)
        self.ei=np.asarray(edge_index);self.us=np.asarray(us);self.ue=np.asarray(ue)

    def apply(self,x):
        delta=x[:,None]-self.starts[None];fraction=np.clip(np.einsum('nsi,si->ns',delta,self.direction)/self.length2,0.,1.)
        pp=self.starts[None]+fraction[:,:,None]*self.direction[None];dist=np.sum((x[:,None]-pp)**2,axis=2);segment=np.argmin(dist,axis=1)
        index=np.arange(len(x));q=pp[index,segment];u=self.us[segment]+fraction[index,segment]*(self.ue[segment]-self.us[segment])
        edge=self.ei[segment];yy=np.empty_like(x);rr=np.empty((len(x),3,3));base=np.empty_like(x)
        for ei,row in enumerate(self.edges):
            ids=np.flatnonzero(edge==ei)
            if not len(ids):continue
            t=u[ids];spline=self.splines[ei];new=spline(t);r=proper_between(self.direction[segment[ids]],spline(t,1))
            left=np.broadcast_to(self.endpoint_rotation[row['source_parent']],r.shape);right=np.broadcast_to(self.endpoint_rotation[row['source_child']],r.shape)
            r=blend_rotations(left,r,smoothstep(t/.12));r=blend_rotations(r,right,1.-smoothstep((1.-t)/.12))
            yy[ids]=new+np.einsum('nij,nj->ni',r,x[ids]-q[ids]);rr[ids]=r;base[ids]=new
        return yy,rr,q,base


def candidate():
    began=time.monotonic();folder=BASE/'right-S2-ordered-chart-contract-v21';c=json.loads((folder/'contract.json').read_text());assert c['construction_ready']
    a=np.load(folder/'literal-ordered-guides-and-roots.npz',allow_pickle=False);family,atlas,parts,reader=read_family();assert family['source_catalogue_sha256']==c['source_catalogue_sha256']
    mapping=CourseMap(c,a);selected=c['selected_source_curves'];coords,graph,curve_nodes,radius=graph_of_curves(selected,a,sections_key='source_sections_array')
    root_nodes=[]
    for row in c['native_arterial_origins_and_parent_ownership']:
        pp=a[row['profile_array']];cc=pp.mean(0);node=int(np.linalg.norm(coords-cc,axis=1).argmin());assert np.linalg.norm(coords[node]-cc)<1e-6;root_nodes.append(node)
    distance=dijkstra(graph,directed=False,indices=np.asarray(root_nodes),min_only=True);assert np.isfinite(distance).all()
    mapped,rotation,q,base=mapping.apply(coords);weight=smoothstep((distance-.004)/.016)
    target=coords+weight[:,None]*(mapped-coords);rotation=blend_rotations(np.tile(np.eye(3),(len(coords),1,1)),rotation,weight)
    out=BASE/'right-S2-terminal-bound-candidate-v22';out.mkdir();saved={};position_maps=defaultdict(dict);rotation_maps=defaultdict(dict);metrics=[];fragment_paths=[]
    source_components=np.load(BASE/'literal-components.npz',allow_pickle=False);audit=json.loads((BASE/'evidence.json').read_text())
    for i,(row,node_ids)in enumerate(zip(selected,curve_nodes)):
        sections=a[row['source_sections_array']];cc=sections.mean(1);r=rotation[node_ids];displacement=target[node_ids]-coords[node_ids]
        new=cc[:,None]+displacement[:,None]+np.einsum('nij,nkj->nki',r,sections-cc[:,None])
        local=sections-cc[:,None];newlocal=new-new.mean(1)[:,None]
        area0=np.linalg.norm(np.cross(local,np.roll(local,-1,axis=1)).sum(1),axis=1)*.5
        area1=np.linalg.norm(np.cross(newlocal,np.roll(newlocal,-1,axis=1)).sum(1),axis=1)*.5
        profile=np.max(np.abs(np.linalg.norm(new[:,:,None]-new[:,None,:],axis=3)-np.linalg.norm(sections[:,:,None]-sections[:,None,:],axis=3)))
        key=row['key'];saved[key+'_source_sections']=sections;saved[key+'_candidate_sections']=new;saved[key+'_proper_rotation']=r;saved[key+'_root_distance_m']=distance[node_ids]
        for si,(oldring,newring)in enumerate(zip(sections,new)):
            for oldpoint,newpoint in zip(oldring,newring):
                pkey=tuple(oldpoint)
                if pkey in position_maps[row['part_id']]:assert np.linalg.norm(position_maps[row['part_id']][pkey]-newpoint)<1e-12
                position_maps[row['part_id']][pkey]=newpoint;rotation_maps[row['part_id']][pkey]=r[si]
        metrics.append(dict(part_id=row['part_id'],component=row['component'],stations=len(cc),
          maximum_relative_section_area_error=float(np.abs(area1/area0-1.).max()),maximum_literal_profile_pairwise_error_m=float(profile),
          maximum_radius_error_m=float(np.abs(np.linalg.norm(newlocal,axis=2)-np.linalg.norm(local,axis=2)).max()),
          maximum_displacement_m=float(np.linalg.norm(new-sections,axis=2).max()),
          source_length_m=float(np.linalg.norm(np.diff(cc,axis=0),axis=1).sum()),candidate_length_m=float(np.linalg.norm(np.diff(new.mean(1),axis=0),axis=1).sum())))
    previous=json.loads((BASE/'takeoff-and-root-controls-v10/evidence.json').read_text());meshes=[];whole=[]
    for pi,pid in enumerate(previous['whole_runtime_qualified_blood_source_ids']):
        mesh=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres');vv=mesh.vertices.copy();attrs={k:x.copy()for k,x in mesh.attributes.items()};changed=pid in position_maps
        if changed:
            for vi,p in enumerate(mesh.vertices):
                pkey=tuple(p)
                if pkey not in position_maps[pid]:continue
                vv[vi]=position_maps[pid][pkey]
                if 'normals'in attrs:attrs['normals'][vi]=np.rint(np.clip(rotation_maps[pid][pkey]@np.asarray(mesh.attributes['normals'][vi],float),-32767,32767)).astype(np.int16)
        candidate=Mesh(vv,mesh.faces,mesh.frame,attrs,dict(part_id=pid,source_geometry_sha256=mesh.geometry_sha,
          contract_sha256=digest_file(folder/'contract.json'),modified_S2_components_only=changed,finite_wall_acceptance=False))
        path=out/f'whole-family-part-{pi:02d}.npz';candidate.save(path);meshes.append(dict(id=pid,mesh=str(path.resolve()),geometry_sha256=candidate.geometry_sha))
        whole.append(dict(id=pid,changed=changed,source_geometry_sha256=mesh.geometry_sha,candidate_geometry_sha256=candidate.geometry_sha,
          source_faces_unchanged=np.array_equal(mesh.faces,candidate.faces),normalSigns_unchanged=all(np.array_equal(mesh.attributes[k],attrs[k])for k in ['normalSigns']if k in attrs),
          source_UV_colors_unchanged=all(np.array_equal(mesh.attributes[k],attrs[k])for k in ['uvs','colors']if k in attrs)))
        if changed:
            for row in selected:
                if row['part_id']!=pid:continue
                ff=a[row['source_faces_array']];used,inv=np.unique(mesh.faces[ff],return_inverse=True)
                fragment=Mesh(vv[used],inv.reshape(-1,3),mesh.frame,{k:x[used]for k,x in attrs.items()},dict(source_part_id=pid,source_component=row['component'],source_face_ids=ff.tolist()))
                fragment_path=out/(row['key']+'-candidate-fragment.npz');fragment.save(fragment_path)
                fragment_paths.append(dict(part_id=pid,component=row['component'],mesh=str(fragment_path.resolve()),geometry_sha256=fragment.geometry_sha))
    root_replay=[]
    for row in c['native_arterial_origins_and_parent_ownership']:
        pp=a[row['profile_array']];new=np.array([position_maps[row['part_id']][tuple(p)]for p in pp]);root_replay.append(dict(part_id=row['part_id'],source_component=row['component'],maximum_origin_profile_error_m=float(np.linalg.norm(new-pp,axis=1).max())))
    np.savez_compressed(out/'literal-proper-ring-transport.npz',**saved,canonical_source_centres=coords,canonical_candidate_centres=target,
                        source_bronchial_projection=q,current_HRA_course_base=base,unprotected_offset_candidate=mapped,
                        proper_rotations=rotation,root_protection_weights=weight)
    offset_before=np.linalg.norm(coords-q,axis=1);offset_after=np.linalg.norm(mapped-base,axis=1)
    result=dict(schema='right-S2-two-origin-terminal-bound-proper-ring-course-candidate/v22',source_catalogue_sha256=c['source_catalogue_sha256'],
      construction_contract_sha256=digest_file(folder/'contract.json'),source_semantic_root_review_sha256=c['source_semantic_root_review_sha256'],
      target_terminal_evidence_sha256=c['current_recipient_terminal_evidence_sha256'],
      candidate_meshes=meshes,selected_fragment_meshes=fragment_paths,source_and_candidate_part_replay=whole,proper_ring_metrics=metrics,
      native_origin_profile_replay=root_replay,
      maximum_original_bronchovascular_offset_length_error_m=float(np.abs(offset_before-offset_after).max()),
      every_literal_ring_calibre_preserved=all(r['maximum_relative_section_area_error']<1e-10 and r['maximum_literal_profile_pairwise_error_m']<1e-12 for r in metrics),
      all_other_52_family_source_geometries_literal=all(r['source_geometry_sha256']==r['candidate_geometry_sha256']for r in whole if not r['changed']),
      arterial_incidence_preserved_independently_of_bronchial_guide_topology=True,
      actual_bronchial_triangle_clearance_pending=True,actual_whole_context_contacts_pending=True,
      actual_selected_surface_self_intersections_pending=True,actual_finite_communicating_wall_pending=True,render_pending=True,
      all_heart_skeleton_current_lung_airway_and_nine_lingular_source_dependencies_fixed=True,
      arrays_sha256=digest_file(out/'literal-proper-ring-transport.npz'),
      geometry_candidate_only=True,catalogue_changed=False,installed=False,anatomical_acceptance=False,
      elapsed_seconds=time.monotonic()-began)
    (out/'candidate-evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),modified_parts=sum(r['changed']for r in whole),selected_source_components=len(metrics),
      all_other_52_sources_literal=result['all_other_52_family_source_geometries_literal'],
      ring_calibre_preserved=result['every_literal_ring_calibre_preserved'],
      maximum_origin_error_m=max(r['maximum_origin_profile_error_m']for r in root_replay),
      evidence_sha256=digest_file(out/'candidate-evidence.json'),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['contract','topology','ordered','candidate'],required=True);args=parser.parse_args()
    if args.stage=='contract':contract()
    elif args.stage=='topology':topology()
    elif args.stage=='ordered':ordered()
    else:candidate()
