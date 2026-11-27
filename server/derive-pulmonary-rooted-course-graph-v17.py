"""Literal source course-node graph and fixed recipient course-chart prefixes.

Coincident centre nodes are construction correspondences, never evidence of
communicating vessel tissue. Current HRA airway interface charts remain actual
opposed boundary edges. No global displacement field is fitted here.
"""
import json
import time
import importlib.util
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,StationTopologyError
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
_spec=importlib.util.spec_from_file_location('pulmonary_course_charts_v12',Path(__file__).with_name('derive-pulmonary-diverged-course-controls-v12.py'))
_helper=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_helper)
intrinsic_profile=_helper.intrinsic_profile

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def literal_graph(curves,arrays,prefix):
    points=[];members=[];offsets=[]
    for i,row in enumerate(curves):
        cc=arrays[row['centres_array']];offsets.append(len(points));points.extend(cc);members.extend((i,j)for j in range(len(cc)))
    points=np.asarray(points);parent=np.arange(len(points));rank=np.zeros(len(points),int)
    def find(a):
        while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
        return int(a)
    for a,b in sorted(cKDTree(points).query_pairs(1e-6)):
        # Arterial/venous crossings are not graph branch nodes.
        if curves[members[a][0]]['system']!=curves[members[b][0]]['system']:continue
        aa,bb=find(a),find(b)
        if aa!=bb:parent[bb]=aa
    labels=np.array([find(i)for i in range(len(points))]);groups=defaultdict(list)
    for i,label in enumerate(labels):groups[int(label)].append(i)
    junctions=[];nodes={};nextnode=0
    for label,ids in groups.items():
        branches=sorted(set(members[i][0]for i in ids));endpoint=any(members[i][1]in [0,len(arrays[curves[members[i][0]]['centres_array']])-1]for i in ids)
        if len(branches)<2 and not endpoint:continue
        key=prefix+f'_node{nextnode}';nextnode+=1;nodes[label]=key;arrays[key+'_literal_centre_representations']=points[ids]
        junctions.append(dict(key=key,curve_station_members=[dict(curve=members[i][0],station=members[i][1],
          part_id=curves[members[i][0]]['part_id'],source_component=curves[members[i][0]]['component'])for i in ids],
          literal_representations_array=key+'_literal_centre_representations',
          maximum_pairwise_centre_difference_m=float(np.linalg.norm(points[ids,None]-points[None,ids],axis=2).max()),
          near_coincident_section_centres_not_shared_surface=True,communicating_tissue_established=False))
    links=[]
    for i,row in enumerate(curves):
        count=len(arrays[row['centres_array']]);stations=[j for j in range(count)if labels[offsets[i]+j]in nodes]
        for a,b in zip(stations,stations[1:]):
            left,right=nodes[labels[offsets[i]+a]],nodes[labels[offsets[i]+b]]
            if left==right:continue
            links.append(dict(curve=i,part_id=row['part_id'],source_component=row['component'],
              source_station_interval=[a,b],source_node=left,target_node=right,
              orientation_is_literal_serial_order_not_physiological_flow=True))
    return dict(nodes=junctions,links=links,source_curves=curves,
      geometric_branch_node_tolerance_m=1e-6,
      actual_wall_junctions_not_established_by_centre_node_incidence=True)


def main():
    began=time.monotonic();family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue'])
    assert digest_file(path)==family['source_catalogue_sha256'];atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    prior=json.loads((BASE/'takeoff-and-root-controls-v10/evidence.json').read_text())
    interfaces=json.loads((BASE/'source-correspondence-contract-v7/contract.json').read_text())['airway_interface_landmarks']
    charts=np.load(BASE/'source-correspondence-contract-v7/literal-fixed-airway-interface-controls.npz',allow_pickle=False)
    out=BASE/'rooted-source-recipient-graph-v17';out.mkdir();arrays={};vessels=[];guides=[];unsupported=[];annular=[]
    source_guides=[r['source_id']for r in prior['named_control_pairs']]+[
      'ZA:Right main bronchus','ZA:Right superior lobar bronchus','ZA:Middle lobar bronchus.r','ZA:Right inferior lobar bronchus',
      'ZA:Left main bronchus','ZA:Left superior lobar bronchus','ZA:Left inferior lobar bronchus']
    for group,ids,records in [('vascular',prior['whole_runtime_qualified_blood_source_ids'],vessels),('bronchial',source_guides,guides)]:
        for pi,pid in enumerate(ids):
            mesh=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres');v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces]
            try:tubes=recover_station_rings(v,f)
            except StationTopologyError as error:
                unsupported.append(dict(id=pid,group=group,reason=str(error),geometry_sha256=mesh.geometry_sha));continue
            for ci,tube in enumerate(tubes):
                key=group+f'_p{pi}_c{ci}';sections=v[np.array(tube['rings'])];cc=sections.mean(1)
                arrays[key+'_sections']=sections;arrays[key+'_centres']=cc;arrays[key+'_source_face_ids']=tube['face_ids']
                records.append(dict(part_id=pid,component=ci,system=parts[pid]['system']if group=='vascular'else'airway-reference',
                  key=key,centres_array=key+'_centres',sections_array=key+'_sections',source_faces_array=key+'_source_face_ids',
                  stations=len(cc),ring_points=sections.shape[1],geometry_sha256=mesh.geometry_sha,
                  source_hidden_reference=group=='bronchial',flow_orientation_unproven=True))
    vascular_graph=literal_graph(vessels,arrays,'vascular');bronchial_graph=literal_graph(guides,arrays,'bronchial')
    # Retain actual finite lingular interface surfaces, including their true
    # boundary loops, without replacing annuli with centre-point holes.
    for pid in prior['whole_runtime_qualified_blood_source_ids']:
        if not parts[pid].get('pairedTissue'):continue
        mesh=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres')
        for gi,group in enumerate(parts[pid].get('tissueDrawGroups',[])):
            if not group.get('neighbor'):continue
            ff=mesh.faces[group['start']//3:(group['start']+group['count'])//3]
            vv,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);ff=inv[ff];key='finite_interface_'+str(len(annular))
            arrays[key+'_faces']=ff;arrays[key+'_vertices']=vv
            try:
                loops=boundary_loops(ff)
                for li,loop in enumerate(loops):arrays[key+f'_boundary_loop_{li}']=vv[loop]
                loop_counts=[len(x)for x in loops]
            except ValueError as error:loop_counts=[]
            annular.append(dict(part_id=pid,neighbor=group['neighbor'],source_index_interval=[group['start'],group['start']+group['count']],
              faces_array=key+'_faces',vertices_array=key+'_vertices',actual_boundary_loop_points=loop_counts,
              actual_tissue_surface_not_new_opening=True))
    # Each segment's entry is its actual named parent interface. Combined
    # donor BI+II is represented as two separate recipient chart sequences.
    targets={}
    for r in prior['named_control_pairs']:
        for chart in r['target_actual_boundary_charts']:
            targets[chart['part_id']]=(r['target_parent_id'],chart)
    for child,parent in [
      ('right_superior_lobar_bronchus','right_main_bronchus'),
      ('right_inferior_lobar_bronchus','right_main_bronchus'),
      ('right_middle_lobar_bronchus','right_inferior_lobar_bronchus'),
      ('left_superior_lobar_bronchus','left_main_bronchus'),
      ('left_inferior_lobar_bronchus','left_main_bronchus')]:
        pid='LOCAL:male:pulmonary:bronchus:'+child;parent='LOCAL:male:pulmonary:bronchus:'+parent
        found=[r for r in interfaces if r['part_id']==pid and parent in r['other_branch_ids']and r['all_edges_have_literal_opposed_other_branch']]
        if len(found)==1:targets[pid]=(parent,found[0])
    recipient=[]
    for ti,(pid,(parent,chart))in enumerate(targets.items()):
        mesh=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres');v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces]
        seed=charts[chart['array_key']];key=f'recipient{ti}';arrays[key+'_literal_parent_entry']=seed
        samples=[];stopped=None
        for level in [.001,.003,.006,.010,.015,.020,.030,.040,.050,.060,.080]:
            pp,proof=intrinsic_profile(v,f,seed,level)
            if pp is None:stopped=dict(edge_distance_m=level,reason=proof);break
            pk=key+f'_course_{len(samples)}';arrays[pk]=pp
            samples.append(dict(literal_chart_array=pk,chart_centre_m=pp.mean(0).tolist(),edge_graph_distance_m=level,proof=proof))
        recipient.append(dict(part_id=pid,parent_id=parent,geometry_sha256=mesh.geometry_sha,
          actual_parent_entry_chart=chart,parent_entry_array=key+'_literal_parent_entry',
          unique_common_prefix_course_charts=samples,first_unresolved_contour=stopped,
          chart_centres_are_initialization_landmarks_not_airway_wall_or_vascular_lumen=True,
          distal_fork_not_replaced_with_an_average_of_child_entries=True))
    np.savez_compressed(out/'literal-source-nodes-and-recipient-course-charts.npz',**arrays)
    result=dict(schema='pulmonary-literal-source-rooted-node-and-recipient-course-graph/v17',
      source_catalogue_sha256=family['source_catalogue_sha256'],current_runtime_evidence_sha256=digest_file(BASE/'actual-runtime-context-v9/runtime-context.json'),
      arrays_sha256=digest_file(out/'literal-source-nodes-and-recipient-course-charts.npz'),
      vascular_source_graph=vascular_graph,original_hidden_bronchial_reference_graph=bronchial_graph,
      current_HRA_recipient_interface_graph=interfaces,current_HRA_recipient_course_prefixes=recipient,
      existing_finite_lingular_interfaces=annular,unsupported_serial_sources=unsupported,
      construction_contract=[
        'Map a donor common anatomical parent before child divergence; never independently interpolate coincident donor nodes to different child ports.',
        'Use fixed recipient common-parent entry and actual course chart sequences; retain original bronchial references separately as semantic witnesses.',
        'Construct one shared vascular junction from declared actual branch mouths before assigning TA98 tissue owners.',
        'Transport each literal profile by a proper rotation along the reconstructed smooth branch course; retain radii and explicit curvature/calibre tolerances.',
        'A centre-node graph is a construction guide; it does not establish an existing communicating wall or lumen.',
        'Root patches and all fixed proximal profiles must be declared before any new finite-wall reconstruction.'],
      global_compressed_field_rejected=True,new_field_trial=False,geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,
      elapsed_seconds=time.monotonic()-began)
    (out/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),vascular_curves=len(vessels),bronchial_reference_curves=len(guides),
      vascular_construction_nodes=len(vascular_graph['nodes']),vascular_curve_intervals=len(vascular_graph['links']),
      recipient_named_prefixes=len(recipient),recipient_actual_prefix_charts=sum(len(r['unique_common_prefix_course_charts'])for r in recipient),
      finite_lingular_interfaces=len(annular),evidence_sha256=digest_file(out/'evidence.json'),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
