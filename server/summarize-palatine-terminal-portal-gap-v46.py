"""Explain the finite v45 graph frontier using saved nodes and saved certificates only."""
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.spatial import cKDTree
from collections import deque
from geometry_tools.common import atomic_json,digest_file

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';I=D/'palatine-surface-portal-graph-v44';O=D/'palatine-surface-portal-reachability-v45';E=I/'evidence.json';e=json.loads(E.read_bytes());T=O/'evidence.json';t=json.loads(T.read_bytes());assert t['completed']and not t['source_to_SG_SP_band_connection_found'];assert digest_file(O/'reachability-progress.npz')==t['progress_sha256'];p=np.load(O/'reachability-progress.npz');q=p['points'];terminal=p['reachable'];pairs={tuple(map(int,a)):bool(b)for a,b in zip(p['edge_pairs'],p['edge_pass'])};proof={tuple(r['local_nodes']):(j,r)for j,r in enumerate(e['certificates'])if r.get('local_nodes')is not None};seeds={int(r['local_node'])for r in e['certificates']if r.get('grid_node')is not None and r['passed']};reached=np.zeros(len(q),bool);reached[list(seeds)]=True;queue=deque(sorted(seeds));neighbors=cKDTree(q).query_ball_point(q,.0015);blocked=[];accepted=[]
    while queue:
        a=queue.popleft()
        for b in neighbors[a]:
            if reached[b]:continue
            pair=(min(a,b),max(a,b));ok=pairs[pair];ix,row=proof[pair]
            if ok:reached[b]=True;queue.append(b);accepted.append(dict(nodes=list(pair),prior_proof_index=ix))
            else:blocked.append(dict(nodes=list(pair),points_m=q[list(pair)].tolist(),prior_proof_index=ix,prior_proof_sha256=digest_file(E),certificate=row['proof']))
    assert np.array_equal(reached,terminal)and len(blocked)==3 and len(accepted)==2
    final_cross_pairs=[]
    for a in np.flatnonzero(~terminal):
        for b in neighbors[a]:
            if terminal[b]:final_cross_pairs.append((int(a),int(b)))
    targets=np.array(t['target_nodes']);reach_ids=np.flatnonzero(terminal);dist,ix=cKDTree(q[reach_ids]).query(q[targets]);gaps=[]
    for target,dd,k in zip(targets,dist,ix):
        other=int(reach_ids[k]);gaps.append(dict(target_node=int(target),target_point_m=q[target].tolist(),nearest_reachable_node=other,nearest_reachable_point_m=q[other].tolist(),point_cloud_gap_m=float(dd),allowed_local_neighbor_distance_m=.0015,no_segment_or_tissue_clearance_inferred=True))
    bank=np.load(D/'global-palatine-actual-radius-v39/literal-surface-query-bank-v40.npz');binding=json.loads((D/'global-palatine-actual-radius-v39/edge-search-evidence-v40.json').read_bytes())
    for row in blocked:
        row['both_nodes_became_reachable_by_other_certified_edges']=bool(terminal[row['nodes']].all())
        cert=row['certificate']
        if cert['surface']=='composite':
            face=int(cert['face']);source_row=binding['source_bindings'][int(bank['source_index'][face])];row['literal_source_witness']=dict(part_id=source_row['part_id'],geometry_sha256=source_row['geometry_sha256'],source_face=int(bank['native_face'][face]),combined_query_face=face)
    result=dict(schema='palatine-terminal-finite-portal-graph-next-action/v46',source_catalogue_sha256=t['source_catalogue_sha256'],terminal_evidence_sha256=digest_file(T),fixed_node_bank_sha256=digest_file(I/'local-graph-progress.npz'),replay_matches_terminal_progress=True,reached_nodes=int(terminal.sum()),unreached_nodes=int((~terminal).sum()),three_blocked_edges_encountered_during_traversal=blocked,two_accepted_bridges_reused=accepted,
      final_reachable_to_unreachable_pairs_within_1p5mm=len(final_cross_pairs),four_target_geometric_gaps=gaps,
      precise_interpretation='The three rejected edges encountered during reachability no longer form the terminal cut: both endpoints of each are reachable through other certified edges. At the terminal state there are zero sampled reachable/unreachable neighbor pairs within 1.5 mm. The four required SG/SP band nodes are geometrically separated from sampled reachable nodes. This is finite surface-sample graph nonconnection, not a verified anatomical wall obstruction.',
      next_principled_action='Inspect the actual opposing SG/SP face neighborhoods between these literal nearest pairs; construct a continuous surface-distance corridor certificate with constrained medial patches/adaptive subdivision of actual source faces. Determine whether the intervening source geometry permits the full fixed radius before adding any connection. No random densification, global field rerun, takeoff shift or radius reduction.',
      anatomy_limits='Open/nonmanifold source occupied volumes and complete oral-fluid exclusion remain unresolved. Even a future graph connection will require a C2 course, emitted wall patency, complete tissue/lumen packing and rendered review.',all_geometry_and_install_held=True,no_new_queries_or_geometry=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'terminal-next-action-v46.json',result);print(json.dumps(dict(reached_nodes=result['reached_nodes'],unreached_nodes=result['unreached_nodes'],final_cross_pairs=len(final_cross_pairs),target_gaps_mm=[r['point_cloud_gap_m']*1000 for r in gaps],three_blocked_edges=[dict(nodes=r['nodes'],certificate=r['certificate'],both_nodes_now_reachable=r['both_nodes_became_reachable_by_other_certified_edges'])for r in blocked],elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
