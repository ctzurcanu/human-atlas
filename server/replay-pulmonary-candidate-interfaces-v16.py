"""Replay literal paired interfaces and protected rims of diagnostic v13."""
import json
import time
from collections import defaultdict
from pathlib import Path
import numpy as np
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader,Mesh

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began=time.monotonic();root=BASE/'rooted-course-initialization-v13';contract=json.loads((root/'contract.json').read_text())
    assert digest_file(Path(contract['source_catalogue']))==contract['source_catalogue_sha256']
    atlas=json.loads(Path(contract['source_catalogue']).read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    candidate=json.loads((root/'candidate-evidence.json').read_text());source={};after={};mapping={};edge_map=defaultdict(list)
    for row in candidate['candidate_meshes']:
        pid=row['id'];s=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres');d=Mesh.load(row['mesh'])
        assert np.array_equal(s.faces,d.faces)
        assert s.geometry_sha==d.source['source_geometry_sha256']
        source[pid]=s;after[pid]=d
        mapping[pid]={tuple(x):y for x,y in zip(s.vertices,d.vertices)}
        for fi,face in enumerate(s.faces):
            for a,b in zip(face,np.roll(face,-1)):
                aa,bb=tuple(s.vertices[a]),tuple(s.vertices[b]);edge_map[tuple(sorted((aa,bb)))].append((pid,fi,int(a),int(b)))
    by_pair={};worst=0.;count=0
    for key,rows in edge_map.items():
        owners=sorted(set(r[0]for r in rows))
        if len(owners)<2:continue
        for i,left in enumerate(owners):
            for right in owners[i+1:]:
                gap=max(np.linalg.norm(mapping[left][p]-mapping[right][p])for p in key)
                group=by_pair.setdefault((left,right),dict(shared_source_edges=0,maximum_endpoint_gap_m=0.))
                group['shared_source_edges']+=1;group['maximum_endpoint_gap_m']=max(group['maximum_endpoint_gap_m'],float(gap))
                worst=max(worst,float(gap));count+=1
    protection=[];arrays=np.load(root/'literal-protected-and-shared-profiles.npz',allow_pickle=False)
    for row in contract['root_profile_protection']:
        pid=row['part_id'];pp=arrays[row['profile_array']];dd=np.array([mapping[pid][tuple(p)]for p in pp])
        protection.append(dict(part_id=pid,profile=row['profile_array'],maximum_coordinate_error_m=float(np.linalg.norm(dd-pp,axis=1).max())))
    rigid=np.load(root/'literal-ring-transport-and-jacobians.npz',allow_pickle=False)
    anchor=rigid['lingular_anchor_source'];target=rigid['lingular_anchor_candidate'];rotation=rigid['lingular_proper_rotation']
    rigid_rows=[]
    for pid in contract['rigid_lingular_closure']:
        exact=(source[pid].vertices-anchor)@rotation.T+target
        rigid_rows.append(dict(part_id=pid,maximum_rigid_replay_error_m=float(np.linalg.norm(exact-after[pid].vertices,axis=1).max())))
    result=dict(schema='pulmonary-diagnostic-candidate-literal-interface-replay/v16',
      source_catalogue_sha256=contract['source_catalogue_sha256'],candidate_evidence_sha256=digest_file(root/'candidate-evidence.json'),
      literal_cross_owner_source_edge_comparisons=count,
      maximum_cross_owner_candidate_endpoint_gap_m=worst,
      cross_owner_pairs=[dict(left=k[0],right=k[1],**v)for k,v in sorted(by_pair.items())],
      protected_rim_coordinates=protection,lingular_rigid_replays=rigid_rows,
      closed_and_open_source_communication_not_inferred_from_this_replay=True,
      rejected_for_installation=True,
      rejection_reason='Coarse named-course residual remains 13.944 mm RMS/34.075 mm max; sampled initialization singular value 0.061 and nonserial area compression remain unacceptable.',
      no_coefficient_trial_authorized=True,catalogue_changed=False,anatomical_acceptance=False,
      elapsed_seconds=time.monotonic()-began)
    out=BASE/'literal-interface-replay-v16';out.mkdir();(out/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),literal_shared_edges=count,maximum_shared_endpoint_gap_m=worst,
       maximum_rim_error_m=max(r['maximum_coordinate_error_m']for r in protection),
       evidence_sha256=digest_file(out/'evidence.json'),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
