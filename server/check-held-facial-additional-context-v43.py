"""Actual held candidate surface crossings with the 80 newly bound physical sources."""
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'global-palatine-actual-radius-v39';R=D/'global-palatine-graph-v33/additional-participation-frame-v37.json';r=json.loads(R.read_bytes());ids=r['retained_indexed_surface_or_tissue_ids'];assert len(ids)==80;source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(source)==r['source_catalogue_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};status={p['id']:p for p in r['rows']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';rows=[];crossings=[];complete=True;branches={}
    for name in ['ascending-palatine','tonsillar','glandular-0','glandular-1']:
        path=D/'feasibility-v2'/(name+'-actual-wall.npz');m=Mesh.load(path);branches[name]=(m,BVH(m),path)
    for pid in ids:
        if time.monotonic()-began>48:complete=False;break
        s=status[pid];m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=r['source_catalogue_sha256']);assert m.geometry_sha==s['registry_geometry_binding']['geometry_sha256'];used=m.vertices[np.unique(m.faces)];lo,hi=used.min(0),used.max(0);tree=None
        for name,(branch,bvh,path)in branches.items():
            if np.any(hi<branch.bounds[0])or np.any(lo>branch.bounds[1]):rows.append(dict(branch=name,part_id=pid,geometry_sha256=m.geometry_sha,complete=True,excluded_by_literal_triangle_bounds=True));continue
            if tree is None:tree=BVH(m)
            result,summary=pair_contacts(branch,m,0.,first_bvh=bvh,second_bvh=tree,max_pairs=250000,max_tests=600000)
            row=dict(branch=name,part_id=pid,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,source_frame=frame,actual_emitted_branch_sha256=digest_file(path),summary=summary,complete=summary['complete_within_pair_budget'],closed_occupied_volume_not_established_here=True)
            if summary['crossing_face_pairs']or summary['coplanar_overlap_face_pairs']:
                selected=np.flatnonzero(np.isin(result['relation'],['crossing','coplanar-overlap']));p=O/('held-'+name+'-extra-context-'+str(len(crossings))+'.npz');atomic_npz(p,**{k:v[selected]for k,v in result.items()});row['actual_crossing_witness_file']=str(p.relative_to(ROOT));row['actual_crossing_witness_sha256']=digest_file(p);crossings.append(row)
            rows.append(row)
    output=dict(schema='held-facial-additional-default-context-surface-audit/v43',source_catalogue_sha256=r['source_catalogue_sha256'],runtime_frame_cache_sha256=digest_file(R),candidate_fixed_v17=True,checked_source_count=len({x['part_id']for x in rows}),expected_source_count=80,branch_source_pairs=len(rows),complete=complete and len({x['part_id']for x in rows})==80 and all(x['complete']for x in rows),actual_crossings=crossings,rows=rows,not_new_global_route_geometry=True,no_occupied_volume_or_anatomical_acceptance=True,no_source_mutation_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'held-v17-additional-context-v43.json',output);print(json.dumps(dict(complete=output['complete'],checked_source_count=output['checked_source_count'],branch_source_pairs=len(rows),actual_crossings=[dict(branch=x['branch'],part_id=x['part_id'],crossing_face_pairs=x['summary']['crossing_face_pairs'],coplanar_overlap=x['summary']['coplanar_overlap_face_pairs'])for x in crossings],elapsed_seconds=output['elapsed_seconds'])))

if __name__=='__main__':main()
