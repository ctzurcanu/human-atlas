"""Conservative whole-source convex separation; no native closure or repair.

A strictly separating plane is checked against every source corner and every
actual emitted conduit corner. Failure to find one is only an unresolved
conservative envelope, never an anatomical collision claim.
"""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import linprog
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=base/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'conservative-envelopes-v21';assert not O.exists();O.mkdir()
    prior=json.loads((D/'native-classifiers-v20/evidence.json').read_bytes());binding=json.loads((D/'source-binding-v19/evidence.json').read_bytes());source=ROOT/binding['snapshot'];expected=binding['snapshot_sha256'];assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();wall=Mesh.load(D/'finite-wall.npz');bank=np.load(D/'surfaces-and-ancestry.npz');ports=json.loads((D/'ports.json').read_bytes());names=ports['component_names'];origin=wall.bounds.mean(0)
    candidate={};candidate_hashes={}
    for ci,name in enumerate(names[1:],1):
        faces=np.concatenate([bank['outer_faces'][bank['outer_component']==ci],bank['inner_faces'][bank['inner_component']==ci],bank[f'annulus_{ci+1}_faces']]);emitted=bank['vertices'][np.unique(faces)]
        fieldpath=D/(name+'-fields.npz');field=np.load(fieldpath)
        # Include the complete original branch profile carrier as a superset of
        # the Boolean-clipped actual emitted corner set, including its root.
        points=np.unique(np.concatenate([emitted,field['fields'].reshape(-1,3)]),axis=0);hull=ConvexHull((points-origin)*1000)
        candidate[name]=(points,(points[hull.vertices]-origin)*1000)
        candidate_hashes[name]=dict(actual_emitted_face_count=len(faces),actual_emitted_corner_points=len(emitted),conservative_total_points=len(points),fields_sha256=digest_file(fieldpath),whole_wall_geometry_sha256=wall.geometry_sha,conservative_superset_includes_complete_unclipped_branch_carrier=True)
    rows=[]
    for item in prior['rows']:
        if item['classification']!='unresolved-physical-source-volume':continue
        if time.monotonic()-began>45:raise RuntimeError('45-second conservative-separation work budget reached')
        meshes=[];source_rows=[]
        for part in item['parts']:
            mesh=reader.read_part(atlas,parts[part['part_id']],frame=wall.frame,catalogue_path=source,catalogue_sha=expected);assert mesh.geometry_sha==part['geometry_sha256'];meshes.append(mesh.vertices[np.unique(mesh.faces)]);source_rows.append(dict(part_id=part['part_id'],geometry_sha256=mesh.geometry_sha))
        allpoints=np.unique(np.concatenate(meshes),axis=0);hull=ConvexHull((allpoints-origin)*1000);hp=(allpoints[hull.vertices]-origin)*1000;results=[]
        for name in sorted(set(x['branch']for x in item['affected_branch_pairs'])):
            points,cp=candidate[name];A=np.r_[np.c_[cp,-np.ones(len(cp)),np.ones(len(cp))],np.c_[-hp,np.ones(len(hp)),np.ones(len(hp))]]
            solved=linprog([0,0,0,0,-1],A_ub=A,b_ub=np.zeros(len(A)),bounds=[(-1,1)]*3+[(None,None),(0,None)],method='highs',options={'time_limit':3.})
            row=dict(branch=name,solver_status=int(solved.status),source_hull_points=len(hp),conduit_hull_points=len(cp),certificate_pass=False)
            if solved.success and np.linalg.norm(solved.x[:3])>0:
                n=solved.x[:3]/np.linalg.norm(solved.x[:3]);source_min=float(((allpoints-origin)@n).min());candidate_max=float(((points-origin)@n).max());gap=source_min-candidate_max
                row.update(plane_normal=n.tolist(),plane_origin_m=origin.tolist(),plane_offset_from_origin_m=(source_min+candidate_max)/2,minimum_all_source_projection_m=source_min,maximum_all_conduit_projection_m=candidate_max,verified_full_point_separation_m=gap,verification_point_counts=[len(allpoints),len(points)],numerical_guard_m=2e-12,certificate_pass=bool(gap>2e-12))
            row['classification']='strict-conservative-envelope-separation'if row['certificate_pass']else'intersecting-or-unseparated-conservative-envelopes-not-an-actual-crossing'
            results.append(row)
        rows.append(dict(source_object_id=item['source_object_id'],full_literal_part_ids=item['part_ids'],source_geometry=source_rows,full_source_surface_points=len(allpoints),source_coordinates_literal=True,source_caps_or_repairs=0,results=results,all_affected_conduits_separated=all(r['certificate_pass']for r in results)))
        atomic_json(O/'progress.json',dict(source_catalogue_sha256=expected,completed_objects=len(rows),rows=rows,elapsed_seconds=time.monotonic()-began),overwrite=True)
    unresolved=[r['source_object_id']for r in rows if not r['all_affected_conduits_separated']];result=dict(schema='conservative-native-facial-envelope-certificates/v21',source_catalogue_sha256=expected,source_snapshot=str(source.relative_to(ROOT)),source_snapshot_sha256=digest_file(source),candidate_geometry_sha256=wall.geometry_sha,candidate_geometry_unchanged=True,conduit_inputs=candidate_hashes,source_object_count=len(rows),certified_source_objects=[r['source_object_id']for r in rows if r['all_affected_conduits_separated']],unique_unresolved_physical_objects=unresolved,passed_pair_certificates=sum(r['certificate_pass']for x in rows for r in x['results']),total_pair_certificates=sum(len(x['results'])for x in rows),intersecting_envelope_is_not_actual_crossing=True,oral_pharyngeal_lumen_exclusion_remains_separate=True,source_aware_gate_pass=False,rows=rows,installed=False,blender_queued=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);assert Mesh.load(D/'finite-wall.npz').geometry_sha==wall.geometry_sha and digest_file(source)==expected;print(json.dumps({k:v for k,v in result.items()if k not in{'rows','conduit_inputs'}}))


if __name__=='__main__':main()
