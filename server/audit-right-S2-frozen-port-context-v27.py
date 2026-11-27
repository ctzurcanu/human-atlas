"""Actual native origin/prefix contact before any further arterial solve."""
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import winding_contains
from geometry_tools.algorithms.tube_sections import arclength

BASE = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OUT = BASE / 'right-S2-frozen-port-context-v27'
FRAME = 'model:male:assembled-metres'


def main():
    began = time.monotonic(); OUT.mkdir()
    f = json.loads((BASE / 'family-contract-v2.json').read_text()); assert digest_file(Path(f['source_catalogue'])) == f['source_catalogue_sha256']
    atlas = json.loads(Path(f['source_catalogue']).read_text()); parts = {p['id']:p for p in atlas['parts']}; reader = AtlasReader()
    contract = json.loads((BASE / 'right-S2-ordered-chart-contract-v21/contract.json').read_text())
    arrays = np.load(BASE / 'right-S2-ordered-chart-contract-v21/literal-ordered-guides-and-roots.npz', allow_pickle=False)
    ids = json.loads((BASE / 'right-S2-parent-first-course-v26/contract.json').read_text())['fixed_airway_ids']
    vv = []; ff = []; owner = []; owner_faces = []; offset = 0
    for oi, pid in enumerate(ids):
        m = reader.read_part(atlas, parts[pid], frame=FRAME)
        vv.extend(m.vertices); ff.extend(m.faces + offset); owner.extend([oi]*len(m.faces)); owner_faces.extend(range(len(m.faces))); offset += len(m.vertices)
    airway = Mesh(np.asarray(vv), np.asarray(ff), FRAME); tree = BVH(airway)
    territory = reader.read_part(atlas, parts['TA98PULMONARY:male:right-superior:posterior'], frame=FRAME)
    owner = np.asarray(owner); owner_faces = np.asarray(owner_faces); saved = {}; rows = []
    for root in contract['native_arterial_origins_and_parent_ownership']:
        row = next(r for r in contract['selected_source_curves'] if r['part_id'] == root['part_id'] and r['component'] == root['component'])
        key = row['key']; sections = arrays[row['source_sections_array']]; cc = sections.mean(1)
        oriented = sections[::-1]; oriented_cc = oriented.mean(1); distance = arclength(oriented_cc)
        last = int(np.searchsorted(distance, .004)); prefix = oriented[:last+1]
        source = reader.read_part(atlas, parts[root['part_id']], frame=FRAME)
        point_set = {tuple(p) for ring in prefix for p in ring}
        face_ids = arrays[row['source_faces_array']]
        keep = np.array([all(tuple(p) in point_set for p in tri) for tri in source.triangles[face_ids]])
        prefix_faces = face_ids[keep]; mesh = Mesh(source.vertices, source.faces[prefix_faces], FRAME)
        path = OUT / (key + '-literal-protected-prefix.npz'); mesh.save(path)
        contacts, summary = pair_contacts(mesh, airway, 0., second_bvh=tree, max_pairs=100000, max_tests=500000)
        for name, value in contacts.items(): saved[key + '_contacts_' + name] = value
        named = []
        wrong = np.isin(contacts['relation'], ['crossing', 'coplanar-overlap'])
        for oi in np.unique(owner[contacts['pairs'][wrong, 1]]):
            where = wrong & (owner[contacts['pairs'][:, 1]] == oi)
            pairs = contacts['pairs'][where]
            named.append(dict(airway_id=ids[int(oi)], crossing_face_pairs=int(where.sum()),
                              witnesses=[dict(native_source_face=int(prefix_faces[a]), airway_source_face=int(owner_faces[b]),
                                              source_point_m=contacts['first_points'][j].tolist() if 'first_points' in contacts else None)
                                         for j,(a,b) in list(zip(np.flatnonzero(where), pairs))[:8]]))
        profile = oriented[0]; near = [tree.closest_point(p) for p in profile]
        # Literal ring-edge crossing witnesses against the complete airway.
        # Strict plane crossing and triangle barycentrics; boundary/coplanar
        # cases remain explicitly unclassified, not silently clear.
        edge_hits = []; boundary_cases = 0; triangles = airway.triangles
        e1 = triangles[:,1]-triangles[:,0]; e2 = triangles[:,2]-triangles[:,0]
        normals = np.cross(e1,e2); nn = np.linalg.norm(normals,axis=1); normals /= nn[:,None]
        for ei,(a,b) in enumerate(zip(profile,np.roll(profile,-1,axis=0))):
            low=np.minimum(a,b); high=np.maximum(a,b)
            overlap=np.all(triangles.max(1)>=low,axis=1)&np.all(triangles.min(1)<=high,axis=1)
            face=np.flatnonzero(overlap); da=np.sum((a-triangles[face,0])*normals[face],axis=1); db=np.sum((b-triangles[face,0])*normals[face],axis=1)
            boundary_cases += int(np.sum((abs(da)<1e-12)|(abs(db)<1e-12)))
            face=face[da*db<0]
            for j in face:
                d0=float((a-triangles[j,0])@normals[j]);d1=float((b-triangles[j,0])@normals[j]);fraction=d0/(d0-d1);p=a+fraction*(b-a)
                basis=np.column_stack([e1[j],e2[j]]);uv=np.linalg.lstsq(basis,p-triangles[j,0],rcond=None)[0]
                if uv.min()>=0 and uv.sum()<=1:
                    edge_hits.append(dict(origin_profile_edge=ei,edge_fraction=float(fraction),airway_id=ids[int(owner[j])],
                                          airway_source_face=int(owner_faces[j]),point_m=p.tolist(),triangle_barycentric=[float(1-uv.sum()),float(uv[0]),float(uv[1])]))
        contained=winding_contains(territory,oriented_cc,distance_backend='bvh')
        inside=np.asarray(contained['inside']); transitions=np.flatnonzero(inside[:-1]!=inside[1:])
        saved[key+'_origin_profile']=profile;saved[key+'_prefix_sections']=prefix;saved[key+'_centre_distances_m']=distance;saved[key+'_S2_centre_inside']=inside
        rows.append(dict(key=key,TA98_construction_identity=root['TA98_construction_identity'],source_parent_part=root['parent_part'],
                         prefix_mesh=str(path.resolve()),prefix_source_face_ids=prefix_faces.tolist(),actual_prefix_length_m=float(distance[last]),
                         prefix_vs_all_27_actual_airways=summary,named_airway_crossings=named,origin_profile_strict_edge_crossing_witnesses=edge_hits,
                         origin_profile_boundary_or_coplanar_edge_cases=boundary_cases,
                         origin_profile_minimum_vertex_to_airway_distance_m=min(x['distance_m']for x in near),
                         origin_profile_vertex_nearest_owners=[ids[int(owner[x['face_id']])]for x in near],
                         source_centres_in_current_S2=inside.tolist(),source_centre_S2_classification_transition_intervals=transitions.tolist(),
                         no_claim_current_S2_source_classification_is_anatomical_entry=True,
                         S2_containment_must_follow_established_actual_segment_entry_not_fixed_arclength=True))
    np.savez_compressed(OUT/'literal-port-and-prefix-witnesses.npz',**saved)
    data=dict(schema='right-A2-native-frozen-port-versus-current-airways/v27',source_catalogue_sha256=f['source_catalogue_sha256'],
              actual_current_airway_ids=ids,actual_current_airway_geometry_sha256=airway.geometry_sha,
              actual_current_S2_geometry_sha256=territory.geometry_sha,rows=rows,
              artificial_20mm_S2_entry_assumption_released=True,source_profile_proximity_not_patency=True,
              geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,
              elapsed_seconds=time.monotonic()-began)
    (OUT/'evidence.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),sha256=digest_file(OUT/'evidence.json'),
                          rows=[dict(key=r['key'],prefix=r['prefix_vs_all_27_actual_airways'],named=r['named_airway_crossings'],
                                     origin_edge_hits=len(r['origin_profile_strict_edge_crossing_witnesses']),S2_transitions=r['source_centre_S2_classification_transition_intervals'])for r in rows],
                          elapsed_seconds=data['elapsed_seconds'])))


if __name__=='__main__':main()
