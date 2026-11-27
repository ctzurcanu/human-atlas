"""Join literal source material siblings; never repair or cap native anatomy."""
from pathlib import Path
import sys,json,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import manifold3d as md
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=base/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'native-classifiers-v20';assert not O.exists();O.mkdir()
    evidence=json.loads((D/'feasibility-v2/evidence.json').read_bytes());plan=json.loads((base/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());binding=json.loads((D/'source-binding-v19/evidence.json').read_bytes());source=ROOT/binding['snapshot'];expected=binding['snapshot_sha256'];assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();unknown=evidence['unverified_native_occupied_volumes'];groups={}
    for row in unknown:
        part=parts[row['part_id']];key=part.get('sourceId')or row['part_id'];groups.setdefault(key,[]).append(row)
    spec=importlib.util.spec_from_file_location('source_classification_v18',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    wall=Mesh.load(D/'finite-wall.npz');origin=wall.bounds.mean(0);parent,parentmode=helper.occupancy_classifier(Mesh.load(ROOT/plan['fixed_parent']['mesh']),True);assert parent is not None
    parent_solid=helper.solid(parent,origin);occupied={}
    for route in plan['routes']:
        z=np.load(D/(route['name']+'-fields.npz'));m=Mesh(z['fields'][0],z['faces'],wall.frame);c,_=helper.occupancy_classifier(m,True);assert c is not None;occupied[route['name']]=helper.solid(c,origin)-parent_solid
    checkpoint=D/'native-classifiers-v19/progress.json';saved=json.loads(checkpoint.read_bytes());assert saved['source_catalogue_sha256']==expected
    rows=saved['rows'];finished={r['source_object_id']for r in rows}
    for sourceid,affected in groups.items():
        if sourceid in finished:continue
        if time.monotonic()-began>48:raise RuntimeError('48-second literal native classifier work budget reached')
        original=parts[affected[0]['part_id']];sourcekey=original.get('sourceId')
        siblings=[p for p in atlas['parts']if p['id'].startswith('ZA:')and not p.get('suppressed')and p.get('sourceId')==sourcekey]if sourcekey else[original]
        assert all(p.get('sourceId')==sourcekey for p in siblings)
        meshes=[reader.read_part(atlas,p,frame=wall.frame,catalogue_path=source,catalogue_sha=expected)for p in siblings];vertices=[];faces=[];offset=0
        for mesh in meshes:vertices.append(mesh.vertices);faces.append(mesh.faces+offset);offset+=len(mesh.vertices)
        raw=np.concatenate(vertices);v,inverse=np.unique(raw,axis=0,return_inverse=True);f=inverse[np.concatenate(faces)];joined=Mesh(v,f,wall.frame);gate=topology(joined)
        row=dict(source_object_id=sourceid,source_membership_basis='Literal installed ZA sourceId shared by material partitions; no semantic-name aggregation.',part_ids=[p['id']for p in siblings],parts=[dict(part_id=p['id'],sourceId=p.get('sourceId'),conceptId=p.get('conceptId'),material=p.get('material'),geometry_sha256=m.geometry_sha)for p,m in zip(siblings,meshes)],source_coordinates_literal=True,source_triangle_corner_coordinates_literal=bool(np.array_equal(v[f],raw[np.concatenate(faces)])),exact_equal_coordinate_vertex_join_only=True,joined_geometry_sha256=joined.geometry_sha,joined_topology=gate,affected_branch_pairs=affected,cap_faces_added=0)
        if gate['closed_oriented_manifold']:
            file=O/('closed-source-'+str(len(rows))+'.npz');joined.save(file);row.update(classification='closed-source-whole',classifier=str(file.relative_to(ROOT)),classifier_sha256=digest_file(file),overlaps=[])
            solid=helper.solid(joined,origin)
            for name in sorted(set(r['branch']for r in affected)):
                overlap=occupied[name]^solid;assert overlap.status()==md.Error.NoError;row['overlaps'].append(dict(branch=name,occupied_overlap_mm3=float(overlap.volume()),empty=bool(overlap.is_empty())))
        elif original.get('system')=='integumentary':
            row.update(classification='actual-open-integumentary-surface',applicability='Actual source triangle exclusion applies; this open skin/hair/reference-region surface establishes no occupied 3D solid. Native full skin containment is separately unverified.')
        elif sourceid=='Pharynx.j':
            row.update(classification='two-triangle-nonvolume-native-reference',applicability='This named source has two triangles and no pharyngeal wall/airspace solid; cannot certify an oral/pharyngeal lumen.')
        else:
            row.update(classification='unresolved-physical-source-volume',applicability='All literal sibling material triangles joined exactly; complete occupied anatomical object still not established. Do not infer empty space or a valid solid from generalized winding.')
        rows.append(row);atomic_json(O/'progress.json',dict(source_catalogue_sha256=expected,completed_source_objects=len(rows),total_source_objects=len(groups),rows=rows,elapsed_seconds=time.monotonic()-began),overwrite=True)
    unresolved=[r['source_object_id']for r in rows if r['classification']=='unresolved-physical-source-volume'];overlap=[dict(source=r['source_object_id'],**x)for r in rows for x in r.get('overlaps',[])if not x['empty']]
    result=dict(schema='literal-source-object-classifiers/v18',source_catalogue_sha256=expected,candidate_geometry_sha256=wall.geometry_sha,unique_native_source_objects=len(rows),resolved_closed_sources=[r['source_object_id']for r in rows if r['classification']=='closed-source-whole'],open_nonvolume_sources=[r['source_object_id']for r in rows if r['classification']in{'actual-open-integumentary-surface','two-triangle-nonvolume-native-reference'}],unresolved_physical_source_objects=unresolved,forbidden_closed_occupied_overlap=overlap,oral_pharyngeal_lumen_full_proof_pending=True,source_aware_occupied_gate_pass=not unresolved and not overlap and False,installed=False,blender_queued=False,rows=rows,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);assert digest_file(source)==expected;print(json.dumps({k:v for k,v in result.items()if k!='rows'}))


if __name__=='__main__':main()
