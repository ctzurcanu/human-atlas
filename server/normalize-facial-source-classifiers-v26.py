"""Classifier-only exact surface normalization; never fill or move anatomy."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import manifold3d as md
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import topology,audit


def oriented(tri):
    t=tuple(map(int,tri));return min(t,t[1:]+t[:1],t[2:]+t[:2])


def exact_components(v,f):
    tri=v[f];zero=np.all(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])==0,axis=1);seen={};kept=[];duplicates=[]
    for i,face in enumerate(f):
        if zero[i]:continue
        key=oriented(face)
        if key in seen:duplicates.append([i,seen[key]])
        else:seen[key]=i;kept.append(i)
    g=f[kept];assert set(oriented(t)for t in f[~zero])==set(oriented(t)for t in g)
    parent=np.arange(len(g));edges={}
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return int(i)
    for i,face in enumerate(g):
        for a,b in zip(face,np.roll(face,-1)):
            edge=tuple(sorted((int(a),int(b))))
            if edge in edges:parent[find(i)]=find(edges[edge])
            else:edges[edge]=i
    labels=np.array([find(i)for i in range(len(g))]);out=[]
    for label in np.unique(labels):
        selected=np.flatnonzero(labels==label);faces=g[selected];used,inverse=np.unique(faces,return_inverse=True);out.append((used,inverse.reshape(-1,3),np.array(kept)[selected]))
    return out,dict(exact_zero_area_source_faces=np.flatnonzero(zero).tolist(),same_orientation_duplicate_source_faces=duplicates,opposite_orientations_cancelled=0,source_oriented_nonzero_surface_set_preserved=True,source_coordinates_moved=0,fill_or_cap_faces_added=0,topological_edge_connected_components=len(out))


def main():
    began=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'classifier-normalization-v26';assert not O.exists();O.mkdir()
    context=json.loads((D/'registered-context-v25/evidence.json').read_bytes());assert context['complete']
    # These known real crossings remain failures, independently of classifiers.
    retained_crossings=context['actual_crossings'];assert len(retained_crossings)==2
    prior=json.loads((D/'native-classifiers-v20/evidence.json').read_bytes());pending=json.loads((D/'conservative-envelopes-v21/evidence.json').read_bytes())['unique_unresolved_physical_objects'];binding=json.loads((D/'source-binding-v19/evidence.json').read_bytes());source=ROOT/binding['snapshot'];expected=binding['snapshot_sha256'];assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();wall=Mesh.load(D/'finite-wall.npz');origin=wall.bounds.mean(0)
    import importlib.util
    spec=importlib.util.spec_from_file_location('normalize_context_v26',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());parent,_=helper.occupancy_classifier(Mesh.load(ROOT/plan['fixed_parent']['mesh']),True);assert parent is not None;parent_solid=helper.solid(parent,origin);occupied={}
    for route in plan['routes']:
        z=np.load(D/(route['name']+'-fields.npz'));m=Mesh(z['fields'][0],z['faces'],wall.frame);c,_=helper.occupancy_classifier(m,True);assert c is not None;occupied[route['name']]=helper.solid(c,origin)-parent_solid
    # A fixture verifies that exact duplicate removal and a vertex-only pinch
    # do not add a surface, change any coordinate or join distinct solids.
    tetra=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float);tf=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]]);tv=np.r_[tetra,-tetra[1:]];ff=np.r_[tf,np.array([[0,4,5],[0,6,4],[0,5,6],[4,6,5]]),tf[:1]];fixture,proof=exact_components(tv,ff);assert len(fixture)==2 and len(proof['same_orientation_duplicate_source_faces'])==1;assert all(topology(Mesh(tv[u],f,wall.frame))['closed_oriented_manifold']for u,f,_ in fixture);atomic_json(O/'exact-operation-fixture.json',proof)
    rows=[];complete=True
    for item in prior['rows']:
        if item['source_object_id']not in pending:continue
        if time.monotonic()-began>45:complete=False;break
        meshes=[];offset=0;vv=[];ff=[]
        for p in item['parts']:
            m=reader.read_part(atlas,parts[p['part_id']],frame=wall.frame,catalogue_path=source,catalogue_sha=expected);assert m.geometry_sha==p['geometry_sha256'];vv.append(m.vertices);ff.append(m.faces+offset);offset+=len(m.vertices);meshes.append(dict(part_id=p['part_id'],geometry_sha256=m.geometry_sha))
        raw=np.concatenate(vv);v,inv=np.unique(raw,axis=0,return_inverse=True);f=inv[np.concatenate(ff)];components,normalization=exact_components(v,f);checks=[];allclosed=True
        for ci,(used,faces,source_faces)in enumerate(components):
            if time.monotonic()-began>48:complete=False;allclosed=False;break
            mesh=Mesh(v[used],faces,wall.frame);gate=topology(mesh);row=dict(component=ci,geometry_sha256=mesh.geometry_sha,topology=gate,source_face_indices=source_faces.tolist(),source_triangle_corner_coordinates_literal=bool(np.array_equal(mesh.triangles,v[f[source_faces]])))
            if gate['closed_oriented_manifold']:
                result=audit(mesh,require_closed=True,max_pairs=1000000,max_tests=12000000);row['audit']=result
                if result['passed']:
                    file=O/('object-'+str(len(rows))+'-component-'+str(ci)+'.npz');mesh.save(file);solid=helper.solid(mesh,origin);row.update(classifier=str(file.relative_to(ROOT)),classifier_sha256=digest_file(file),overlaps=[])
                    for name in sorted(set(r['branch']for r in item['affected_branch_pairs'])):
                        overlap=occupied[name]^solid;assert overlap.status()==md.Error.NoError;row['overlaps'].append(dict(branch=name,overlap_mm3=float(overlap.volume()),empty=bool(overlap.is_empty())))
                else:allclosed=False
            else:allclosed=False
            checks.append(row)
        rows.append(dict(source_object_id=item['source_object_id'],source_members=meshes,normalization=normalization,components=checks,entire_physical_source_has_valid_complete_components=allclosed,partial_closed_components_do_not_resolve_remaining_open_components=True));atomic_json(O/'progress.json',dict(rows=rows,elapsed_seconds=time.monotonic()-began),overwrite=True)
        if not complete:break
    unresolved=[r['source_object_id']for r in rows if not r['entire_physical_source_has_valid_complete_components']];result=dict(schema='exact-classifier-surface-normalization/v26',source_catalogue_sha256=expected,candidate_geometry_sha256=wall.geometry_sha,complete=complete,rows=rows,exact_original_surface_set_preserved=True,whole_model_source_arrays_unchanged=True,classifier_only=True,real_v25_crossings_retained=retained_crossings,unresolved_physical_sources=unresolved,oral_pharyngeal_lumen_exclusion_pending=True,source_aware_gate_pass=False,packing_performed=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);assert digest_file(source)==expected and Mesh.load(D/'finite-wall.npz').geometry_sha==wall.geometry_sha;print(json.dumps({k:v for k,v in result.items()if k!='rows'}))


if __name__=='__main__':main()
