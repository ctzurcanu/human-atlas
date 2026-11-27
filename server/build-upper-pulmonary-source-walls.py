"""Construct the two native left anterior arterial trees before naming children.

Source connectivity and coincident station centres define the branch graph.
No anatomical branch-name assignment or catalogue installation is performed.
The radial wall fraction is an illustrative finite wall, not histology.
"""
from pathlib import Path
import json
import sys
import time

import manifold3d as md
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,CommonPositionBank
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim,continuity_audit

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-upper-pulmonary-subbranches-v1'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-anterior-arterial-walls-v1'
ORIGIN=np.array([.0625,1.4375,.03125])
FRACTION=.12


def construct(source,profiles,components,joints,external,folder):
    folder.mkdir()
    v0,f0=source.vertices,source.faces
    rings={i:profiles[i]['rings']for i in components}
    centres={i:v0[rings[i]].mean(1)for i in components}
    native_faces={i:profiles[i]['source_face_ids']for i in components}
    for child,station,parent,where in joints:
        assert np.linalg.norm(centres[child][station]-centres[parent][where])<1e-7
    layers=[];proof=[]
    for inner in [False,True]:
        native=v0.copy()
        if inner:
            for i in components:
                for ring,centre in zip(rings[i],centres[i]):
                    native[ring]=centre+(v0[ring]-centre)*(1-FRACTION)
        inputs={};solids=[]
        for i in components:
            v=native.copy();f=f0[native_faces[i]].copy();count=len(f)
            temporary_external=[];temporary_internal=[]
            root=next((station for child,station,parent,where in joints if child==i),None)
            for loop in boundary_loops(f):
                centre=len(v);v=np.r_[v,v[loop].mean(0)[None]]
                ids=np.arange(len(f),len(f)+len(loop))
                f=np.r_[f,[[centre,int(loop[(j+1)%len(loop)]),int(loop[j])]for j in range(len(loop))]]
                internal=root is not None and set(loop)==set(rings[i][root])
                (temporary_internal if internal else temporary_external).extend(ids.tolist())
            solid,oid=make_solid((v-ORIGIN)*1000,f)
            assert solid.status()==md.Error.NoError
            inputs[oid]=dict(component=i,source_count=count,internal=temporary_internal,external=temporary_external)
            solids.append(solid)
        joined=md.Manifold.batch_boolean(solids,md.OpType.Add)
        assert joined.status()==md.Error.NoError
        raw=joined.to_mesh64();v=np.asarray(raw.vert_properties)[:,:3]/1000+ORIGIN
        f=np.asarray(raw.tri_verts,int);owners,ancestry=provenance(raw)
        internal=np.array([int(face)in inputs[int(owner)]['internal']for owner,face in zip(owners,ancestry)])
        external_caps=np.array([int(face)in inputs[int(owner)]['external']for owner,face in zip(owners,ancestry)])
        cap_result=dict(inner=inner,exposed_internal_classifier_faces=int(internal.sum()))
        proof.append(cap_result);atomic_json(folder/'classifier-evidence.json',proof,overwrite=True)
        assert not internal.any(),'A disconnected origin has an exposed temporary classifier cap'
        bank=CommonPositionBank(native,tolerance=2e-12)
        mapping=np.array([bank.insert(point)for point in v])
        sheet=Mesh(bank.positions,mapping[f[~external_caps]],source.frame)
        check=continuity_audit(sheet.vertices,sheet.faces,expected_mouths=len(external),expected_components=1)
        atomic_json(folder/('inner-continuity.json'if inner else'outer-continuity.json'),check)
        assert check['topological_continuity_pass']
        sheet.save(folder/('inner.npz'if inner else'outer.npz'))
        layers.append(sheet)
    outer,inner=layers;offset=len(outer.vertices)
    faces=[outer.faces,inner.faces[:,::-1]+offset]
    outer_loops=boundary_loops(outer.faces);inner_loops=boundary_loops(inner.faces);ports=[]
    for component,station in external:
        expected=v0[rings[component][station]];centre=expected.mean(0)
        a=min(outer_loops,key=lambda loop:np.linalg.norm(outer.vertices[loop].mean(0)-centre))
        b=min(inner_loops,key=lambda loop:np.linalg.norm(inner.vertices[loop].mean(0)-centre))
        target=centre+(outer.vertices[a]-centre)*(1-FRACTION)
        order=np.argmin(np.linalg.norm(inner.vertices[b][None]-target[:,None],axis=2),axis=1)
        assert len(a)==len(b)==len(set(order))==12
        assert np.linalg.norm(inner.vertices[b[order]]-target,axis=1).max()<1e-11
        assert all(np.any(np.all(outer.vertices[a]==point,axis=1))for point in expected)
        faces.append(annular_rim(a,b[order]+offset))
        ports.append(dict(component=component,station=station,literal_native_outer_controls=True))
    wall=Mesh(np.r_[outer.vertices,inner.vertices],np.concatenate(faces),source.frame)
    wall.save(folder/'finite-wall.npz')
    check=audit(wall,require_closed=True,expected_components=1,expected_genus=[len(external)-1])
    atomic_json(folder/'audit.json',check)
    assert check['passed'],check['failures']
    original={tuple(sorted(map(tuple,triangle)))for triangle in source.triangles}
    return dict(components=components,joints=joints,external_ports=ports,faces=len(wall.faces),
                literal_outer_faces=sum(tuple(sorted(map(tuple,triangle)))in original for triangle in outer.triangles),
                finite_wall_geometry_sha256=wall.geometry_sha,raw_audit_passed=True)


def main():
    started=time.monotonic()
    source=Mesh.load(BASE/'left-anterior-source.npz')
    inspection=json.loads((BASE/'inspection.json').read_bytes())
    assert inspection['source_catalogue_sha256']==digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')
    assert not OUT.exists(),'Retain previous candidate evidence'
    OUT.mkdir()
    profiles={i:dict(np.load(BASE/f'left-anterior-component-{i}.npz',allow_pickle=False))for i in range(6)}
    plan=dict(source_geometry_sha256=source.geometry_sha,source_arrays=[digest_file(BASE/f'left-anterior-component-{i}.npz')for i in range(6)],
              source_catalogue_sha256=inspection['source_catalogue_sha256'],inventory='TA98-CIRCULATION.md',
              planned_authority_terms=['A12.2.01.205','A12.2.01.206','A12.2.01.207'],
              adjoining_interfaces=['two native left pulmonary arterial origins','six native distal source profiles'],
              all_outside_anatomy_fixed=True,wall_fraction=FRACTION,histological_thickness_proven=False,
              branch_identity_assignment_accepted=False,installed=False)
    atomic_json(OUT/'plan.json',plan)
    first=construct(source,profiles,[0,1],[(1,0,0,60)],[(0,0),(0,72),(1,24)],OUT/'native-tree-0')
    second=construct(source,profiles,[2,3,4,5],[(3,0,2,60),(4,0,2,24),(5,0,4,24)],
                     [(2,0),(2,84),(3,24),(4,48),(5,36)],OUT/'native-tree-2')
    result=dict(schema='native-left-anterior-arterial-finite-wall-candidate/v1',trees=[first,second],
                elapsed_seconds=time.monotonic()-started,anatomy_accepted=False,installed=False,
                remaining=['actual source/branch identity review','fixed-host context','Float32 bank replay','native render and section','catalogue installation'])
    atomic_json(OUT/'construction.json',result)
    print(json.dumps(result))


if __name__=='__main__':
    main()
