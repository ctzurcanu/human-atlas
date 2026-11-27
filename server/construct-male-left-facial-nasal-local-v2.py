"""One bounded local side-opening trial; does not replace any live source."""
from pathlib import Path
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import manifold3d as md
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.finite_tissue_boolean import make_solid, provenance, CommonPositionBank
from geometry_tools.algorithms.tube_continuity import boundary_loops, finite_wall_cells, annular_rim, continuity_audit
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-nasal-v1'
OUT=BASE/'local-wall-v2'
EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'


def main():
    start=time.monotonic()
    assert digest_file(SOURCE)==EXPECTED
    assert not OUT.exists()
    OUT.mkdir()
    original=json.loads((ROOT/'reports/ta98-male/facial-arterial-family-contract-v1.json').read_bytes())
    row=next(x for x in original['fixed_sources'] if x['role']=='angular-left')
    mesh=Mesh.load(BASE/'angular-left.npz')
    components=row['source_station_components']
    fields=[];carriers=[];rings_all=[];source_faces=[]
    # Only a physically contiguous local candidate is reconstructed. The other
    # four angular components and both far parent stretches are retained outside
    # this packet and are not falsely counted as part of a complete replacement.
    for ci in [0,1,2]:
        rings=np.asarray(components[ci]['rings'],int)
        if ci==0:rings=rings[180:211]
        ids=np.unique(rings)
        selected=np.isin(mesh.faces,ids).all(1)
        if ci:selected &= np.isin(np.arange(len(mesh.faces)),components[ci]['source_faces'])
        oldfaces=np.flatnonzero(selected)
        remap=np.full(len(mesh.vertices),-1,int);remap[ids]=np.arange(len(ids))
        faces=remap[mesh.faces[oldfaces]]
        localrings=remap[rings]
        outer=mesh.vertices[ids]
        inner=outer.copy()
        centers=outer[localrings].mean(1)
        for ring,center in zip(localrings,centers):
            d=outer[ring]-center
            inner[ring]-=0.000045*d/np.linalg.norm(d,axis=1)[:,None]
        field=np.array([outer,inner])
        cell=finite_wall_cells(field,faces)
        atomic_json(OUT/f'component-{ci}-prisms.json',cell)
        atomic_npz(OUT/f'component-{ci}-fields.npz',fields=field,faces=faces,rings=localrings,source_face_ids=oldfaces,source_vertex_ids=ids)
        assert cell['finite_wall_cells_pass'],(ci,cell)
        fields.append(field);carriers.append(faces);rings_all.append(localrings);source_faces.append(oldfaces)
    ports=[(0,0),(0,len(rings_all[0])-1),(1,len(rings_all[1])-1),(2,len(rings_all[2])-1)]
    plan=dict(source_catalogue_sha256=EXPECTED,term_candidate='A12.2.05.028',
              local_source_components=[0,1,2],parent_station_interval=[180,210],
              native_outer_coordinates_fixed=True,illustrative_wall_thickness_m=0.000045,
              source_seven_component_carrier_relabelled=False,partial_local_candidate_only=True,
              septal_term='A12.2.05.027',septal_status='Not constructed: left anterior septal mucoperichondrial/airway boundary absent.',
              context_unmodified=True,installed=False,installable=False)
    atomic_json(OUT/'construction-plan.json',plan)
    origin=fields[0][0].mean(0)
    sheets=[];ancestries=[];classifiers=[]
    for boundary in [0,1]:
        solids=[];ownership={}
        for ci,(field,faces,rings) in enumerate(zip(fields,carriers,rings_all)):
            vertices=field[boundary].copy();triangles=faces.copy();caps={}
            for loop in boundary_loops(faces):
                index=len(vertices);vertices=np.r_[vertices,vertices[loop].mean(0)[None]]
                cap=np.array([[index,int(loop[(j+1)%len(loop)]),int(loop[j])] for j in range(len(loop))])
                external=any(ci==pc and set(loop)==set(rings[ps]) for pc,ps in ports)
                for fi in range(len(triangles),len(triangles)+len(cap)):caps[fi]=external
                triangles=np.r_[triangles,cap]
            solid,identifier=make_solid((vertices-origin)*1000,triangles)
            solids.append(solid);ownership[identifier]=dict(component=ci,caps=caps)
        joined=md.Manifold.batch_boolean(solids,md.OpType.Add)
        assert joined.status()==md.Error.NoError
        raw=joined.to_mesh64()
        vertices=np.asarray(raw.vert_properties)[:,:3]/1000+origin
        triangles=np.asarray(raw.tri_verts,int)
        owners,faceids=provenance(raw)
        external=np.array([ownership[int(o)]['caps'].get(int(f)) is True for o,f in zip(owners,faceids)])
        internal=np.array([ownership[int(o)]['caps'].get(int(f)) is False for o,f in zip(owners,faceids)])
        report=dict(boundary=boundary,retained_internal_caps=int(internal.sum()),retained_external_cap_faces=int(external.sum()))
        atomic_json(OUT/f'boundary-{boundary}-caps.json',report)
        assert not internal.any(),report
        bank=CommonPositionBank(np.concatenate([x[boundary] for x in fields]),tolerance=2e-12)
        mapping=np.array([bank.insert(p) for p in vertices])
        sheet=Mesh(bank.positions,mapping[triangles[~external]],mesh.frame)
        sheet.save(OUT/f'boundary-{boundary}-sheet.npz')
        continuity=continuity_audit(sheet.vertices,sheet.faces,expected_mouths=4,expected_components=1)
        atomic_json(OUT/f'boundary-{boundary}-continuity.json',continuity)
        assert continuity['topological_continuity_pass'],continuity
        classifier=Mesh(bank.positions,mapping[triangles],mesh.frame)
        classifier.save(OUT/f'boundary-{boundary}-classifier.npz')
        sheets.append(sheet);classifiers.append(classifier)
        ancestries.append(dict(component=np.array([ownership[int(o)]['component'] for o in owners[~external]]),source_face=faceids[~external]))
    outer,inner=sheets
    offset=len(outer.vertices);vertices=np.r_[outer.vertices,inner.vertices]
    faces=[outer.faces,inner.faces[:,::-1]+offset]
    oloops=boundary_loops(outer.faces);iloops=boundary_loops(inner.faces)
    for ci,station in ports:
        profile=fields[ci][:,rings_all[ci][station]]
        ol=min(oloops,key=lambda x:np.linalg.norm(outer.vertices[x].mean(0)-profile[0].mean(0)))
        il=min(iloops,key=lambda x:np.linalg.norm(inner.vertices[x].mean(0)-profile[1].mean(0)))
        order=np.argmin(np.linalg.norm(outer.vertices[ol][:,None]-profile[0][None],axis=2),axis=1)
        perm=np.argmin(np.linalg.norm(profile[1][order][:,None]-inner.vertices[il][None],axis=2),axis=1)
        assert len(set(perm))==len(il)==12
        assert np.linalg.norm(inner.vertices[il[perm]]-profile[1][order],axis=1).max()<1e-11
        faces.append(annular_rim(ol,il[perm]+offset))
    wall=Mesh(vertices,np.concatenate(faces),mesh.frame)
    wall.save(OUT/'finite-wall.npz')
    atomic_npz(OUT/'surface-ancestry.npz',outer_component=ancestries[0]['component'],outer_source_face=ancestries[0]['source_face'],
               inner_component=ancestries[1]['component'],inner_source_face=ancestries[1]['source_face'])
    wallcheck=audit(wall,require_closed=True,expected_components=1,expected_genus=[3])
    bloodcheck=audit(classifiers[1],require_closed=True,expected_components=1,expected_genus=[0])
    atomic_json(OUT/'wall-audit.json',wallcheck);atomic_json(OUT/'blood-audit.json',bloodcheck)
    assert wallcheck['passed'],wallcheck['failures']
    assert bloodcheck['passed'],bloodcheck['failures']
    contexts=[]
    for name in ['nasal-septal-cartilage','alar-cartilage-left','lateral-septal-process-left','nasal-mucosa']:
        fixed=Mesh.load(BASE/(name+'.npz'))
        arrays,proof=pair_contacts(wall,fixed,0.,max_pairs=20000,max_tests=600000)
        contexts.append(dict(name=name,geometry_sha256=fixed.geometry_sha,proof=proof))
        atomic_npz(OUT/(name+'-witnesses.npz'),**arrays)
    atomic_json(OUT/'fixed-context.json',dict(checks=contexts))
    result=dict(**plan,local_wall_pass=True,blood_continuous=True,open_blood_ports=4,
                self_crossings=wallcheck['self_intersections']['defect_pairs'],
                forbidden_surface_crossings=sum(x['proof']['crossing_face_pairs'] for x in contexts),
                context_complete=all(x['proof']['complete_within_pair_budget'] for x in contexts),
                recipient_finite_wall_ownership_verified=False,airway_exclusion_verified=False,
                emitted_float32_bank_checked=False,full_carrier_ancestry_partition_checked=False,
                actual_native_render_reviewed=False,geometry_sha256=wall.geometry_sha,elapsed_seconds=time.monotonic()-start)
    atomic_json(OUT/'gates.json',result)
    assert digest_file(SOURCE)==EXPECTED
    print(json.dumps(result))


if __name__=='__main__':main()
