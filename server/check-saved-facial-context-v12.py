"""Recover verified source-face lineage, then check saved native context.

Candidate vertices/faces are copied literally. No Boolean or anatomical
geometry is rebuilt. An output face must lie within its original carrier
triangle before it receives a component/source-face label.
"""
from pathlib import Path
import sys,json,time,runpy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim


def ancestry(sheet,fields,carriers,boundary):
    points=[];triangles=[];components=[];local=[];offset=0
    for component,(field,face)in enumerate(zip(fields,carriers)):
        points.append(field[boundary]);triangles.append(face+offset);offset+=len(field[boundary]);components.extend([component]*len(face));local.extend(range(len(face)))
    points=np.concatenate(points);triangles=np.concatenate(triangles);centers=sheet.triangles.mean(1)
    d2,ids,_=igl.point_mesh_squared_distance(centers,points,triangles);source=points[triangles[ids]];q=sheet.triangles-source[:,0,None]
    first,second=source[:,1]-source[:,0],source[:,2]-source[:,0];a=(first*first).sum(1);b=(first*second).sum(1);c=(second*second).sum(1);denominator=a*c-b*b
    assert np.all(denominator>0),'Degenerate original carrier triangle'
    d=np.einsum('fvi,fi->fv',q,first);e=np.einsum('fvi,fi->fv',q,second)
    u=(c[:,None]*d-b[:,None]*e)/denominator[:,None];v=(a[:,None]*e-b[:,None]*d)/denominator[:,None];weights=np.stack([1-u-v,u,v],axis=-1)
    projected=np.einsum('fvk,fki->fvi',weights,source);error=np.linalg.norm(projected-sheet.triangles,axis=2).max(1)
    invalid=(error>8e-12)|(weights.min((1,2))<-.000003)
    proof=dict(output_faces=len(ids),maximum_centroid_surface_distance_m=float(np.sqrt(d2.max())),maximum_full_triangle_projection_error_m=float(error.max()),minimum_barycentric_weight=float(weights.min()),invalid_source_face_witnesses=int(invalid.sum()))
    if invalid.any():
        proof['failed_face_examples']=np.flatnonzero(invalid)[:20].tolist();return None,None,proof
    return np.asarray(components,int)[ids],np.asarray(local,int)[ids],proof


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';directory=base/'male-left-facial-cervical-finite-network-free-c2-v10';route_dir=base/'male-left-facial-cervical-routes-free-c2-v10'
    gate=json.loads((directory/'self-contact-budget-v11/evidence.json').read_bytes());assert gate['complete_wall_and_blood_audits_pass'];plan=json.loads((route_dir/'routes.json').read_bytes());source=ROOT/plan['source_catalogue'];expected=plan['source_catalogue_sha256'];assert digest_file(source)==expected
    assert not(directory/'surfaces-and-ancestry.npz').exists()
    wall=Mesh.load(directory/'finite-wall.npz');outer=Mesh.load(directory/'outer-sheet.npz');inner=Mesh.load(directory/'inner-sheet.npz');n=len(outer.vertices)
    assert np.array_equal(wall.vertices,np.r_[outer.vertices,inner.vertices]);assert np.array_equal(wall.faces[:len(outer.faces)],outer.faces)
    assert np.array_equal(wall.faces[len(outer.faces):len(outer.faces)+len(inner.faces)],inner.faces[:,::-1]+n)
    parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);rings=np.array(plan['fixed_parent']['source_station_components'][0]['rings']);centers=parent.vertices[rings].mean(1);native_inner=parent.vertices.copy()
    for ring,center in zip(rings,centers):
        delta=parent.vertices[ring]-center;native_inner[ring]-=plan['illustrative_parent_wall_m']*delta/np.linalg.norm(delta,axis=1)[:,None]
    fields=[np.array([parent.vertices,native_inner])];carriers=[parent.faces];allrings=[rings]
    for route in plan['routes']:
        z=np.load(directory/(route['name']+'-fields.npz'));fields.append(z['fields']);carriers.append(z['faces']);allrings.append(np.arange(len(z['fields'][0])).reshape(-1,12))
    oc,of,op=ancestry(outer,fields,carriers,0);ic,iff,ip=ancestry(inner,fields,carriers,1)
    proof=dict(source_catalogue_sha256=expected,candidate_geometry_sha256=wall.geometry_sha,geometry_mutated=False,outer_original_carrier_witnesses=op,inner_original_carrier_witnesses=ip,
               parent_inner_field_scope='Exact replay of declared illustrative radial offset only for original-carrier witness validation; saved candidate geometry is not reconstructed',elapsed_seconds=time.monotonic()-began)
    if oc is None or ic is None:
        atomic_json(directory/'lineage-recovery-evidence-v12.json',proof);print(json.dumps(proof));return
    ports=[(0,0),(0,len(rings)-1)]+[(i,len(allrings[i])-1)for i in range(1,len(fields))]
    outer_loops=boundary_loops(outer.faces);inner_loops=boundary_loops(inner.faces);annuli=[];portrows=[]
    for component,station in ports:
        expected_section=fields[component][:,allrings[component][station]]
        first=min(outer_loops,key=lambda loop:np.linalg.norm(outer.vertices[loop].mean(0)-expected_section[0].mean(0)))
        second=min(inner_loops,key=lambda loop:np.linalg.norm(inner.vertices[loop].mean(0)-expected_section[1].mean(0)))
        order=np.argmin(np.linalg.norm(expected_section[0][None]-outer.vertices[first][:,None],axis=2),axis=1);target=expected_section[1][order]
        permutation=np.argmin(np.linalg.norm(inner.vertices[second][None]-target[:,None],axis=2),axis=1)
        assert len(set(permutation))==len(second)==12 and np.linalg.norm(inner.vertices[second[permutation]]-target,axis=1).max()<1e-11
        annuli.append(annular_rim(first,second[permutation]+n));portrows.append(dict(component=component,source_station=station,outer_indices=first.tolist(),inner_indices=(second[permutation]+n).tolist(),recipient=plan['routes'][component-1]['recipient']if component else None))
    assert np.array_equal(np.concatenate(annuli),wall.faces[len(outer.faces)+len(inner.faces):])
    atomic_npz(directory/'surfaces-and-ancestry.npz',vertices=wall.vertices,outer_faces=outer.faces,inner_faces=inner.faces[:,::-1]+n,outer_component=oc,outer_source_face=of,inner_component=ic,inner_source_face=iff,frame=np.array(wall.frame),**{f'annulus_{i}_faces':a for i,a in enumerate(annuli)})
    atomic_json(directory/'ports.json',dict(ports=portrows,component_names=['facial-left']+[r['name']for r in plan['routes']],component_terms=['A12.2.05.020']+[r['term']for r in plan['routes']]))
    proof.update(all_face_lineage_verified=True,saved_annular_faces_literal=True,metadata_only_recovery=True,elapsed_seconds=time.monotonic()-began);atomic_json(directory/'lineage-recovery-evidence-v12.json',proof)
    remaining=min(48.,max(1.,55-(time.monotonic()-began)))
    sys.argv=['check-facial-cervical-feasibility.py','--route-version','free-c2-v10','--network-version','free-c2-v10','--occupied-volume','--work-seconds',str(remaining)]
    runpy.run_path(str(ROOT/'server/check-facial-cervical-feasibility.py'),run_name='__main__')
    result=json.loads((directory/'feasibility-v2/evidence.json').read_bytes());atomic_json(directory/'native-context-stage-v12.json',dict(source_catalogue_sha256=expected,candidate_geometry_sha256=wall.geometry_sha,source_face_lineage_verified=True,context_evidence='feasibility-v2/evidence.json',numerical_feasibility_pass=result['numerical_feasibility_pass'],geometry_mutated=False,installed=False,blender_queued=False,elapsed_seconds=time.monotonic()-began))
    assert digest_file(source)==expected


if __name__=='__main__':main()
