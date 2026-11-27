"""One finite left lingual family with an actual native ECA side opening."""
from pathlib import Path
import sys,time,json,argparse
import numpy as np
import manifold3d as md
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames,path_tangents,centripetal_curve
from geometry_tools.algorithms.tube_continuity import boundary_loops,loft_faces,annular_rim,finite_wall_cells,continuity_audit
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,CommonPositionBank
from geometry_tools.audit import audit
from geometry_tools.contact import self_contact_summary

args=argparse.ArgumentParser();args.add_argument('--version',default='v3');args.add_argument('--side',choices=['left','right'],default='left');args.add_argument('--output-version');args=args.parse_args()
T=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';P=B/('male-'+args.side+'-lingual-routes-'+args.version);plan=json.loads((P/'routes.json').read_text());D=B/('male-'+args.side+'-lingual-finite-network-'+(args.output_version or args.version));D.mkdir(exist_ok=False);eca=Mesh.load(ROOT/plan['fixed_native_carotid']['mesh']);frame=eca.frame;origin=np.array([.025 if args.side=='left'else-.025,1.625,.035]);parent_rings=np.array(plan['fixed_native_carotid']['source_station_components'][0]['rings']);parent_centers=eca.vertices[parent_rings].mean(1);station=plan['source_origin_station'];delta=eca.vertices[parent_rings[station]]-parent_centers[station];_,_,vt=np.linalg.svd(delta,full_matrices=False);normal=vt[-1]
if normal@(parent_centers[station+1]-parent_centers[station-1])<0:normal=-normal
ref=delta[0]-normal*(delta[0]@normal);ref/=np.linalg.norm(ref);basis=np.c_[ref,np.cross(normal,ref),normal];profile=delta@basis;profile[:,2]=0;profile/=np.linalg.norm(profile,axis=1).mean()
routes=plan['routes'];guides=[np.load(r['guide'])for r in routes];paths=[z['points'].copy()for z in guides];radii=[z['outer_radii'].copy()for z in guides];walls=[np.full(len(p),r['illustrative_wall_m'])for p,r in zip(paths,routes)]
# Deep lingual is the continuous terminal course, with one shared true polygon
# at the change of ownership; its start is not a small overlapping blind tube.
main_tangent=path_tangents(paths[0]);paths[-1],_=centripetal_curve([paths[0][-1],paths[0][-1]+main_tangent[-1]*.003,paths[-1][-1]],chord_tolerance_m=.00001,max_step_m=.0005);radii[-1]=np.linspace(radii[0][-1],radii[-1][-1],len(paths[-1]));walls[-1]=np.linspace(walls[0][-1],walls[-1][-1],len(paths[-1]))
fields=[];carriers=[];rings=[];names=['external-carotid-'+args.side]+[r['name']for r in routes];terms=[None]+[r['term']for r in routes];parent_inner=eca.vertices.copy()
for ring,c in zip(parent_rings,parent_centers):
    d=eca.vertices[ring]-c;parent_inner[ring]=eca.vertices[ring]-.00035*d/np.linalg.norm(d,axis=1)[:,None]
fields.append(np.array([eca.vertices,parent_inner]));carriers.append(eca.faces);rings.append(parent_rings)
for j,(path,radius,wall)in enumerate(zip(paths,radii,walls)):
    tangent=path_tangents(path)
    if j==len(paths)-1:tangent[0]=main_tangent[-1]
    frames=rotation_minimizing_frames(path,tangent,initial_reference=ref);section=np.einsum('kij,nj->kni',frames,profile);outer=path[:,None]+radius[:,None,None]*section;inner=path[:,None]+(radius-wall)[:,None,None]*section
    if j==len(paths)-1:
        outer[0]=fields[1][0][rings[1][-1]];inner[0]=fields[1][1][rings[1][-1]]
    f=loft_faces(len(path),12);p=np.array([outer.reshape(-1,3),inner.reshape(-1,3)]);check=finite_wall_cells(p,f);atomic_json(D/(names[j+1]+'-paired-field.json'),check);atomic_npz(D/(names[j+1]+'-fields.npz'),fields=p,faces=f,guide_points=path,outer_radii=radius,inner_radii=radius-wall,frame=np.array(frame));
    assert check['finite_wall_cells_pass'],(names[j+1],check['invalid_counts_by_interval'])
    fields.append(p);carriers.append(f);rings.append(np.arange(len(path)*12).reshape(-1,12))
external=[(0,0),(0,len(rings[0])-1)]+[(i,len(rings[i])-1)for i in range(2,len(fields))];assert len(external)==7
contract=dict(schema='left-lingual-finite-wall-construction/v1',route_contract_sha256=digest_file(P/'routes.json'),source_catalogue_sha256=plan['source_catalogue_sha256'],
    owned_term_ids=sorted(set(x for x in terms if x)),actual_tube_names=names,external_blood_ports=external,
    native_parent_coordinates_fixed=True,native_parent_opening_required=True,all_neck_skeleton_and_tongue_sources_fixed=True,
    shared_lingual_deep_polygon=True,source_profile='Actual source ECA polygon, normalized and calibrated to the primary-study lingual calibre; no private reference surface copied',
    clinically_verified_wall_thickness=False,terminal_perfusion_beds_unmodeled=True,anatomical_acceptance=False,installed=False)
atomic_json(D/'construction-plan.json',contract)
inputs={};sheets=[];ancestry=[];blood=None
for inner in [False,True]:
    solids=[]
    for i,(field,f,rr)in enumerate(zip(fields,carriers,rings)):
        v=field[int(inner)].copy();ff=f.copy();caps=[]
        for k,loop in enumerate(boundary_loops(ff)):
            # Opposed shared classifier polygons must receive bit-identical
            # fan centres, independent of their cyclic traversal direction.
            polygon=v[loop];order=np.lexsort((polygon[:,2],polygon[:,1],polygon[:,0]));point=len(v);v=np.r_[v,polygon[order].mean(0)[None]];new=np.array([[point,int(loop[(j+1)%len(loop)]),int(loop[j])]for j in range(len(loop))]);ids=np.arange(len(ff),len(ff)+len(new));ff=np.r_[ff,new];is_external=any(i==ci and set(loop)==set(rr[si])for ci,si in external);caps.extend((int(fi),is_external)for fi in ids)
        solid,oid=make_solid((v-origin)*1000,ff);assert solid.status()==md.Error.NoError;inputs[oid]=dict(component=i,source_faces=len(f),caps=dict(caps));solids.append(solid)
    joined=md.Manifold.batch_boolean(solids,md.OpType.Add);assert joined.status()==md.Error.NoError;raw=joined.to_mesh64();v=np.asarray(raw.vert_properties)[:,:3]/1000+origin;f=np.asarray(raw.tri_verts,int);owner,sourceface=provenance(raw);external_cap=np.array([inputs[int(o)]['caps'].get(int(fi))is True for o,fi in zip(owner,sourceface)]);internal_cap=np.array([inputs[int(o)]['caps'].get(int(fi))is False for o,fi in zip(owner,sourceface)]);atomic_json(D/(('inner'if inner else'outer')+'-classifier-caps.json'),dict(retained_internal_caps=int(internal_cap.sum()),retained_internal_by_component={names[c]:int(sum(inputs[int(o)]['component']==c for o in owner[internal_cap]))for c in range(len(names))},retained_external_caps=int(external_cap.sum())));assert not internal_cap.any(),('retained_internal_classifier',inner,int(internal_cap.sum()))
    bank=CommonPositionBank(np.concatenate([p[int(inner)]for p in fields]),tolerance=2e-12);mapping=np.array([bank.insert(point)for point in v]);sheet=Mesh(bank.positions,mapping[f[~external_cap]],frame);label='inner'if inner else'outer';sheet.save(D/(label+'-sheet.npz'));check=continuity_audit(sheet.vertices,sheet.faces,expected_mouths=7,expected_components=1);atomic_json(D/(label+'-continuity.json'),check);assert check['topological_continuity_pass'];sheets.append(sheet);ancestry.append(dict(component=np.array([inputs[int(o)]['component']for o in owner[~external_cap]]),sourceface=sourceface[~external_cap]))
    if inner:blood=Mesh(bank.positions,mapping[f],frame);blood.save(D/'closed-blood-classifier.npz')
outer,inner=sheets;offset=len(outer.vertices);vertices=np.r_[outer.vertices,inner.vertices];faces=[outer.faces,inner.faces[:,::-1]+offset];ol=boundary_loops(outer.faces);il=boundary_loops(inner.faces);ports=[]
for component,station in external:
    expected=fields[component][:,rings[component][station]];center=expected[0].mean(0);first=min(ol,key=lambda r:np.linalg.norm(outer.vertices[r].mean(0)-center));second=min(il,key=lambda r:np.linalg.norm(inner.vertices[r].mean(0)-expected[1].mean(0)));order=np.argmin(np.linalg.norm(expected[0][None]-outer.vertices[first][:,None],axis=2),axis=1);target=expected[1][order];permutation=np.argmin(np.linalg.norm(inner.vertices[second][None]-target[:,None],axis=2),axis=1);assert len(set(permutation))==len(second)==12;assert np.linalg.norm(inner.vertices[second[permutation]]-target,axis=1).max()<1e-11;faces.append(annular_rim(first,second[permutation]+offset));ports.append(dict(component=component,source_station=station,outer_indices=first.tolist(),inner_indices=(second[permutation]+offset).tolist(),terminal_inside_recipient=component>0))
wall=Mesh(vertices,np.concatenate(faces),frame);wall.save(D/'finite-wall.npz');check=audit(wall,require_closed=True,expected_components=1,expected_genus=[6],max_pairs=200000,max_tests=2000000);atomic_json(D/'wall-audit.json',check);assert check['passed'],check['failures'];bloodcheck=audit(blood,require_closed=True,expected_components=1,expected_genus=[0]);atomic_json(D/'blood-audit.json',bloodcheck);assert bloodcheck['passed']
atomic_npz(D/'surfaces-and-ancestry.npz',vertices=vertices,outer_faces=outer.faces,inner_faces=inner.faces[:,::-1]+offset,
    outer_component=ancestry[0]['component'],outer_source_face=ancestry[0]['sourceface'],inner_component=ancestry[1]['component'],inner_source_face=ancestry[1]['sourceface'],frame=np.array(frame),**{f'annulus_{i}_faces':ff for i,ff in enumerate(faces[2:])})
atomic_json(D/'ports.json',dict(ports=ports,component_names=names,component_terms=terms));report=dict(schema='left-lingual-five-term-finite-network/v1',source_catalogue_sha256=plan['source_catalogue_sha256'],geometry_sha256=wall.geometry_sha,term_ids=contract['owned_term_ids'],actual_branch_tubes=6,
    finite_wall_passed=True,wall_faces=len(wall.faces),blood_volume_connected=True,external_blood_ports=7,internal_caps=0,lumen_disks=0,
    native_ECA_parent_coordinates_literal=True,source_profile_preserved_by_proper_section_transport=True,ownership_partition_pending=True,
    full_tissue_oral_cavity_check_pending=True,anatomical_acceptance=False,installable=False,elapsed_seconds=time.monotonic()-T)
atomic_json(D/'evidence.json',report);print(json.dumps(report))
