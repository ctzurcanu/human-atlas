"""Native SG/SP coarse corridor analysis only: no guide, tube or takeoff edit."""
from pathlib import Path
import sys,json,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader


def main():
    began=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'global-corridor-analysis-v31';assert not O.exists();O.mkdir()
    binding=json.loads((D/'source-binding-v19/evidence.json').read_bytes());source=ROOT/binding['snapshot'];expected=binding['snapshot_sha256'];assert digest_file(source)==expected;atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());frame=plan['frame'];proof=json.loads((B/'male-left-facial-cervical-routes-palatine-active-c2-v28/literal-source-rebase.json').read_bytes())
    spec=importlib.util.spec_from_file_location('global_portal_occupancy_v31',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    fixed={}
    for row in proof['fixed_blocker_geometry']:
        pid=row['part_id']
        if pid in parts:m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=expected);allow=parts[pid].get('system')in{'arterial','venous'}
        else:
            extra=next(r for r in plan['extra_fixed_meshes']if r['part_id']==pid);assert digest_file(ROOT/extra['mesh'])==extra['mesh_sha256'];m=Mesh.load(ROOT/extra['mesh']);allow=True
        assert m.geometry_sha==row['geometry_sha256'];classifier,mode=helper.occupancy_classifier(m,allow);fixed[pid]=(m,classifier,mode)
    sg=fixed['TA98TONGUE:male:styloglossus:left'][0];sp=fixed['ZA:Stylopharyngeus muscle.l'][0]
    # A new global domain comes from the actual muscle extents, not offsets
    # around the rejected v17 guide or dense private-reference coordinates.
    low=np.maximum(sg.bounds[0],sp.bounds[0])-.002;high=np.minimum(sg.bounds[1],sp.bounds[1])+.002;step=.001
    axes=[np.arange(low[k],high[k]+step/2,step)for k in range(3)];q=np.stack(np.meshgrid(*axes,indexing='ij'),-1).reshape(-1,3)
    da,fa,pa=igl.point_mesh_squared_distance(q,sg.vertices,sg.faces);db,fb,pb=igl.point_mesh_squared_distance(q,sp.vertices,sp.faces);da=np.sqrt(da);db=np.sqrt(db);na=(q-pa)/np.maximum(da[:,None],1e-30);nb=(q-pb)/np.maximum(db[:,None],1e-30);opposed=np.einsum('ij,ij->i',na,nb)<-.3
    selected=np.flatnonzero((da<.004)&(db<.004)&opposed);points=q[selected];clearance=np.full(len(points),np.inf);safe=np.ones(len(points),bool);radius=.0005501;margin=.00007;source_rows=[];complete=True
    for pid,(mesh,classifier,mode)in fixed.items():
        if time.monotonic()-began>45:complete=False;break
        if pid in {'ZA:Palatopharyngeus muscle.l','ZA:Soft palate'}:source_rows.append(dict(part_id=pid,declared_intended_recipient_interface=True));continue
        delta=np.maximum(mesh.bounds[0]-points,0)+np.maximum(points-mesh.bounds[1],0);box=np.linalg.norm(delta,axis=1);take=np.flatnonzero(box<.004);clearance=np.minimum(clearance,box-radius)
        if not len(take):source_rows.append(dict(part_id=pid,geometry_sha256=mesh.geometry_sha,queried=0,safely_distant_by_native_AABB=True));continue
        squared,faces,near=igl.point_mesh_squared_distance(points[take],mesh.vertices,mesh.faces);distance=np.sqrt(squared)
        if classifier is not None:inside=abs(igl.winding_number(classifier.vertices,classifier.faces,points[take]))>.99;distance[inside]*=-1
        # Replace the AABB lower bound by the exact value at queried points.
        # Keep an independent exact min for ranking; no candidate is constructed.
        bad=distance<radius+margin;safe[take[bad]]=False;source_rows.append(dict(part_id=pid,geometry_sha256=mesh.geometry_sha,queried=len(take),rejected=int(bad.sum()),classifier=mode,minimum_signed_or_unsigned_distance_m=float(distance.min())))
    portals=[]
    if complete:
        valid=np.flatnonzero(safe)
        for y in axes[1]:
            ids=valid[abs(points[valid,1]-y)<step/4]
            if not len(ids):continue
            # Prefer the center of the observed opposing SG/SP surfaces.
            score=np.minimum(da[selected[ids]],db[selected[ids]]);i=int(ids[np.argmax(score)]);ix=int(selected[i]);witness=[]
            for pid in ['TA98TONGUE:male:styloglossus:left','ZA:Stylopharyngeus muscle.l','ZA:Stylohyoid ligament.l','ZA:Stylomandibular ligament.l','ZA:Superior pharyngeal constrictor.l','ZA:Facial vein.l','ZA:Hypoglossal nerve (XII).l']:
                m,c,mode=fixed[pid];d2,faces,near=igl.point_mesh_squared_distance(points[i:i+1],m.vertices,m.faces);distance=float(np.sqrt(d2[0]));inside=bool(abs(igl.winding_number(c.vertices,c.faces,points[i:i+1])[0])>.99)if c is not None else None
                witness.append(dict(part_id=pid,geometry_sha256=m.geometry_sha,distance_m=distance,inside_valid_closed_source=inside,face=int(faces[0]),source_point_m=near[0].tolist(),classifier=mode))
            portals.append(dict(point_m=points[i].tolist(),opposing_muscle_surface_points_m=[pa[ix].tolist(),pb[ix].tolist()],local_opposing_surface_separation_m=float(np.linalg.norm(pa[ix]-pb[ix])),source_profile_radius_m=radius,minimum_requested_clearance_m=margin,clears_all_93_at_this_discrete_center=True,source_witnesses=witness))
    parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);rings=np.array(plan['fixed_parent']['source_station_components'][0]['rings']);centers=parent.vertices[rings].mean(1);takeoffs=[]
    for station in range(25):
        point=centers[station];witness=[]
        for pid in ['TA98TONGUE:male:styloglossus:left','ZA:Stylopharyngeus muscle.l','ZA:Stylohyoid ligament.l','ZA:Stylomandibular ligament.l']:
            m,_,_=fixed[pid];d2,face,near=igl.point_mesh_squared_distance(point[None],m.vertices,m.faces);witness.append(dict(part_id=pid,distance_m=float(np.sqrt(d2[0])),face=int(face[0]),source_point_m=near[0].tolist()))
        takeoffs.append(dict(native_station=station,point_m=point.tolist(),held_tonsillar_root_station=16,near_held_tonsillar_root=abs(station-16)<=2,native_parent_unchanged=True,source_witnesses=witness))
    result=dict(schema='native-global-palatine-portal-analysis/v31',source_catalogue_sha256=expected,source_registry_revision='31fab36e6897adaec9c8ddca201261607902203d75c64f3b214e57d329f460c2',held_v17_geometry_sha256=Mesh.load(D/'finite-wall.npz').geometry_sha,analysis_complete=complete,domain_basis='Actual SG/SP overlapping source extents plus2mm, one-mm coarse samples, points with opposing nearest surface vectors and within4mm of both actual muscles; not a rejected-guide offset corridor.',domain_bounds_m=[low.tolist(),high.tolist()],global_grid_points=len(q),between_muscle_surface_samples=len(points),safe_discrete_samples=int(safe.sum())if complete else None,all_fixed_objects=len(fixed),fixed_source_queries=source_rows,coarse_portal_landmarks=portals,native_origin_station_measurements=takeoffs,global_centreline_or_mesh_constructed=False,continuous_clearance_or_patency_not_proved=True,open_source_occupied_volumes_remain_unknown=True,complete_oral_fluid_exclusion_unresolved=True,source_points_moved=False,new_takeoff_selected=False,primary_reference_interpretation=[dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC5951712/',finding='Cadaveric Fig4D supports APA in prestyloid interval between SG/SP. Text places ECA between SG/SH and stylomandibular ligament overSH/partofECA; this supports approaching the interval through another side of the finite styloid apparatus rather than using a local ligament-crossing correction.'),dict(url='https://pubmed.ncbi.nlm.nih.gov/24510490/',qualification='Abstract uses APA for ascending palatine and aPA for ascending pharyngeal; the12.6±3.9mm sentence literally says aPA. Do not impose this as a palatine target without resolving fulltext abbreviation.'),dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11298333/',finding='Cadaveric Fig2 describes APA from proximal facial course crossingSG/SP towardpharynx; crosses doesnotalone prove muscle perforation or authorize changing fixed muscle meshes.')],installed=False,packing_performed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);assert digest_file(source)==expected;print(json.dumps({k:result[k]for k in ['analysis_complete','global_grid_points','between_muscle_surface_samples','safe_discrete_samples','all_fixed_objects','global_centreline_or_mesh_constructed','elapsed_seconds']}));print(json.dumps(dict(portal_count=len(portals),portal_y_range_m=[portals[0]['point_m'][1],portals[-1]['point_m'][1]]if portals else None)))


if __name__=='__main__':main()
