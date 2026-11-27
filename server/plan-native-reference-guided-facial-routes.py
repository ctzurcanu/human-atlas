"""One bounded native corridor series; private reference topology only.

The rejected 10.97-mm-RMS registration is never applied to route controls.
Actual fixed native muscles/veins define an occupied 3-D search domain.
"""
from pathlib import Path
import sys,json,time,heapq,importlib.util,itertools
from copy import deepcopy
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter1d
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.tube_sections import centripetal_curve

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion'
OUT=BASE/'male-left-facial-cervical-routes-native-domain-v7'


def lattice_path(safe,start,stop,shape,coordinates,step,began):
    tree=cKDTree(coordinates[safe]);ids=np.flatnonzero(safe)
    gaps,nearest=tree.query(np.array([start,stop]))
    if gaps.max()>.0025:
        return None,dict(reason='No source-calibre-safe endpoint connector within 2.5 mm',endpoint_gaps_m=gaps.tolist())
    first,last=ids[nearest]
    dims=np.array(shape);stride=np.array([shape[1]*shape[2],shape[2],1]);neighbors=[(np.array(d),np.linalg.norm(d)*step)for d in itertools.product([-1,0,1],repeat=3)if d!=(0,0,0)]
    distance={int(first):0.};previous={};queue=[(float(np.linalg.norm(coordinates[first]-coordinates[last])),0.,int(first))];visited=0
    while queue:
        _,cost,node=heapq.heappop(queue)
        if cost!=distance.get(node):continue
        if node==last:
            path=[node]
            while node!=first:node=previous[node];path.append(node)
            return np.r_[np.asarray(start)[None],coordinates[path[::-1]],np.asarray(stop)[None]],dict(visited=visited,endpoint_gaps_m=gaps.tolist(),nodes=len(path))
        visited+=1
        if visited%256==0 and time.monotonic()-began>46:
            return None,dict(reason='46-second internal proposal budget reached',visited=visited)
        cell=np.array(np.unravel_index(node,shape))
        for delta,length in neighbors:
            other=cell+delta
            if np.any(other<0)or np.any(other>=dims):continue
            target=int(node+delta@stride)
            if not safe[target]:continue
            # Prevent diagonal edges from cutting an occupied voxel corner.
            if any(not safe[int(node+np.array([delta[k]if j==k else 0 for j in range(3)])@stride)]for k in range(3)if delta[k]):continue
            trial=cost+length
            if trial<distance.get(target,np.inf):
                distance[target]=trial;previous[target]=node
                heapq.heappush(queue,(trial+float(np.linalg.norm(coordinates[target]-coordinates[last])),trial,target))
    return None,dict(reason='Fixed native occupied domain disconnected',visited=visited,endpoint_gaps_m=gaps.tolist())


def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir()
    old_file=BASE/'male-left-facial-cervical-routes-anatomical-v5/routes.json'
    plan=deepcopy(json.loads(old_file.read_bytes()));source=ROOT/plan['source_catalogue'];expected=plan['source_catalogue_sha256']
    assert digest_file(source)==expected
    evidence=json.loads((BASE/'male-left-facial-cervical-finite-network-anatomical-v5/feasibility-v2/evidence.json').read_bytes())
    native_contract=json.loads((BASE/'male-left-facial-cervical-routes-anatomical-v5/construction-contract.json').read_bytes())
    reference=BASE/'male-facial-family-inputs-v1/reference-coarse-facial-corridors-v6.json'
    reference_evidence=json.loads(reference.read_bytes());assert reference_evidence['source_catalogue_sha256']==expected
    spec=importlib.util.spec_from_file_location('facial_volume',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);centers=np.asarray(plan['fixed_parent']['source_station_components'][0]['centers_m']);rings=np.asarray(plan['fixed_parent']['source_station_components'][0]['rings'])
    parent_closed,_=helper.occupancy_classifier(parent,True);assert parent_closed is not None
    fixed_ids=sorted(set(r['part_id']for r in evidence['forbidden_crossings'])|{'ZA:Facial vein.l','ZA:Stylopharyngeus muscle.l','ZA:Submandibular duct.l','ZA:Oropharynx','ZA:Nasopharynx'})
    fixed=[]
    for pid in fixed_ids:
        m=reader.read_part(atlas,parts[pid],frame=parent.frame,catalogue_path=source,catalogue_sha=expected)
        classifier,mode=helper.occupancy_classifier(m,parts[pid].get('system')in{'arterial','venous'})
        fixed.append((pid,m,classifier,mode))
    interval=np.array(native_contract['sg_sp_interval_witness']['authored_interval_point_m'])
    top0,top1=[np.array(w['authored_clearance_waypoint_m'])for w in native_contract['native_upper_margin_witnesses']]
    target={r['name']:np.array(r['recipient']['point_m'])for r in plan['routes']}
    authored={
        'ascending-palatine':dict(station=8,controls=[centers[8],[.034,1.635,.004],interval,top0+[.002,.001,0],top1+[0,.001,0],[.004,1.674,.006],target['ascending-palatine']]),
        'tonsillar':dict(station=16,controls=[centers[16],[.026,1.631,.014],[.022,1.645,.014],[.022,1.655,.022],target['tonsillar']]),
        'glandular-0':dict(station=32,preserve_passing_guide=True),
        'glandular-1':dict(station=48,controls=[centers[48],[.034,1.614,.050],[.029,1.614,.051],target['glandular-1']])}
    contract=deepcopy(native_contract);contract.update(schema='native-source-relative-facial-construction/v7',predecessor_contract_sha256=digest_file(ROOT/plan['construction_contract']),
        reference_relationship_evidence=str(reference.relative_to(ROOT)),reference_relationship_evidence_sha256=digest_file(reference),
        rejected_reference_registration_used_for_coordinates=False,private_geometry_copied=False,
        parent_departure_policy='Native stations 8,16,32,48 preserve source parent controls; changed branch departures use native anatomy, unchanged source-calibre transported sections and real old recipients.',
        native_authored_controls={k:dict(parent_station=v['station'],controls_m=[np.asarray(x).tolist()for x in v.get('controls',[])],passing_old_guide_retained=v.get('preserve_passing_guide',False))for k,v in authored.items()},
        proposal_method='A 0.7-mm three-dimensional native lattice within an 8-mm anatomically authored route domain; exact nearest triangles and closed-source winding propose clearance. Actual emitted surfaces and full occupied volumes remain mandatory.',
        source_calibres_not_reduced=True,limited_tonsillar_constrictor_entry_radius_m=.010,installed=False)
    atomic_json(OUT/'construction-contract.json',contract);atomic_json(ROOT/'reports/ta98-male/facial-native-domain-construction-contract-v7.json',contract)
    rows=[];failed=False
    for route in plan['routes']:
        name=route['name'];cfg=authored[name];station=cfg['station'];route['parent_station']=station
        if cfg.get('preserve_passing_guide'):
            route['guide_retained_from_passing_native_context']=True
            rows.append(dict(name=name,passing_predecessor_guide_literal=True));continue
        baseline,_=centripetal_curve(cfg['controls'],chord_tolerance_m=.00003,max_step_m=.0007)
        step=.0007;margin=.008;low=baseline.min(0)-margin;high=baseline.max(0)+margin
        axes=[np.arange(low[k],high[k]+step/2,step)for k in range(3)];shape=tuple(map(len,axes));q=np.stack(np.meshgrid(*axes,indexing='ij'),-1).reshape(-1,3)
        domain=cKDTree(baseline).query(q)[0]<=margin
        points=q[domain];ids=np.flatnonzero(domain);safe=np.ones(len(points),bool)
        source_delta=parent.vertices[rings[station]]-centers[station]
        factor=np.linalg.norm(source_delta,axis=1).max()/np.linalg.norm(source_delta,axis=1).mean()
        radius=route['illustrative_radius_m']*factor+.00010
        pd=np.sqrt(igl.point_mesh_squared_distance(points,parent.vertices,parent.faces)[0]);pi=abs(igl.winding_number(parent_closed.vertices,parent_closed.faces,points))>.99
        source_inside=pi&(pd>radius)
        blockers=[]
        for pid,mesh,classifier,mode in fixed:
            if time.monotonic()-began>42:
                failed=True;blockers.append(dict(part_id=pid,reason='42-second internal proposal-query budget reached'));break
            selected=np.flatnonzero(safe&np.all(points>=mesh.bounds[0]-radius,axis=1)&np.all(points<=mesh.bounds[1]+radius,axis=1))
            if not len(selected):continue
            p=points[selected];distance=np.sqrt(igl.point_mesh_squared_distance(p,mesh.vertices,mesh.faces)[0]);bad=distance<=radius
            if classifier is not None:
                within=np.all(p>=classifier.bounds[0],axis=1)&np.all(p<=classifier.bounds[1],axis=1)
                bad[within]|=abs(igl.winding_number(classifier.vertices,classifier.faces,p[within]))>.99
            if name=='tonsillar'and pid=='ZA:Superior pharyngeal constrictor.l':bad&=np.linalg.norm(p-target[name],axis=1)>.010
            bad&=~source_inside[selected]
            safe[selected[bad]]=False;blockers.append(dict(part_id=pid,geometry_sha256=mesh.geometry_sha,queried=len(selected),rejected=int(bad.sum()),closed_occupancy_mode=mode))
        mask=np.zeros(len(q),bool);mask[ids]=safe
        atomic_npz(OUT/(name+'-occupied-domain.npz'),axes_x=axes[0],axes_y=axes[1],axes_z=axes[2],domain_mask=domain,source_calibre_safe_mask=mask,frame=np.array(parent.frame))
        report=dict(name=name,parent_station=station,source_profile_enclosing_radius_m=radius-.0001,clearance_guard_m=.0001,grid_spacing_m=step,search_domain_m=margin,proposal_points=len(points),safe_points=int(safe.sum()),fixed_boundaries=blockers)
        if failed:
            rows.append(report);break
        # Palatine must use the actual SG/SP interval and pass superior to the
        # two actual constrictor source witnesses; other paths use medial MP.
        anchors=[centers[station],interval,top0+[.002,.001,0],top1+[0,.001,0],target[name]]if name=='ascending-palatine'else[centers[station],target[name]]
        paths=[];segments=[]
        for start,stop in zip(anchors[:-1],anchors[1:]):
            curve,gate=lattice_path(mask,start,stop,shape,q,step,began);segments.append(gate)
            if curve is None:failed=True;break
            paths.extend(curve if not paths else curve[1:])
        report['anatomical_anchor_segments']=segments
        if failed:rows.append(report);break
        path=np.asarray(paths);path=path[np.r_[True,np.linalg.norm(np.diff(path,axis=0),axis=1)>1e-10]]
        smoothed=gaussian_filter1d(path,1.3,axis=0,mode='nearest');delta=smoothed-path;length=np.linalg.norm(delta,axis=1);delta*=np.minimum(1.,.00015/np.maximum(length,1e-15))[:,None];path+=delta;path[0]=centers[station];path[-1]=target[name]
        path,_=centripetal_curve(path,chord_tolerance_m=.00001,max_step_m=.0004)
        file=OUT/(name+'-guide.npz');atomic_npz(file,points=path,outer_radii=np.full(len(path),route['illustrative_radius_m']),inner_radii=np.full(len(path),route['illustrative_radius_m']-.00009),frame=np.array(parent.frame))
        route.update(guide=str(file.relative_to(ROOT)),guide_sha256=digest_file(file),stations=len(path),length_m=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum()))
        report.update(proposed_route_length_m=route['length_m'],proposed_stations=len(path));rows.append(report)
    result=dict(schema='native-reference-guided-facial-corridor/v7',source_catalogue_sha256=expected,stage=1,maximum_stages=3,stopped_before_finite_construction=failed,
                all_source_calibre_corridors_found=not failed,routes=rows,precise_donor_registration_used=False,installed=False,elapsed_seconds=time.monotonic()-began)
    if not failed:
        plan.update(construction_contract=str((OUT/'construction-contract.json').relative_to(ROOT)),construction_contract_sha256=digest_file(OUT/'construction-contract.json'),predecessor_route_contract=str(old_file.relative_to(ROOT)),predecessor_route_contract_sha256=digest_file(old_file))
        atomic_json(OUT/'routes.json',plan)
    assert digest_file(source)==expected
    atomic_json(OUT/'stage1-evidence.json',result)
    print(json.dumps(result))


if __name__=='__main__':main()
