"""Finite read-only renal source-mouth/host audit; no construction or joins."""
import argparse,itertools,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit,topology
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_continuity import boundary_loops,continuity_audit
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39'
LEFT=['ZA:Left renal artery','ZA:Anterior branch of renal artery.l','ZA:Posterior branch of renal artery.l','ZA:Left renal vein']
NETWORKS=['ZA:Intrarenal arteries of left kidney','ZA:Intrarenal veins of left kidney']
HOSTS=['ZA:Abdominal aorta','ZA:Inferior vena cava (abdominal part)','RENALDETAIL:male:left:capsule','ZA:Kidney.l']
def read(p):return json.loads(Path(p).read_text())
def require(x,msg):
    if not x:raise ValueError(msg)
def ring_area(p):return float(np.linalg.norm(np.cross(p-p.mean(0),np.roll(p-p.mean(0),-1,axis=0)).sum(0))/2)
def ring_match(a,b):
    if len(a)!=len(b):return dict(comparable=False,first_controls=len(a),second_controls=len(b),literal_shared_mouth=False)
    choices=[]
    for reverse in [False,True]:
        for phase in range(len(b)):
            order=np.roll(np.arange(len(b))[::-1]if reverse else np.arange(len(b)),phase);gaps=np.linalg.norm(a-b[order],axis=1);choices.append((float(gaps.max()),reverse,phase,order,gaps))
    gap,reverse,phase,order,gaps=min(choices,key=lambda r:r[:3]);aa,ab=ring_area(a),ring_area(b)
    return dict(comparable=True,maximum_control_gap_m=gap,minimum_control_gap_m=float(gaps.min()),rms_control_gap_m=float(np.sqrt(np.mean(gaps*gaps))),reverse=reverse,phase=phase,second_control_order=order.tolist(),first_area_m2=aa,second_area_m2=ab,relative_area_difference_first=abs(ab-aa)/aa,exact_common_controls=len(set(map(tuple,a))&set(map(tuple,b))),literal_shared_mouth=bool(np.array_equal(a,b[order])),join_inferred_from_gap=False)
def main(out):
    began=time.monotonic();require(out.is_relative_to(WORK)and not out.exists(),'Use a new local renal audit directory');require(digest_file(CAT)==EXPECTED,'Current male revision changed; explicitly rebind audit');out.mkdir(parents=True);atlas=read(CAT);records={r['id']:r for r in atlas['parts']};reader=AtlasReader(ROOT);meshes={};parts=[];mouths={};bindings={str(CAT):EXPECTED,str(Path(__file__).resolve()):digest_file(__file__)}
    for i,pid in enumerate(LEFT+NETWORKS+HOSTS):
        require(time.monotonic()-began<45,'Source audit approaching deadline');m=reader.read_part(atlas,records[pid],frame='model:local-ta98:assembled-metres',catalogue_path=CAT,catalogue_sha=EXPECTED);meshes[pid]=m;bindings[m.source['bank_path']]=m.source['bank_sha256'];path=out/('source-'+str(i).zfill(2)+'.npz');m.save(path);bindings[str(path)]=digest_file(path);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];topo=topology(m);loops=[];err=None
        try:loops=boundary_loops(f)
        except ValueError as error:err=str(error)
        mouths[pid]=[v[r]for r in loops];arrays=dict(physical_vertices=v,physical_faces=f,original_chart_to_physical_vertex=inv);arrays.update({'mouth_'+str(j):r for j,r in enumerate(loops)});serial=[];serial_failure=None
        if pid in LEFT+NETWORKS:
            try:
                for k,c in enumerate(recover_station_rings(v,f)):
                    aligned,proof=align_station_correspondence(v,f[c['face_ids']],c['rings']);arrays['component_'+str(k)+'_source_faces']=c['face_ids'];arrays['component_'+str(k)+'_rings']=np.array(aligned);serial.append(dict(component=k,stations=len(aligned),ring_controls=len(aligned[0]),all_original_faces_replayed=all(r['sourceFacesReplayedExactly']for r in proof),source_interval_proofs=proof))
            except ValueError as error:serial_failure=str(error)
        witness=out/('rings-'+str(i).zfill(2)+'.npz');atomic_npz(witness,**arrays);row=dict(part_id=pid,name=records[pid]['name'],source_record=records[pid],geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,mesh=str(path),mesh_sha256=digest_file(path),bank=m.source['bank_path'],bank_sha256=m.source['bank_sha256'],topology=topo,bounds_m=m.bounds.tolist(),mouths=[dict(index=j,controls=len(r),center_m=v[r].mean(0).tolist(),area_m2=ring_area(v[r]),bounds_m=[v[r].min(0).tolist(),v[r].max(0).tolist()])for j,r in enumerate(loops)],boundary_failure=err,ordered_source_cylinders=serial,serial_recovery_failure=serial_failure,ring_witnesses=str(witness),ring_witnesses_sha256=digest_file(witness),finite_wall_established=bool(records[pid].get('pairedTissue')and records[pid].get('sectionSolid')))
        if pid in LEFT:row['actual_open_sheet_self_audit']=audit(m,require_closed=False,max_pairs=250000,max_tests=2000000)
        parts.append(row);atomic_json(out/('part-'+str(i).zfill(2)+'.json'),row)
    relationp=ROOT/'reports/structural-dependencies/current/male-anatomical-edges.json';relations=read(relationp);require(relations['sourceSha256']==EXPECTED,'Declared relation snapshot is stale');bindings[str(relationp)]=digest_file(relationp);declared=[r for r in relations['edges']if r['from']in LEFT+NETWORKS or r['to']in LEFT+NETWORKS];pairs=[]
    for a,b in itertools.combinations(LEFT[:3],2):
        variants=[dict(first_mouth=i,second_mouth=j,**ring_match(ra,rb))for i,ra in enumerate(mouths[a])for j,rb in enumerate(mouths[b])];best=min(variants,key=lambda r:r.get('maximum_control_gap_m',float('inf')));arrays,contact=pair_contacts(meshes[a],meshes[b],0.,max_pairs=250000,max_tests=2000000);wp=out/('fork-'+str(len(pairs))+'.npz');atomic_npz(wp,**arrays);pairs.append(dict(first=a,second=b,all_mouth_comparisons=variants,closest_explicit_mouth_pair=best,actual_surface_contact=contact,contact_witnesses=str(wp),contact_witnesses_sha256=digest_file(wp),anatomical_join_established=False))
    host_checks=[]
    for pid,host in [('ZA:Left renal artery','ZA:Abdominal aorta'),('ZA:Left renal vein','ZA:Inferior vena cava (abdominal part)')]+[(p,'RENALDETAIL:male:left:capsule')for p in LEFT]:
        require(time.monotonic()-began<45,'Host audit approaching deadline');m,h=meshes[pid],meshes[host];tree=BVH(h);all_queries=[]
        for i,ring in enumerate(mouths[pid]):
            points=[]
            for j,q in enumerate(ring):
                r=tree.closest_point(q);points.append(dict(control=j,source_point=q.tolist(),host_point=r['point'].tolist(),host_face=int(r['face_id']),host_barycentric=r['barycentric'].tolist(),gap_m=float(r['distance_m'])))
            all_queries.append(dict(mouth=i,minimum_gap_m=min(q['gap_m']for q in points),maximum_gap_m=max(q['gap_m']for q in points),median_gap_m=float(np.median([q['gap_m']for q in points])),all_controls=points))
        arr,contact=pair_contacts(m,h,0.,second_bvh=tree,max_pairs=250000,max_tests=2000000);wp=out/('host-'+str(len(host_checks))+'.npz');atomic_npz(wp,**arr);host_checks.append(dict(first=pid,second=host,first_geometry_sha256=m.geometry_sha,second_geometry_sha256=h.geometry_sha,mouth_to_actual_host_surface=all_queries,existing_host_open_mouth_comparisons=[dict(source_mouth=i,host_mouth=j,**ring_match(a,b))for i,a in enumerate(mouths[pid])for j,b in enumerate(mouths[host])],surface_contact=contact,witnesses=str(wp),witnesses_sha256=digest_file(wp),crossings_do_not_certify_anatomical_entry=True,anatomical_interface_accepted=False))
    archive=[];v350p=ROOT/'reports/renal-vascular-master-replay-v350.json';v350=read(v350p);bindings[str(v350p)]=digest_file(v350p)
    for rp,pid in [('reports/renal-arterial-refined-course-v337.json',NETWORKS[0]),('reports/renal-transported-vein-clearance-v349.json',NETWORKS[1])]:
        rp=ROOT/rp;r=read(rp);npz=ROOT/r['candidate'];require(digest_file(npz)==r['candidateSha256'],'Historical renal archive changed');z=np.load(npz);old=Mesh(z['candidate_v'],z['candidate_f'],meshes[pid].frame);bindings[str(rp)]=digest_file(rp);bindings[str(npz)]=digest_file(npz);archive.append(dict(report=str(rp),report_sha256=digest_file(rp),archive=str(npz),archive_sha256=digest_file(npz),current_part_id=pid,current_geometry_sha256=meshes[pid].geometry_sha,archive_geometry_sha256=old.geometry_sha,actual_triangles_literal=np.array_equal(meshes[pid].triangles,old.triangles),current_faces=len(meshes[pid].faces),archive_faces=len(old.faces),old_gate_reusable=old.geometry_sha==meshes[pid].geometry_sha))
    # Anatomical direction is explicit target scope, distinct from cached UI
    # before/after order. No geometric edge is created from a near mouth.
    flow=[dict(first='ZA:Abdominal aorta',second=LEFT[0],role='arterial parent'),dict(first=LEFT[0],second=LEFT[1],role='anterior arterial branch'),dict(first=LEFT[0],second=LEFT[2],role='posterior arterial branch'),dict(first=NETWORKS[1],second=LEFT[3],role='renal venous tributary'),dict(first=LEFT[3],second='ZA:Inferior vena cava (abdominal part)',role='renal venous drainage')]
    current_parts={p['part_id']:p for p in parts};port_centers=np.array([current_parts[p]['mouths'][-1]['center_m']for p in LEFT]);report=dict(schema='male-left-renal-current-source-graph-pilot/v1',catalogue=str(CAT),catalogue_sha256=EXPECTED,source_bindings=bindings,parts=parts,explicit_anatomical_flow_edges=flow,current_cached_relationships=declared,fork_pairs=pairs,host_entries=host_checks,hilar_reference=dict(explicit_hilum_mesh_or_region_present=False,renal_capsule_bounds_m=meshes[HOSTS[2]].bounds.tolist(),native_whole_kidney_bounds_m=meshes[HOSTS[3]].bounds.tolist(),source_terminal_mouth_center_envelope_m=[port_centers.min(0).tolist(),port_centers.max(0).tolist()],terminal_mouth_envelope_is_not_hilum_segmentation=True),historical_v350=dict(path=str(v350p),sha256=digest_file(v350p),ordinary_male_catalogue_sha256=v350['ordinaryCatalogueSha256']['male'],current_catalogue_matches=False,reported_arterial_components=v350['arterialFlowComponents'],reported_missing_arterial_connections=v350['missingArterialConnections'],current_physical_comparison=archive,geometry_pass_reused=False),next_construction_scope=dict(first=LEFT[0],second=LEFT[1],third=LEFT[2],kind='One explicit trunk/anterior/posterior renal Y-fork; shared finite wall and open blood-space topology must be constructed together, not a triple coincident ring weld',authorized=False),skeleton_modified=False,source_geometry_modified=False,openings_capped=False,proximity_welding=False,actual_wall_lumen_patency_proven=False,qualifying_candidate=False,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    require(all(digest_file(p)==s for p,s in bindings.items()),'Bound renal source changed during audit');atomic_json(out/'evidence.json',report);print(json.dumps(dict(output=str(out),parts=len(parts),fork_gaps_m=[p['closest_explicit_mouth_pair']['maximum_control_gap_m']for p in pairs],literal_shared_fork_mouths=sum(p['closest_explicit_mouth_pair']['literal_shared_mouth']for p in pairs),historical_geometry_reusable=[p['old_gate_reusable']for p in archive],elapsed_seconds=report['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
