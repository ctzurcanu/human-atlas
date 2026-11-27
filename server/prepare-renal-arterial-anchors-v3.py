"""Actual local corticomedullary anchors and proximal ports; no new tissue meshes."""
import json,sys,time
from pathlib import Path
from collections import defaultdict
import numpy as np
import point_cloud_utils as pcu
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components,dijkstra
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.surface_rays import first_surface_hits
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.contact import pair_contacts
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/renal-blood-frontier-v1'
DEST=BASE/'arterial-course-v2/a-local-tree-v2'

def sd(mesh,p):
    q=np.ascontiguousarray(p);single=len(q)==1
    # This installed PCU wrapper returns a numerically wrong single-point SDF;
    # the retained single-versus-batch witness proves the discrepancy. Keep the
    # identical query in the correctly working batch path, without a tolerance.
    if single:q=np.repeat(q,2,axis=0)
    result=np.atleast_1d(pcu.signed_distance_to_mesh(q,np.ascontiguousarray(mesh.vertices),np.ascontiguousarray(mesh.faces,dtype=np.int32))[0])
    return result[:1]if single else result

def main():
    start=time.monotonic();assert not DEST.exists();DEST.mkdir();ev=json.loads((BASE/'evidence.json').read_text());assert digest_file(ev['catalogue'])==ev['catalogue_sha256'];bind={r['part_id']:r for r in ev['source_bindings']}
    def load(pid):
        b=bind[pid];assert digest_file(b['mesh'])==b['mesh_sha256'];return Mesh.load(b['mesh'])
    cortex=load('RENALDETAIL:male:left:cortex');pyr=load('RENALDETAIL:male:left:pyramid:a');s=np.load(BASE/'shared-cortex-pyramid-a.npz');tri=s['pyramid_triangles'];centers=tri.mean(1);normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);areas=np.linalg.norm(normals,axis=1)/2;normals/=2*areas[:,None]
    verts,iv=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);faces=iv.reshape(-1,3);edge_faces=defaultdict(list)
    for i,f in enumerate(faces):
        for a,b in zip(f,np.roll(f,-1)):edge_faces[tuple(sorted([int(a),int(b)]))].append(i)
    links=[(a,b)for ids in edge_faces.values()for a in ids for b in ids if a!=b];ee=np.array(links);weights=np.linalg.norm(centers[ee[:,0]]-centers[ee[:,1]],axis=1)
    graph=coo_matrix((weights,ee.T),shape=(len(tri),len(tri))).tocsr();nc,lab=connected_components(graph,directed=False);cc=max(range(nc),key=lambda i:areas[lab==i].sum());ids=np.flatnonzero(lab==cc)
    rays,ray_summary=first_surface_hits(cortex,centers[ids],normals[ids],minimum_m=1e-7,maximum_m=.025,time_limit=5,max_triangle_tests=200000)
    valid=rays['hit']&(rays['distance_m']>.0008)
    # Validate the measured first-hit interval enters the actual finite cortex.
    midpoint=centers[ids]+normals[ids]*rays['distance_m'][:,None]*.5;midpoint[~rays['hit']]=centers[ids][~rays['hit']];inside=sd(cortex,midpoint)<0;valid&=inside
    retained=ids[valid];assert len(retained)>2,'No connected finite cortical patch with a useful normal interval'
    restricted=graph[retained][:,retained];n,l=connected_components(restricted,directed=False);largest=max(range(n),key=lambda i:areas[retained[l==i]].sum());sel=np.flatnonzero(l==largest);actual=retained[sel];gg=restricted[sel][:,sel]
    # Two exact graph sweeps define a local surface arc on this connected patch,
    # rather than an angular bin or a global-coordinate plane.
    dd=dijkstra(gg,directed=False,indices=0);a=int(np.argmax(dd));dd,prev=dijkstra(gg,directed=False,indices=a,return_predecessors=True);b=int(np.argmax(dd));chain=[b]
    while chain[-1]!=a:
        nxt=int(prev[chain[-1]]);assert nxt>=0;chain.append(nxt)
    chain=chain[::-1];patch_ids=actual[chain]
    rayindex={int(k):j for j,k in enumerate(ids)}
    def anchor(k):
        j=rayindex[int(k)];length=float(rays['distance_m'][j]);offset=min(.0006,length/2);base=centers[k];normal=normals[k];arcuate=base+normal*offset
        return dict(shared_patch_face_index=int(k),pyramid_source_face_id=int(s['pyramid_face_ids'][k]),cortex_source_face_id=int(s['cortex_face_ids'][k]),barycentric=[1/3]*3,
                    actual_contact_m=base.tolist(),pyramid_outward_to_cortex_normal=normal.tolist(),arcuate_center_candidate_m=arcuate.tolist(),arcuate_offset_m=offset,
                    first_cortex_boundary_m=rays['points'][j].tolist(),first_cortex_boundary_face=int(rays['face_ids'][j]),first_cortex_boundary_barycentric=rays['barycentric'][j].tolist(),
                    finite_cortex_normal_interval_m=length,arcuate_center_inside_cortex=bool(sd(cortex,arcuate[None])[0]<0),whole_swept_clearance_proven=False)
    arc=[anchor(int(k))for k in patch_ids];arcpoints=np.array([x['arcuate_center_candidate_m']for x in arc]);arclength=np.r_[0,np.cumsum(np.linalg.norm(np.diff(arcpoints,axis=0),axis=1))];j=int(np.argmin(abs(arclength-arclength[-1]/2)));takeoff=arc[j]
    base=np.array(takeoff['actual_contact_m']);normal=np.array(takeoff['pyramid_outward_to_cortex_normal']);end=base+normal*takeoff['finite_cortex_normal_interval_m']*.8
    cortical=dict(takeoff_arc_index=j,source_contact=takeoff,centerline_start_m=takeoff['arcuate_center_candidate_m'],centerline_end_m=end.tolist(),end_inside_actual_cortex=bool(sd(cortex,end[None])[0]<0),
                  terminal_is_open_unresolved_cortical_recipient=True,no_glomerulus_or_microbed_claim=True)
    # Measure actual current fork support and retain literal upstream profiles
    # if a distal end is unsuitable; never freeze an invalid end by convention.
    rp=ROOT/'reports/renal-ureter-fluid-v317.json';fr=json.loads(rp.read_text());archive=ROOT/fr['candidate'];assert digest_file(archive)==fr['candidateSha256'];fz=np.load(archive);fluid=Mesh(fz['urinePositions'],fz['urineFaces'],pyr.frame)
    wallpids=['RENALDETAIL:male:left:collecting-wall:'+c for c in ['adventitia','muscularis','mucosa']];walls={pid:load(pid)for pid in wallpids};forkrows=[];contacts=[]
    for i,pid in enumerate(['ZA:Left renal artery','ZA:Anterior branch of renal artery.l','ZA:Posterior branch of renal artery.l']):
        assert time.monotonic()-start<45,'Anchor stage nearing deadline';m=load(pid);port=np.load(BASE/f'fork-port-{i}.npz');pp=np.r_[port['outer'],port['inner']];ud=sd(fluid,pp)
        row=dict(part_id=pid,legacy_fluid_port_signed_min_m=float(ud.min()),legacy_fluid_port_inside_controls=int(np.sum(ud<0)),legacy_fluid_volume_is_unreplayed_reference=True,current_wall_port_containment={})
        for wp,wm in walls.items():
            wd=sd(wm,pp);row['current_wall_port_containment'][wp]=dict(minimum_m=float(wd.min()),inside_controls=int(np.sum(wd<0)))
            rr,ss=pair_contacts(m,wm,0.,max_pairs=20000,max_tests=400000);assert not ss['truncated'];p=DEST/f'fork-{i}-collecting-{wp.rsplit(":",1)[1]}.npz';atomic_npz(p,**rr);contacts.append(dict(fork=pid,wall=wp,summary=ss,witnesses=str(p),sha256=digest_file(p)))
        if i:
            outer=m.source['record']['tissueDrawGroups'][0];sfids=np.arange(outer['start']//3,(outer['start']+outer['count'])//3);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);ff=inv[m.faces];sequences=recover_station_rings(v,ff[sfids]);assert len(sequences)==1
            seq=sequences[0];rings,proof=align_station_correspondence(v,ff[sfids][seq['face_ids']],seq['rings']);pr=v[np.array(rings)]
            neighbor=next(g for g in m.source['record']['tissueDrawGroups']if g.get('neighbor')=='ZA:Left renal artery');junction=m.vertices[m.faces[neighbor['start']//3:(neighbor['start']+neighbor['count'])//3]].reshape(-1,3).mean(0)
            if np.linalg.norm(pr[-1].mean(0)-junction)<np.linalg.norm(pr[0].mean(0)-junction):pr=pr[::-1]
            rc=pr.mean(1);legacy=sd(fluid,pr.reshape(-1,3)).reshape(len(pr),12);block=np.zeros((len(pr),12),bool)
            for wm in walls.values():block|=sd(wm,pr.reshape(-1,3)).reshape(len(pr),12)<0
            bad=np.any(block|(legacy<0),axis=1);firstbad=int(np.flatnonzero(bad)[0])if bad.any()else None
            p=DEST/f'fork-{i}-ordered-outer-profiles.npz';atomic_npz(p,profiles=pr,centers=rc,legacy_fluid_distance_m=legacy,current_collecting_tissue_inside=block,original_outer_source_faces=sfids)
            row.update(ordered_outer_profiles=str(p),ordered_outer_profiles_sha256=digest_file(p),stations=len(pr),first_sampled_blocked_station=firstbad,
                last_sampled_clear_station_before_first_block=(firstbad-1 if firstbad is not None and firstbad>0 else None),source_profile_clearance_not_continuous_gate=True)
        forkrows.append(row)
    patchout=DEST/'local-corticomedullary-arc.npz';atomic_npz(patchout,actual_shared_patch_face_indices=patch_ids,pyramid_actual_source_face_ids=s['pyramid_face_ids'][patch_ids],cortex_actual_source_face_ids=s['cortex_face_ids'][patch_ids],contact_points=centers[patch_ids],outward_normals=normals[patch_ids],arcuate_centerline_candidates=arcpoints,cortical_radiate_end=end)
    report=dict(schema='renal-A-local-arterial-anchor-contract/v3',catalogue_sha256=ev['catalogue_sha256'],source_evidence_sha256=digest_file(BASE/'evidence.json'),
        anatomical_route='A new interlobar trunk in the tissue column beside pyramid A feeds one local arcuate path outside the actual A cortex-medulla patch; a cortical-radiate daughter runs outward inside the measured finite cortical interval. No arterial loop or segmental anastomosis is inferred.',
        patch=dict(total_shared_faces=len(tri),connected_components=nc,selected_component=cc,selected_component_faces=len(ids),selected_component_area_m2=float(areas[ids].sum()),retained_connected_faces=len(actual),surface_arc_length_m=float(arclength[-1]),ray_summary=ray_summary),
        arc_anchors=arc,cortical_radiate=cortical,profile_source=ev['external_fork_ports'],current_fork=forkrows,actual_fork_collecting_contacts=contacts,
        route_attachment_not_inferred_from_nearest_patch=True,column_entry_corridor_requires_root_volume_grid=True,
        exact_patch_archive=str(patchout),exact_patch_sha256=digest_file(patchout),new_geometry_constructed=False,installed=False,elapsed_seconds=time.monotonic()-start)
    atomic_json(DEST/'anchors.json',report);print(json.dumps(dict(output=str(DEST/'anchors.json'),seconds=report['elapsed_seconds'],patch=report['patch'],arc_endpoints=[arc[0],arc[-1]],cortical_radiate=cortical,current_fork=forkrows)))
if __name__=='__main__':main()
