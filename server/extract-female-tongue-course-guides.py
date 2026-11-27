"""Actual native female extrinsic section guides; no source mesh replacement.

Every loop is a literal source/plane intersection with source face and
barycentric witnesses. Interior section centroids are fitting references,
not new anatomical landmarks, physical caps or copied replacement meshes.
"""
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from tongue_whole_input import WORK

BASE=WORK/'female-native-tongue-baseline-v1/baseline.json'

def bound(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=digest_file(p))

def endpoint_witness(m,vertex):
    candidates=np.argwhere(m.faces==vertex);face,corner=map(int,candidates[0]);bary=np.eye(3)[corner]
    return dict(source_face=face,barycentric=bary.tolist(),point_m=m.vertices[vertex].tolist(),source_vertex=int(vertex))

def section(m,axis,level,basis):
    # Exact-coordinate welding follows the source topology audit. Two source
    # UV vertices at the same coordinate do not create a fabricated boundary.
    vertices,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];dist=vertices@axis-level
    if np.min(np.abs(dist))<1e-12:
        level+=3.7e-10;dist=vertices@axis-level
    events={};segments=[]
    for face,tri in enumerate(f):
        d=dist[tri];ends=[]
        for k in range(3):
            j=(k+1)%3
            if d[k]*d[j]>=0:continue
            edge=tuple(sorted([int(tri[k]),int(tri[j])]))
            if edge not in events:
                aa,bb=edge;t=float(dist[aa]/(dist[aa]-dist[bb]));point=vertices[aa]*(1-t)+vertices[bb]*t
                # Reconstruct this event from this actual source face/corner.
                ak=int(np.flatnonzero(tri==aa)[0]);bk=int(np.flatnonzero(tri==bb)[0]);bary=np.zeros(3);bary[ak]=1-t;bary[bk]=t
                events[edge]=dict(id=len(events),point=point,source_face=int(face),barycentric=bary)
            ends.append(events[edge]['id'])
        if len(ends)==2:segments.append(ends)
        elif ends:raise ValueError('Section intersects an unresolved source vertex')
    if not segments:return dict(plane_offset_m=float(level),loops=[],failure='no section segments')
    rows=sorted(events.values(),key=lambda r:r['id']);adj={i:[]for i in range(len(rows))}
    for aa,bb in segments:adj[aa].append(bb);adj[bb].append(aa)
    if any(len(x)!=2 for x in adj.values()):return dict(plane_offset_m=float(level),loops=[],failure='not all source segment endpoints have degree two')
    pending=set(adj);loops=[]
    while pending:
        loop=[min(pending)];last=None
        while True:
            now=loop[-1];nxt=next(q for q in adj[now]if q!=last)
            if nxt==loop[0]:break
            if nxt in loop:raise ValueError('Source loop repeats before closing')
            loop.append(nxt);last=now
        pending.difference_update(loop);p=np.array([rows[i]['point']for i in loop]);xy=p@basis.T
        cross=xy[:,0]*np.roll(xy[:,1],-1)-np.roll(xy[:,0],-1)*xy[:,1];area=cross.sum()/2
        if abs(area)<1e-16:return dict(plane_offset_m=float(level),loops=[],failure='vanishing section area')
        centroid=np.sum((xy+np.roll(xy,-1,axis=0))*cross[:,None],axis=0)/(6*area)
        center=centroid@basis+level*axis
        loops.append(dict(positions=p,source_face=np.array([rows[i]['source_face']for i in loop]),
            barycentric=np.array([rows[i]['barycentric']for i in loop]),centroid_m=center,
            area_m2=abs(float(area)),source_boundary_points=len(loop)))
    return dict(plane_offset_m=float(level),loops=loops,failure=None)

def main(output,only=None):
    began=time.monotonic()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new bounded guide stage')
    baseline=json.loads(BASE.read_text());sources={r['part_id']:r for r in baseline['sources']};output.mkdir()
    bones={pid:Mesh.load(r['mesh']['path'])for pid,r in sources.items()if r['role']=='bone'};bone_trees={pid:BVH(m)for pid,m in bones.items()}
    tongue=Mesh.load(sources['LOCAL:female:tongue:1089']['mesh']['path']);tt=BVH(tongue)
    wanted=['LOCAL:female:l_styloglossus:243','LOCAL:female:r_styloglossus:256',
        'LOCAL:female:genioglossus:239','LOCAL:female:l_hyoglossus:242','LOCAL:female:r_hyoglossus:255']
    if only:
        if any(pid not in wanted for pid in only):raise ValueError('Only inventoried native extrinsics are eligible')
        wanted=[pid for pid in wanted if pid in only]
    report=[]
    for number,pid in enumerate(wanted):
        row=sources[pid]
        if digest_file(row['mesh']['path'])!=row['mesh']['sha256']:raise ValueError('Frozen native source changed')
        m=Mesh.load(row['mesh']['path']);center=m.vertices.mean(0);_,_,basis3=np.linalg.svd(m.vertices-center,full_matrices=False);axis=basis3[0]
        if 'styloglossus'in pid:
            proximal=int(np.argmax(m.vertices[:,1]));distal=int(np.argmax(m.vertices[:,2]));interpretation='Superior/posterior styloid-origin side to anterior tongue side; source cortical footprint remains provisional.'
        elif 'hyoglossus'in pid:
            proximal=int(np.argmin(m.vertices[:,1]));distal=int(np.argmax(m.vertices[:,2]));interpretation='Inferior hyoid side toward anterior/lateral tongue; sparse sheet, not one cylindrical fibre bundle.'
        else:
            proximal=int(np.argmax(m.vertices[:,2]));distal=int(np.argmin(m.vertices[:,2]));interpretation='Geometric anterior-to-posterior section sequence of one bilateral fan; not a unique anatomical fascicle or two separately owned GG muscles.'
        if (m.vertices[distal]-m.vertices[proximal])@axis<0:axis=-axis
        b0=basis3[1];b0-=axis*(b0@axis);b0/=np.linalg.norm(b0);b1=np.cross(axis,b0);basis=np.array([b0,b1])
        projection=m.vertices@axis;low,high=projection.min(),projection.max();start=int(np.argmin(projection));end=int(np.argmax(projection));margin=min(.0005,(high-low)*.025)
        stations=np.linspace(low+margin,high-margin,25);sections=[];allp=[];allface=[];allbary=[];loopoffsets=[0];sectionloopoffsets=[0];course=[];course_sections=[];trunk_open=True
        for level in stations:
            s=section(m,axis,level,basis);rows=[]
            for loop in s['loops']:
                rows.append(dict(centroid_m=loop['centroid_m'].tolist(),area_m2=loop['area_m2'],source_boundary_points=loop['source_boundary_points'],point_start=loopoffsets[-1],point_stop=loopoffsets[-1]+len(loop['positions'])))
                allp.extend(loop['positions']);allface.extend(loop['source_face']);allbary.extend(loop['barycentric']);loopoffsets.append(loopoffsets[-1]+len(loop['positions']))
            sectionloopoffsets.append(len(loopoffsets)-1)
            if len(rows)!=1:trunk_open=False
            if len(rows)==1 and trunk_open:
                course.append(np.array(rows[0]['centroid_m']));course_sections.append(len(sections))
            sections.append(dict(plane_offset_m=s['plane_offset_m'],loops=rows,failure=s['failure']))
        curve=np.array(course);tangents=np.zeros_like(curve)
        if len(curve)>1:
            tangents[0]=curve[1]-curve[0];tangents[-1]=curve[-1]-curve[-2];tangents[1:-1]=curve[2:]-curve[:-2]
            tangents/=np.maximum(np.linalg.norm(tangents,axis=1),1e-30)[:,None]
        path=output/f'course-{number:02d}.npz';atomic_npz(path,section_boundary_points=np.array(allp),source_face=np.array(allface),source_barycentric=np.array(allbary),
            loop_point_offsets=np.array(loopoffsets),section_loop_offsets=np.array(sectionloopoffsets),
            course_centroids=curve,course_tangents=tangents,course_section_indices=np.array(course_sections),plane_axis=axis,plane_basis=basis,
            plane_offsets=np.array([s['plane_offset_m']for s in sections]),frame=np.array(m.frame))
        top=endpoint_witness(m,proximal);last=endpoint_witness(m,distal);axisstart=endpoint_witness(m,start);axisend=endpoint_witness(m,end)
        targets=[]
        for site_name,site in [('anatomical-origin-side',top),('axis-start',axisstart),('axis-end',axisend)]:
            q=np.array(site['point_m']);candidates=[]
            for target,b in bone_trees.items():
                h=b.closest_point(q);candidates.append(dict(target_part_id=target,target_source=sources[target]['mesh'],distance_m=float(h['distance_m']),
                    source_face=int(h['face_id']),barycentric=h['barycentric'].tolist(),point_m=h['point'].tolist()))
            h=tt.closest_point(q);targets.append(dict(site=site_name,nearest_bone=min(candidates,key=lambda r:r['distance_m']),
                nearest_tongue=dict(distance_m=float(h['distance_m']),source_face=int(h['face_id']),barycentric=h['barycentric'].tolist(),point_m=h['point'].tolist())))
        residual=max((float(np.linalg.norm(np.array(b)@m.triangles[int(f)]-np.array(q)))for q,f,b in zip(allp,allface,allbary)),default=0)
        record=dict(part_id=pid,source=row['mesh'],source_geometry_sha256=m.geometry_sha,guide=bound(path),
            origin_side_surface_witness=top,tongue_side_surface_witness=last,axis_start=axisstart,axis_end=axisend,
            start_tangent=None if len(tangents)==0 else tangents[0].tolist(),end_tangent=None if len(tangents)==0 else tangents[-1].tolist(),
            section_axis=axis.tolist(),section_basis=basis.tolist(),sections=sections,
            maximum_source_barycentric_residual_m=residual,closed_sections=sum(bool(s['loops'])for s in sections),
            multiple_loop_sections=sum(len(s['loops'])>1 for s in sections),source_target_witnesses=targets,
            course_centroid_scope='Continuous initial sequence of single real loops only; no centroid is invented between disconnected branch loops.',
            course_complete_through_all_sections=len(course)==len(sections),
            interpretation=interpretation,source_section_only=True,no_cap_added=True,no_mesh_copy_or_fit=True)
        atomic_json(output/f'course-{number:02d}.json',record);report.append(dict(part_id=pid,record=bound(output/f'course-{number:02d}.json'),guide=bound(path),
            origin_side_surface_witness=top,start_tangent=record['start_tangent'],end_tangent=record['end_tangent'],
            closed_sections=record['closed_sections'],multiple_loop_sections=record['multiple_loop_sections'],maximum_source_barycentric_residual_m=residual))
        # SG is first so its actual endpoints are available to the fitter
        # immediately, before the broad GG/HG section sequences finish.
        atomic_json(output/f'progress-{number+1:02d}.json',dict(baseline=bound(BASE),guides=report,elapsed_seconds=time.monotonic()-began))
    final=dict(schema='native-female-extrinsic-course-guides/v1',baseline=bound(BASE),guides=report,installed=False,
        source_geometry_modified=False,no_replacement_meshes=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'guides.json',final);print(json.dumps(final),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--only',nargs='+');a=p.parse_args();main(a.output.resolve(),a.only)
