"""Fit an adjoining palm artery and original-detail finger branch as one route.

The two surfaces must traverse one connected recipient compartment. A clear
child on the other side of a fascia is insufficient. Original faces, branches,
source section gauges and unrelated source components remain present.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse.csgraph import shortest_path

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.commands.route_hand_axis import fixed_context_indices
from geometry_tools.commands.couple_hand_arteries import rotation
from geometry_tools.algorithms.source_profile_axis import source_profile_graph, complete_source_sections
from geometry_tools.algorithms.layered_tube_corridor import layered_tube_corridor,route_station_correspondence
from geometry_tools.algorithms.feature_normals import indexed_feature_normals, inward_sheet_clearance
from geometry_tools.algorithms.tube_axis_map import axis_transport
from geometry_tools.algorithms.tube_sections import arclength, rotation_minimizing_frames, unit
from geometry_tools.open3d_backend import library


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('spec');args=p.parse_args()
    spec_path=Path(args.spec).resolve();spec=json.loads(spec_path.read_text());start=time.monotonic()
    out=Path(spec['output']).resolve()
    if out.exists():raise ValueError('Use a new saved connected-branch stage')
    keys=['source','candidate','profiles','sourceAxis','courseMatches','bones','skin','recipientContext']
    paths={k:Path(spec[k]).resolve()for k in keys};data={k:json.loads(v.read_text())for k,v in paths.items()}
    source={r['partId']:r for r in data['source']['rows']};rows={r['partId']:r for r in data['candidate']['rows']}
    parent,child,root_peer=[spec[k]for k in ['parent','child','rootPeer']]
    axis_record=data['sourceAxis']
    if not axis_record['complete']or axis_record['partId']!=child or axis_record['component']!=spec['component']:
        raise ValueError('One complete explicitly owned source branch required')
    if (sha(axis_record['axisPath'])!=axis_record['axisSha256']or
        axis_record['sourceMeshSha256']!=source[child]['sha256']):raise ValueError('Source branch axis changed')
    with np.load(axis_record['axisPath'],allow_pickle=False)as z:
        child_axis=z['axis'].copy();child_faces=z['source_face_ids'].copy()
    for pid in [parent,child,root_peer]:
        for r in [source[pid],rows[pid]]:
            if sha(r['path'])!=r['sha256']:raise ValueError('Assembly source changed')
        if not np.array_equal(Mesh.load(source[pid]['path']).faces,Mesh.load(rows[pid]['path']).faces):
            raise ValueError('Original source face ancestry changed')
    old_parent=Mesh.load(source[parent]['path']);old_child=Mesh.load(source[child]['path'])
    current_parent=Mesh.load(rows[parent]['path']);current_child=Mesh.load(rows[child]['path'])
    def graph(pid,mesh):
        r=next(r for r in data['profiles']['rows']if r['partId']==pid)
        if sha(r['path'])!=r['sha256']:raise ValueError('Source section bank changed')
        with np.load(r['path'],allow_pickle=False)as z:pairing=z['retained_source_triangle_pairs'].copy()
        return source_profile_graph(mesh.vertices,mesh.faces,triangle_pairing=pairing,extend_annuli=True)
    pg=graph(parent,old_parent);cg=graph(child,old_child)
    pairs=data['source']['pairs'];join=next(r for r in pairs if {r['firstPart'],r['secondPart']}=={parent,child})
    root=next(r for r in pairs if {r['firstPart'],r['secondPart']}=={parent,root_peer})
    def ids(pair,pid):return np.array(pair['firstRaw']if pair['firstPart']==pid else pair['secondRaw'],int)
    root_raw,join_raw,child_join=ids(root,parent),ids(join,parent),ids(join,child)
    if not np.array_equal(old_parent.vertices[join_raw],old_child.vertices[child_join]):
        raise ValueError('The source artery and digit are not an actual shared mouth')
    if not np.allclose(child_axis[0],old_child.vertices[child_join].mean(0),rtol=0,atol=1e-12):
        raise ValueError('Declared source branch axis begins at another mouth')
    if len(pg['boundary_nodes'])!=2 or np.max(np.diff(pg['graph'].indptr))>2:
        raise ValueError('A branched palm parent requires a full network fit')
    def node(raw):
        matches=[n for n,r in zip(pg['boundary_nodes'],pg['boundary_raw'])if set(r)==set(raw)]
        if len(matches)!=1:raise ValueError('Unmatched actual parent mouth')
        return matches[0]
    root_node,end_node=node(root_raw),node(join_raw)
    _,pred=shortest_path(pg['graph'],directed=False,indices=root_node,return_predecessors=True)
    order=[end_node]
    while order[-1]!=root_node:
        k=int(pred[order[-1]])
        if k<0 or k in order:raise ValueError('Invalid original parent section chain')
        order.append(k)
    order=order[::-1]
    if len(order)!=len(pg['centres']):raise ValueError('Unvisited source parent section detail')
    parent_axis=pg['centres'][order]
    root_points=Mesh.load(rows[root_peer]['path']).vertices[ids(root,root_peer)]
    parent_guide=np.array([current_parent.vertices[pg['profiles'][i]].mean(0)for i in order])
    parent_guide[0]=root_points.mean(0)
    match=next(r for r in data['courseMatches']['courseMatches']if r['sourcePartId']==child and r['sourceComponent']==spec['component'])
    child_guide=np.array(match['targetAxis'])
    # The native distal guide may begin slightly proximal to the transported
    # donor end. Drop only those redundant GUIDE stations, never source faces.
    control=np.r_[parent_guide[:-1],child_guide]
    control_s=arclength(control);port_phase=float(control_s[len(parent_guide)-1]/control_s[-1])
    samples=np.linspace(0,control_s[-1],int(np.ceil(control_s[-1]/.001))+1)
    guide=PchipInterpolator(control_s,control)(samples)
    parent_radii=np.array([np.linalg.norm(old_parent.vertices[pg['profiles'][i]]-pg['centres'][i],axis=1).max()for i in order])
    child_radii=np.array([np.linalg.norm(old_child.vertices[child_join]-child_axis[0],axis=1).max(),
        *[r['maximumContourRadius_m']for r in axis_record['proof']['contourWitnesses']]])
    child_radii=np.r_[child_radii,child_radii[-1]]
    phases=np.r_[arclength(parent_axis)/arclength(parent_axis)[-1]*port_phase,
        port_phase+np.linspace(0,1,len(child_radii))[1:]*(1-port_phase)]
    radius=np.interp(samples/samples[-1],phases,np.r_[parent_radii,child_radii[1:]])
    side=rows[parent]['side'];fixed=[r for r in data['bones']['rows']if r['partId'].endswith('.l'if side=='left'else'.r')]
    for r in [*fixed,data['skin'],*data['recipientContext']['rows']]:
        if sha(r['path'])!=r['sha256']:raise ValueError('Fixed recipient changed')
    bone,skin,context,scene=fixed_context_indices(tuple((r['path'],r['sha256'])for r in fixed),
        (data['skin']['path'],data['skin']['sha256']),
        tuple((r['path'],r['sha256'],r['closedOrientedManifold'],r['partId'])for r in data['recipientContext']['rows']),.00002)
    def gap(q):
        result=[]
        for first in range(0,len(q),24000):
            v=q[first:first+24000];hit,_=bone.closest_points(v)
            tensor=library().core.Tensor(((v-bone.origin)/bone.scale).astype(np.float32))
            inside=bone.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5
            value=hit['distance_m']*np.where(inside,-1.,1.)-.0002
            hit,_=skin.closest_points(v);normal,rim=indexed_feature_normals(skin,hit)
            _,clearance,_=inward_sheet_clearance(v,hit,normal,rim)
            result.append(np.minimum(np.minimum(value,clearance-.0002),context.clearance(v)))
        return np.concatenate(result)
    if spec.get('savedRoute'):
        saved=spec['savedRoute'];saved_path=Path(saved['report'])
        if sha(saved_path)!=saved['reportSha256']:raise ValueError('Saved route evidence changed')
        previous=json.loads(saved_path.read_text());proof=previous['connectedBranchRoute']['routeProof']
        if not proof['complete']:raise ValueError('An incomplete route cannot be reused')
        if previous['inputs']!={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()}:
            raise ValueError('Saved route belongs to different actual source/context inputs')
        file=Path(saved['arrays'])
        if sha(file)!=saved['arraysSha256']:raise ValueError('Saved route arrays changed')
        with np.load(file,allow_pickle=False)as z:
            if not np.array_equal(z['guide'],guide)or not np.array_equal(z['radius'],radius):
                raise ValueError('Saved guide or original source radii do not replay literally')
            route=z['route'].copy()
    else:
        route,proof=layered_tube_corridor(guide,radius,gap,segment_clear=scene.segment_clear,
            fixed_prefix=1,fixed_distal=False,grid_m=.0005,half_width_m=.006,
            allow_lateral_detours=True,wall_seconds=16.)
    out.mkdir();np.savez(out/'route.npz',guide=guide,radius=radius,route=np.empty((0,3))if route is None else route)
    evidence=dict(parentPart=parent,childPart=child,childComponent=spec['component'],routeProof=proof,
        wholeConnectedRoute=True,portGuidePhase=port_phase,sourceFacesRetained=True,
        skeletonChanged=False,installed=False,installationAllowed=False,
        inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},sourceCodeSha256=sha(__file__))
    if route is None:
        if 'closestRegionPoints_m'in proof:
            q=np.array(proof['closestRegionPoints_m']);evidence['barrierWitnesses']=context.witnesses(np.r_[q,q.mean(0)[None]])
        evidence['elapsedSeconds']=time.monotonic()-start;write_json(out/'route-failure.json',evidence)
        print(json.dumps(dict(route=proof,elapsedSeconds=evidence['elapsedSeconds'])));return
    def sampled(axis):
        s=arclength(axis);curve=PchipInterpolator(s/s[-1],axis);t=np.linspace(0,1,192)
        return curve(t),unit(curve.derivative()(t))
    (op,opt),(oc,oct)=sampled(parent_axis),sampled(child_axis)
    old=np.r_[op[:-1],oc]
    route_s=arclength(route);curve=PchipInterpolator(route_s/route_s[-1],route)
    target_phase=np.r_[np.linspace(0,port_phase,len(op))[:-1],np.linspace(port_phase,1,len(oc))]
    station_s,station_proof=route_station_correspondence(guide,route,grid_m=.0005,half_width_m=.006)
    mapped_phase=PchipInterpolator(samples/samples[-1],station_s/route_s[-1])(target_phase)
    evidence['sourceStationCorrespondence']=station_proof
    evidence['sourcePortPhaseCorrection_m']=float(np.interp(port_phase,samples/samples[-1],station_s)-port_phase*route_s[-1])
    evidence['savedRouteReused']=spec.get('savedRoute')
    new=curve(mapped_phase);new[0]=root_points.mean(0)
    # Preserve the shared-frame construction by default. Independent digital
    # roll is an explicit experimental option: its v902 replay regressed actual
    # body/self checks and must not become an implicit replacement.
    npart,nchild=new[:len(op)],new[len(op)-1:]
    target_t=unit(curve.derivative()(mapped_phase));pt,ct=target_t[:len(op)].copy(),target_t[len(op)-1:].copy()
    cbefore=rotation_minimizing_frames(oc,oct,initial_reference=np.array([0.,0.,1.]))
    cafter=rotation_minimizing_frames(nchild,ct,initial_reference=np.array([0.,0.,1.]))
    cr=cafter@np.swapaxes(cbefore,1,2)
    pbefore=rotation_minimizing_frames(op,opt);root_R=rotation(old_parent.vertices[root_raw],root_points)
    pt[0]=root_R@opt[0];pt[-1]=cr[0]@opt[-1]
    pafter=rotation_minimizing_frames(npart,pt,initial_reference=root_R@pbefore[0,:,0],
        final_reference=cr[0]@pbefore[-1,:,0])
    pr=pafter@np.swapaxes(pbefore,1,2);pr[0]=root_R;pr[-1]=cr[0]
    if spec.get('independentDigitalFrame',False):
        poses={parent:(op,npart,pr),child:(oc,nchild,cr)}
        evidence['frameConstruction']='Explicit experimental independent digital frame; requires independent acceptance'
    else:
        old_t=unit(np.gradient(old,axis=0));new_t=unit(np.gradient(new,axis=0))
        before=rotation_minimizing_frames(old,old_t);new_t[0]=root_R@old_t[0]
        after=rotation_minimizing_frames(new,new_t,initial_reference=root_R@before[0,:,0])
        rotations=after@np.swapaxes(before,1,2);rotations[0]=root_R
        poses={parent:(old,new,rotations),child:(old,new,rotations)}
        evidence['frameConstruction']='Continuous shared route frame with source-station correspondence retained'
    child_raw=np.unique(old_child.faces[child_faces]);new_meshes={};section_counts={}
    for pid,m,base,g,selected_raw in [(parent,old_parent,current_parent,pg,np.arange(len(old_parent.vertices))),
                                    (child,old_child,current_child,cg,child_raw)]:
        old_part,new_part,part_R=poses[pid]
        q=base.vertices.copy();q[selected_raw],raw_R=axis_transport(m.vertices[selected_raw],old_part,new_part,part_R)
        _,physical=np.unique(m.vertices,axis=0,return_inverse=True);selected_mask=np.zeros(len(q),bool);selected_mask[selected_raw]=True
        sections=complete_source_sections(g,physical);count=0
        for section in sections:
            if not selected_mask[section].any():continue
            if not selected_mask[section].all():raise ValueError('A source section crosses the declared component scope')
            c=m.vertices[section].mean(0);posed,rot=axis_transport(c[None],old_part,new_part,part_R)
            for raw in section:
                aliases=physical==physical[raw];q[aliases]=(m.vertices[aliases]-c)@rot[0].T+posed[0]
            count+=1
        section_counts[pid]=count
        new_meshes[pid]=(m,base,q,physical,selected_mask)
    pm,pb,pq,pi,ps=new_meshes[parent];cm,cb,cq,ci,cs=new_meshes[child]
    for raw,point in zip(root_raw,root_points):pq[pi==pi[raw]]=point
    for a,b in zip(join_raw,child_join):cq[ci==ci[b]]=pq[a]
    output=copy.deepcopy(data['candidate']);reference={r['partId']:r for r in data['source']['rows']}
    for pid,(m,base,q,physical,selected_mask)in new_meshes.items():
        n=np.zeros((physical.max()+1,3));t=q[m.faces];fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        np.add.at(n,physical[m.faces.ravel()],np.repeat(fn,3,axis=0));n=n[physical]
        n/=np.linalg.norm(n,axis=1)[:,None];attrs={k:v.copy()for k,v in base.attributes.items()}
        attrs['normals'][selected_mask]=np.rint(n[selected_mask]*32767).astype('<i2')
        attrs['normalSigns'][selected_mask]=1
        final=Mesh(q,m.faces,m.frame,attrs,m.source);file=out/('parent.npz'if pid==parent else'child.npz');final.save(file)
        row=next(r for r in output['rows']if r['partId']==pid)
        row.update(path=str(file),sha256=sha(file),geometrySha256=final.geometry_sha,
            sourcePath=reference[pid]['path'],sourceSha256=reference[pid]['sha256'])
    output.update(connectedBranchRoute=evidence,complete=False,installed=False,installationAllowed=False,
        portsRequireRecoupling=True,pending=['Upstream radial interfaces, actual all-recipient/containment/source checks and rendered sections'],
        inputs=evidence['inputs'],sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-start)
    if time.monotonic()-start>22:raise TimeoutError('22-second native connected-branch stage')
    write_json(out/'report.json',output)
    print(json.dumps(dict(route=proof,sourceSectionCounts=section_counts,sourceFacesRetained=True,
        elapsedSeconds=output['elapsedSeconds'])))


if __name__=='__main__':main()
