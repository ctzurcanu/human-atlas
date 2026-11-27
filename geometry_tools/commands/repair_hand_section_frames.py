"""Correct source section frames on the existing saved palm courses.

Station centres and declared mouth coordinates are retained. Complete source
section polygons use proper transported rotations; original branch cores are
reconstructed with bounded ARAP. No skeleton or source branch is replaced.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,complete_source_sections
from geometry_tools.algorithms.profile_graph_transport import transport_linked_profiles
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores
from geometry_tools.algorithms.port_section_graph import shortest_rotation
from geometry_tools.commands.couple_hand_arteries import rotation


def scoped_route_frames(args,paths,data,out,started):
    """Orient literal source closures on a saved branch-aware centre route."""
    from scipy.sparse import coo_matrix
    from geometry_tools.algorithms.source_control_clearance import source_control_basis
    from geometry_tools.algorithms.tree_section_corridor import source_closure_tree
    from geometry_tools.algorithms.port_section_graph import pose_section_graph
    sr=next(r for r in data['source']['rows']if r['partId']==args.part)
    row=next(r for r in data['candidate']['rows']if r['partId']==args.part)
    for r in [sr,row]:
        if sha(r['path'])!=r['sha256']:raise ValueError('Source route changed')
    source,current=Mesh.load(sr['path']),Mesh.load(row['path']);original=Mesh.load(sr['sourcePath'])
    if sha(sr['sourcePath'])!=sr['sourceSha256']or not np.array_equal(original.faces,source.faces)or not np.array_equal(source.faces,current.faces):
        raise ValueError('Literal source face ancestry changed')
    bank=next(r for r in data['profiles']['rows']if r['partId']==args.part)
    if sha(bank['path'])!=bank['sha256']:raise ValueError('Source pairing changed')
    with np.load(bank['path'],allow_pickle=False)as z:pairing=z['retained_source_triangle_pairs'].copy()
    g=source_profile_graph(original.vertices,original.faces,triangle_pairing=pairing,extend_annuli=True)
    _,physical=np.unique(original.vertices,axis=0,return_inverse=True);sections=complete_source_sections(g,physical)
    W,_,_,_=source_control_basis(original.vertices,original.faces,sections,physical,micro_coherence_m=.00002)
    owner,edges,original_groups=source_closure_tree(original.faces,np.argmax(W,axis=1));count=int(owner.max())+1
    before=np.array([source.vertices[owner==k].mean(0)for k in range(count)])
    after=np.array([current.vertices[owner==k].mean(0)for k in range(count)])
    delta=after-before
    if np.max(abs(source.vertices+delta[owner]-current.vertices))>1e-10:
        raise ValueError('Saved route does not preserve its original rigid source closures')
    a,b=edges.T;length=np.linalg.norm(before[a]-before[b],axis=1)
    graph=coo_matrix((np.r_[length,length],(np.r_[a,b],np.r_[b,a])),shape=(count,count)).tocsr()
    _,R,proof=pose_section_graph(dict(centres=before,graph=graph),{i:after[i]for i in range(count)},{})
    q=after[owner]+np.einsum('nij,nj->ni',R[owner],source.vertices-before[owner])
    gauge=max(float(np.abs(np.linalg.norm(q[s,None]-q[None,s],axis=2)-np.linalg.norm(source.vertices[s,None]-source.vertices[None,s],axis=2)).max())for s in sections)
    if gauge>1e-12:raise ValueError('A literal source section calibre changed')
    _,inv=np.unique(physical,return_inverse=True);normal=np.zeros((inv.max()+1,3));t=q[source.faces]
    fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);np.add.at(normal,inv[source.faces.ravel()],np.repeat(fn,3,axis=0));normal=normal[inv]
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    final=Mesh(q,source.faces,source.frame,dict(current.attributes,normals=np.rint(normal*32767).astype('<i2'),normalSigns=np.ones(len(q),'i1')),current.source)
    out.mkdir();file=out/'framed-part.npz';final.save(file);report=copy.deepcopy(data['candidate'])
    r=next(r for r in report['rows']if r['partId']==args.part);r.update(path=str(file),sha256=sha(file),geometrySha256=final.geometry_sha)
    for pair in report['pairs']:
        if args.part in [pair['firstPart'],pair['secondPart']]:pair.update(literalCoordinatesAndEffectiveNormals=False,packedCoordinatesAndEffectiveNormalWordsLiteral=False,requiresRecoupling=True)
    report.update(sectionFrameRepair=dict(partId=args.part,rigidClosures=count,originalSourceGroups=original_groups,sourceSections=len(sections),calibreError_m=gauge,frames=proof),
                  sourceFacesLiteral=True,sourceSectionGaugeError_m=gauge,portsRequireRecoupling=True,complete=False,installed=False,installationAllowed=False,
                  inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
    if report['elapsedSeconds']>18:raise TimeoutError('18-second source-route frame stage')
    write_json(out/'report.json',report);print(json.dumps(dict(closures=count,sections=len(sections),sourceFacesRetained=True,elapsedSeconds=report['elapsedSeconds'])))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','candidate','profiles','output']:
        p.add_argument('--'+k,required=True)
    p.add_argument('--part',help='Only orient this existing source branch route; adjoining mouths require recoupling')
    args=p.parse_args(); started=time.monotonic()
    paths={k:Path(getattr(args,k)).resolve()for k in ['source','candidate','profiles']}
    data={k:json.loads(v.read_text())for k,v in paths.items()}
    out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new frame-repair stage')
    if args.part:
        scoped_route_frames(args,paths,data,out,started)
        return
    originals={r['partId']:r for r in data['source']['rows']}
    rows={r['partId']:r for r in data['candidate']['rows']}
    result={pid:Mesh.load(r['path'])for pid,r in rows.items()}
    evidence=[]
    for changed in data['candidate']['palmParentCoupling']:
        pid=changed['partId']; row=rows[pid]; sr=originals[pid]
        for r in [row,sr]:
            if sha(r['path'])!=r['sha256']:raise ValueError('Source binding changed')
        source,current=Mesh.load(sr['path']),result[pid]
        if not np.array_equal(source.faces,current.faces):raise ValueError('Face ancestry changed')
        bank=next((r for r in data['profiles']['rows']if r['partId']==pid),None)
        if bank:
            if sha(bank['path'])!=bank['sha256']:raise ValueError('Recorded section bank changed')
            with np.load(bank['path'],allow_pickle=False)as z:pairing=z['retained_source_triangle_pairs'].copy()
        else:pairing=np.arange(len(source.faces)//2*2).reshape(-1,2)
        g=source_profile_graph(source.vertices,source.faces,triangle_pairing=pairing,extend_annuli=True)
        _,rep,inv=np.unique(source.vertices,axis=0,return_index=True,return_inverse=True)
        q=g['centres'].copy(); R=np.repeat(np.eye(3)[None],len(q),axis=0)
        for node in range(len(q)):
            raw=rep[g['raw_node'][rep]==node]
            if len(raw):
                q[node]=current.vertices[raw].mean(0)
                if len(raw)>=3:R[node]=rotation(source.vertices[raw],current.vertices[raw])
        protected=set(); boundary={}
        for node,loop in zip(g['boundary_nodes'],g['boundary_raw']):
            q[node]=current.vertices[loop].mean(0)
            R[node]=rotation(source.vertices[loop],current.vertices[loop]); protected.add(node)
            for raw in loop:boundary[int(inv[raw])]=current.vertices[raw].copy()
        axis_before=[]
        for i,ring in enumerate(g['profiles']):
            R[i]=rotation(source.vertices[ring],current.vertices[ring])
            if i in protected:continue
            neighbours=g['graph'].indices[g['graph'].indptr[i]:g['graph'].indptr[i+1]]
            if len(neighbours)==2:
                dx=g['centres'][neighbours[1]]-g['centres'][neighbours[0]]
                dy=q[neighbours[1]]-q[neighbours[0]]
                axis_before.append(float((R[i]@dx)@dy/np.linalg.norm(dx)/np.linalg.norm(dy)))
                R[i]=shortest_rotation(R[i]@dx,dy)@R[i]
        a,b=g['graph'].nonzero(); links=np.array([(i,j)for i,j in zip(a,b)if i<j],int).reshape(-1,2)
        R,frames=transport_linked_profiles(g['centres'],q,R,links,protected=protected)
        keys=np.array(sorted(boundary),int);raw=rep[keys];targets=np.array([boundary[int(k)]for k in keys])
        chosen=[i for i,s in enumerate(g['profiles'])if not np.isin(inv[s],keys).any()]
        profiles=[g['profiles'][i]for i in chosen]
        fitted=transport_profile_cores(source.vertices,source.faces,profiles,q[chosen],R[chosen],
                    physical_ids=inv,prescribed_vertices=raw,prescribed_positions=targets,
                    core_method='arap',core_wall_seconds=.4)
        v=fitted['vertices']
        for k,target in zip(keys,targets):v[inv==k]=target
        for i,section in zip(chosen,profiles):
            if np.linalg.norm(v[section].mean(0)-q[i])>1e-12:raise ValueError('A saved section centre moved')
        normals=np.zeros((len(rep),3));t=v[source.faces];fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        np.add.at(normals,inv[source.faces.ravel()],np.repeat(fn,3,axis=0));normals=normals[inv]
        normals/=np.linalg.norm(normals,axis=1)[:,None]
        attrs=dict(current.attributes,normals=np.rint(normals*32767).astype('<i2'),normalSigns=np.ones(len(v),'i1'))
        result[pid]=Mesh(v,source.faces,source.frame,attrs,current.source)
        evidence.append(dict(partId=pid,sections=len(chosen),savedCentresRetained=True,
                    reversedFramesBefore=sum(x<0 for x in axis_before),frames=frames,
                    maximumMovement_m=float(np.linalg.norm(v-current.vertices,axis=1).max())))
        if time.monotonic()-started>22:raise TimeoutError('22-second frame repair stage')
    attrs={pid:{k:v.copy()for k,v in m.attributes.items()}for pid,m in result.items()}
    pairs=copy.deepcopy(data['candidate']['pairs'])
    for link in pairs:
        pa,pb=link['firstPart'],link['secondPart'];a,b=np.array(link['firstRaw']),np.array(link['secondRaw'])
        if not np.array_equal(result[pa].vertices[a],result[pb].vertices[b]):raise ValueError('Declared mouth moved')
        canonical=link['canonicalPart'];ids=a if canonical==pa else b
        words=attrs[canonical]['normals'][ids].astype(np.int32)*attrs[canonical]['normalSigns'][ids,None]
        for pid,indices in [(pa,a),(pb,b)]:
            for index,word in zip(indices,words):
                mask=np.all(result[pid].vertices==result[pid].vertices[index],axis=1)
                attrs[pid]['normals'][mask]=word;attrs[pid]['normalSigns'][mask]=1
    out.mkdir(); emitted=[]
    for pid,m in result.items():
        final=Mesh(m.vertices,m.faces,m.frame,attrs[pid],m.source);path=out/f'part-{len(emitted):02}.npz';final.save(path)
        emitted.append(dict(rows[pid],path=str(path),sha256=sha(path),geometrySha256=final.geometry_sha))
    report=copy.deepcopy(data['candidate']);report.update(rows=emitted,pairs=pairs,
        sectionFrameRepair=evidence,constructionComplete=True,complete=False,installed=False,installationAllowed=False,
        inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json',report)
    print(json.dumps(dict(parents=len(evidence),reversedFramesBefore=sum(x['reversedFramesBefore']for x in evidence),elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
