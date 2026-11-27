"""Route explicitly selected source components through finite section corridors."""
import argparse, copy, json, time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,complete_source_sections
from geometry_tools.algorithms.source_control_clearance import source_control_basis
from geometry_tools.algorithms.tree_section_corridor import tree_section_corridor
from geometry_tools.commands.check_hand_network import combine


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['candidate','profiles','binding','bones','skin-reference','part','components','output']:
        p.add_argument('--'+k,required=True)
    p.add_argument('--recipient-context')
    p.add_argument('--screen-spans',action='store_true')
    p.add_argument('--retain-source-ports',action='store_true')
    p.add_argument('--guide-route')
    p.add_argument('--grid-spacing-mm',type=float,default=1.)
    p.add_argument('--local-radius-mm',type=float)
    args=p.parse_args();started=time.monotonic()
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in
           ['candidate','profiles','binding','bones','skin-reference']}
    data={k:json.loads(v.read_text())for k,v in paths.items()};out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new source routing stage')
    row=next(r for r in data['candidate']['rows']if r['partId']==args.part)
    if sha(row['path'])!=row['sha256']:raise ValueError('Proposal changed')
    mesh=Mesh.load(row['path']);source=Mesh.load(row['sourcePath'])
    guide=None
    if args.guide_route:
        path=Path(args.guide_route).resolve();saved=json.loads(path.read_text());paths['guide-route']=path
        if not saved.get('finiteSourceRoutes')or not all(r['complete']for r in saved['finiteSourceRoutes']):raise ValueError('Saved route did not complete its native algorithm stage')
        gr=next(r for r in saved['rows']if r['partId']==args.part)
        if sha(gr['path'])!=gr['sha256']:raise ValueError('Saved route guide changed')
        gm=Mesh.load(gr['path'])
        if not np.array_equal(gm.faces,mesh.faces):raise ValueError('Guide source ancestry changed')
        guide=gm.vertices
    if sha(row['sourcePath'])!=row['sourceSha256']or not np.array_equal(source.faces,mesh.faces):
        raise ValueError('Source ancestry changed')
    bank=next(r for r in data['profiles']['rows']if r['partId']==args.part)
    if sha(bank['path'])!=bank['sha256']:raise ValueError('Section bank changed')
    with np.load(bank['path'],allow_pickle=False)as z:pairing=z['retained_source_triangle_pairs'].copy()
    graph=source_profile_graph(source.vertices,source.faces,triangle_pairing=pairing,extend_annuli=True)
    _,physical=np.unique(source.vertices,axis=0,return_inverse=True)
    sections=complete_source_sections(graph,physical)
    W,_,_,_=source_control_basis(source.vertices,source.faces,sections,physical,micro_coherence_m=.00002)
    labels=np.argmax(W,axis=1)
    bones=[r for r in data['bones']['rows']if r['partId'].endswith('.l'if row['side']=='left'else'.r')]
    skin=data['skin-reference']
    for r in [*bones,skin]:
        if sha(r['path'])!=r['sha256']:raise ValueError('Fixed context changed')
    bone=SurfaceIndex(combine([Mesh.load(r['path'])for r in bones]));skin_index=SurfaceIndex(Mesh.load(skin['path']))
    recipient=[]
    if args.recipient_context:
        path=Path(args.recipient_context).resolve();context=json.loads(path.read_text());paths['recipient-context']=path
        for fixed in context['rows']:
            if sha(fixed['path'])!=fixed['sha256']:raise ValueError('Fixed recipient context changed')
            recipient.append((SurfaceIndex(Mesh.load(fixed['path'])),fixed['closedOrientedManifold'],fixed['partId']))
    q=mesh.vertices.copy();proofs=[];scope=set(map(int,args.components.split(',')))
    port_raw=[]
    if args.retain_source_ports:
        for pair in data['candidate']['pairs']:
            if pair['firstPart']==args.part:port_raw.extend(pair['firstRaw'])
            if pair['secondPart']==args.part:port_raw.extend(pair['secondRaw'])
    port_raw=np.unique(np.asarray(port_raw,int))
    span_context=SurfaceIndex(combine([bone.mesh,skin_index.mesh,*[index.mesh for index,_,_ in recipient]]))if args.screen_spans else None
    for owner in data['binding']['rows']:
        if owner['partId']!=args.part or owner['physicalComponent']not in scope:continue
        if not owner['anatomicallyReviewed']:raise ValueError('Unreviewed source owner')
        fids=np.array(owner['sourceFaceIds'],int);raw=np.unique(mesh.faces[fids]);remap=np.full(len(q),-1,int);remap[raw]=np.arange(len(raw))
        moved,proof=tree_section_corridor(q[raw],remap[mesh.faces[fids]],labels[raw],bone,skin_index,
                                        recipient_surfaces=recipient,span_context=span_context,
                                        guide_vertices=None if guide is None else guide[raw],grid_spacing_m=args.grid_spacing_mm*.001,
                                        protected_raw=remap[port_raw][remap[port_raw]>=0],
                                        local_radius_m=None if args.local_radius_mm is None else args.local_radius_mm*.001)
        proof['physicalComponent']=owner['physicalComponent'];proofs.append(proof)
        if proof['complete']:q[raw]=moved
        if time.monotonic()-started>20:raise TimeoutError('20-second source-corridor stage')
    _,inv=np.unique(source.vertices,axis=0,return_inverse=True);normal=np.zeros((inv.max()+1,3));t=q[mesh.faces]
    fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);np.add.at(normal,inv[mesh.faces.ravel()],np.repeat(fn,3,axis=0));normal=normal[inv]
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    final=Mesh(q,mesh.faces,mesh.frame,dict(mesh.attributes,normals=np.rint(normal*32767).astype('<i2'),normalSigns=np.ones(len(q),'i1')),mesh.source)
    out.mkdir();file=out/'routed-part.npz';final.save(file)
    report=copy.deepcopy(data['candidate']);r=next(r for r in report['rows']if r['partId']==args.part)
    r.update(path=str(file),sha256=sha(file),geometrySha256=final.geometry_sha)
    for pair in report.get('pairs',[]):
        if args.part in [pair['firstPart'],pair['secondPart']]:
            pair.update(literalCoordinatesAndEffectiveNormals=False,packedCoordinatesAndEffectiveNormalWordsLiteral=False,requiresRecoupling=True)
    report.update(finiteSourceRoutes=proofs,complete=False,installed=False,installationAllowed=False,
                  portsRequireRecoupling=True,sourceAnnularPropagation=True,
                  inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},
                  sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json',report);print(json.dumps(dict(routes=proofs,elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
