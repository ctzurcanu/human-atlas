"""Replay actual source-cohort triangles, sections and fixed recipient contacts.

Scoped evidence never accepts the whole hand. Moved ports remain an explicit
parent recoupling obligation, even when every selected surface clears context.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts,self_contact_summary
from geometry_tools.course_io import sha,write_json
from geometry_tools.commands.check_hand_network import combine
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,complete_source_sections


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['before','candidate','profiles','binding','bones','skin-reference','recipient-context','part','components','output']:
        p.add_argument('--'+k,required=True)
    args=p.parse_args();start=time.monotonic()
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in
        ['before','candidate','profiles','binding','bones','skin-reference','recipient-context']}
    data={k:json.loads(v.read_text())for k,v in paths.items()};out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new source-cohort audit report')
    before=next(r for r in data['before']['rows']if r['partId']==args.part)
    row=next(r for r in data['candidate']['rows']if r['partId']==args.part)
    for r in [before,row]:
        if sha(r['path'])!=r['sha256']:raise ValueError('Source/candidate arrays changed')
    old,m=Mesh.load(before['path']),Mesh.load(row['path'])
    if not np.array_equal(m.faces,old.faces):raise ValueError('Source face ancestry changed')
    chosen=set(map(int,args.components.split(',')));owners=[r for r in data['binding']['rows']if r['partId']==args.part and r['physicalComponent']in chosen]
    if {r['physicalComponent']for r in owners}!=chosen or not all(r['anatomicallyReviewed']for r in owners):
        raise ValueError('Explicit reviewed component ownership required')
    fids=np.unique(np.concatenate([r['sourceFaceIds']for r in owners]));raw=np.unique(m.faces[fids]);remap=np.full(len(m.vertices),-1,int);remap[raw]=np.arange(len(raw))
    part=Mesh(m.vertices[raw],remap[m.faces[fids]],m.frame)
    fixed=list(data['recipient-context']['rows'])
    bones=[r for r in data['bones']['rows']if r['partId'].endswith('.l'if row['side']=='left'else'.r')]
    for r in [*fixed,*bones,data['skin-reference']]:
        if sha(r['path'])!=r['sha256']:raise ValueError('Fixed recipient arrays changed')
    fixed+=[dict(partId='Fixed '+row['side']+' hand bones',mesh=combine([Mesh.load(r['path'])for r in bones])),
            dict(partId='Visible skin',path=data['skin-reference']['path'])]
    contacts=[]
    for r in fixed:
        reference=r['mesh']if 'mesh'in r else Mesh.load(r['path'])
        if reference.frame!=m.frame:raise ValueError('Fixed recipient and vessel have different coordinate frames')
        hit,proof=pair_contacts(part,reference,0.,max_tests=300000,max_pairs=50000)
        if not proof['complete_within_pair_budget']:raise ValueError('Incomplete recipient triangle replay')
        wrong=np.flatnonzero(np.isin(hit['relation'],['crossing','coplanar-overlap']))
        contacts.append(dict(partId=r['partId'],pairs=len(wrong),
            examples=[dict(sourceFace=int(fids[hit['pairs'][i,0]]),recipientFace=int(hit['pairs'][i,1]))for i in wrong[:12]]))
        if time.monotonic()-start>20:raise TimeoutError('20-second native component-recipient audit')
    self_selected,_=self_contact_summary(part);self_whole,_=self_contact_summary(m)
    if not self_selected['complete_within_pair_budget']or not self_whole['complete_within_pair_budget']:
        raise ValueError('Incomplete source self-contact replay')
    bank=next(r for r in data['profiles']['rows']if r['partId']==args.part)
    if sha(bank['path'])!=bank['sha256']:raise ValueError('Source profile pairing changed')
    with np.load(bank['path'],allow_pickle=False)as z:pairing=z['retained_source_triangle_pairs'].copy()
    graph=source_profile_graph(old.vertices,old.faces,triangle_pairing=pairing,extend_annuli=True)
    physical=np.unique(old.vertices,axis=0,return_inverse=True)[1]
    sections=complete_source_sections(graph,physical)
    gauge=max((float(np.abs(np.linalg.norm(m.vertices[s][:,None]-m.vertices[s][None,:],axis=2)-
                            np.linalg.norm(old.vertices[s][:,None]-old.vertices[s][None,:],axis=2)).max())for s in sections),default=0.)
    unselected=np.setdiff1d(np.arange(len(m.vertices)),raw)
    literal=np.array_equal(m.vertices[unselected],old.vertices[unselected])and all(
        np.array_equal(a[unselected],old.attributes[k][unselected])for k,a in m.attributes.items()if len(a)==len(m.vertices))
    peers={r['partId']:r for r in data['candidate']['rows']};moved_ports=[]
    for pair in data['candidate']['pairs']:
        if args.part not in [pair['firstPart'],pair['secondPart']]:continue
        keys=['first','second']if pair['firstPart']==args.part else ['second','first']
        raw_ids=np.array(pair[keys[0]+'Raw']);peer=Mesh.load(peers[pair[keys[1]+'Part']]['path']);other_ids=np.array(pair[keys[1]+'Raw'])
        gap=np.linalg.norm(m.vertices[raw_ids]-peer.vertices[other_ids],axis=1).max()
        if gap>0:moved_ports.append(dict(peer=pair[keys[1]+'Part'],maximumGap_m=float(gap),sourceRaw=raw_ids.tolist()))
    record=dict(complete=True,partId=args.part,components=sorted(chosen),sourceFaces=len(fids),contacts=contacts,
        ownComponentSelfPairs=self_selected['defect_pairs'],wholePartSelfPairs=self_whole['defect_pairs'],
        sourceFacesLiteral=True,unselectedArraysLiteral=bool(literal),sourceSectionGaugeError_m=gauge,
        sourceSectionCount=len(sections),movedPorts=moved_ports,parentsRequireRecoupling=bool(moved_ports),
        selectedSurfacesClear=not any(r['pairs']for r in contacts)and self_selected['defect_pairs']==0,
        continuousContainmentNotEstablished=True,accepted=False,installed=False,installationAllowed=False,
        inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},sourceCodeSha256=sha(__file__),
        elapsedSeconds=time.monotonic()-start)
    write_json(out,record)
    compact={k:v for k,v in record.items()if k not in ['inputs','sourceCodeSha256','contacts']}
    compact.update(fixedSourcesChecked=len(contacts),nonzeroContacts=[r for r in contacts if r['pairs']])
    print(json.dumps(compact))


if __name__=='__main__':main()
