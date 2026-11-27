"""Repair a complete hand as one declared, source-connected control bank.

Shared source mouths move together. Source rings retain congruence, proximal
forearm controls remain fixed, and per-owner context counts cannot regress.
This command never installs anatomy.
"""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.spatial import BVH
from geometry_tools.commands.check_hand_network import combine
from geometry_tools.algorithms.source_profile_axis import source_profile_graph
from geometry_tools.algorithms.monotone_uncross import monotone_uncross
from geometry_tools.source_profile_io import load_source_pairing


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','candidate','profiles','bones','skin-reference','side','output']:
        p.add_argument('--'+k,required=True)
    p.add_argument('--body-clearance',action='store_true',
                   help='Reduce actual fixed-bone contacts while preserving the complete source-connected hand')
    p.add_argument('--coupled-sections',action='store_true',
                   help='Allow opposite half-motions of two contacting original source sections')
    p.add_argument('--fixed-parts',default='',
                   help='Comma-separated complete source part IDs whose coordinates must remain literal')
    args=p.parse_args();started=time.monotonic()
    if args.side not in ['left','right']:raise ValueError('Explicit hand side required')
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in['source','candidate','profiles','bones','skin-reference']}
    data={k:json.loads(v.read_text())for k,v in paths.items()};out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new shared-closure stage')
    rows=[r for r in data['candidate']['rows']if r['side']==args.side]
    source_rows={r['partId']:r for r in data['source']['rows']};original=[];current=[]
    for r in rows:
        s=source_rows[r['partId']]
        if sha(r['path'])!=r['sha256']or sha(s['path'])!=s['sha256']:raise ValueError('Source/proposal changed')
        a,b=Mesh.load(s['path']),Mesh.load(r['path'])
        if not np.array_equal(a.faces,b.faces):raise ValueError('Source face ancestry changed')
        original.append(a);current.append(b)
    source=combine(original);mesh=combine(current)
    _,rep,physical=np.unique(source.vertices,axis=0,return_index=True,return_inverse=True)
    if not np.array_equal(mesh.vertices,mesh.vertices[rep][physical]):
        raise ValueError('Source alias/declared-mouth coordinates must be literal before whole-closure repair')
    owner=np.concatenate([np.full(len(m.faces),i,int)for i,m in enumerate(current)])
    sections=set();offset=0
    for r,s in zip(rows,original):
        bank=next((v for v in data['profiles']['rows']if v['partId']==r['partId']),None)
        pairing=load_source_pairing(bank,s.faces)
        graph=source_profile_graph(s.vertices,s.faces,triangle_pairing=pairing,extend_annuli=True)
        for cycle in graph['profiles']+graph['boundary_raw']:
            sections.add(tuple(sorted(set(physical[offset+cycle].tolist()))))
        offset+=len(s.vertices)
    groups=[]
    for section in map(set,sorted(sections)):
        touching=[i for i,g in enumerate(groups)if g.intersection(section)]
        for i in reversed(touching):section.update(groups.pop(i))
        groups.append(section)
    groups=[sorted(g)for g in groups]
    fixed_parts=set(filter(None,args.fixed_parts.split(',')))
    if not fixed_parts.issubset({r['partId']for r in rows}):
        raise ValueError('Every fixed part must belong to the selected hand closure')
    fixed_raw=[];offset=0
    for row,m in zip(rows,current):
        if row['partId']in fixed_parts:fixed_raw.extend(range(offset,offset+len(m.vertices)))
        offset+=len(m.vertices)
    protected=np.unique(np.r_[physical[source.vertices[:,1]>=1.08],
                              physical[np.asarray(fixed_raw,dtype=int)]])
    bones=[r for r in data['bones']['rows']if r['partId'].endswith('.l'if args.side=='left'else'.r')]
    skin=data['skin-reference']
    for r in [*bones,skin]:
        if sha(r['path'])!=r['sha256']:raise ValueError('Fixed context changed')
    if not skin['defaultHiddenExcluded']:raise ValueError('Actual active skin required')
    contexts=[BVH(combine([Mesh.load(r['path'])for r in bones])),BVH(Mesh.load(skin['path']))]
    if args.body_clearance:
        from geometry_tools.algorithms.monotone_context_clearance import monotone_context_clearance
        from geometry_tools.open3d_backend import SurfaceIndex
        q,proof=monotone_context_clearance(mesh,physical,[rep[g]for g in groups],protected,
              contexts[0],SurfaceIndex(contexts[0].mesh),contexts[1:],
              wall_seconds=12.,maximum_trials=300,maximum_accepted_steps=2,
              connected_contact_groups=True)
    else:
        q,proof=monotone_uncross(mesh,source,physical,protected,contexts,
              wall_seconds=12.,maximum_trials=120 if args.coupled_sections else 300,preserve_context_counts=True,
              rigid_groups=groups,rigid_motion_m=.001,maximum_accepted_steps=2 if args.coupled_sections else 6,
              face_owners=owner,defer_full_replay=True,coupled_groups=args.coupled_sections)
    for g in groups:
        raw=rep[g]
        if np.max(np.abs(np.linalg.norm(q[raw][:,None]-q[raw][None],axis=2)-
                         np.linalg.norm(mesh.vertices[raw][:,None]-mesh.vertices[raw][None],axis=2)))>1e-12:
            raise ValueError('Original source-ring or shared-mouth calibre changed')
    if fixed_raw and not np.array_equal(q[fixed_raw],mesh.vertices[fixed_raw]):
        raise ValueError('A fixed finger vessel moved during connector repair')
    normals=np.zeros((len(rep),3));t=q[mesh.faces];fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    np.add.at(normals,physical[mesh.faces.ravel()],np.repeat(fn,3,axis=0));normals=normals[physical]
    length=np.linalg.norm(normals,axis=1)
    if (length<=1e-20).any():raise ValueError('Undefined source normal')
    normals=np.rint(normals/length[:,None]*32767).astype('<i2')
    out.mkdir();changed={};offset=0
    for r,m in zip(rows,current):
        sl=slice(offset,offset+len(m.vertices));offset+=len(m.vertices)
        final=Mesh(q[sl],m.faces,m.frame,dict(m.attributes,normals=normals[sl],normalSigns=np.ones(len(m.vertices),'i1')),m.source)
        file=out/f'part-{len(changed):02}.npz';final.save(file)
        changed[r['partId']]=dict(r,path=str(file),sha256=sha(file),geometrySha256=final.geometry_sha)
    report=copy.deepcopy(data['candidate']);report.update(rows=[changed.get(r['partId'],r)for r in report['rows']],
          sharedClosureRepair=dict(side=args.side,proof=proof,rigidSourceGroups=len(groups),
                                   fixedProximalControls=len(protected),
                                   fixedParts=sorted(fixed_parts),fixedPartCoordinatesLiteral=True,
                                   objective='fixed-bone-clearance'if args.body_clearance else 'source-crossings',
                                   perOwnerContextNonincrease=not args.body_clearance),
          complete=False,installed=False,installationAllowed=False,skeletonChanged=False,
          inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},
          sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
    for pair in report['pairs']:
        if args.side+'_'in pair['firstPart']:
            pair.update(literalCoordinatesAndEffectiveNormals=True,coordinatesLiteral=True,
                        effectiveNormalWordsRequireSync=False,requiresRecoupling=False,
                        packedCoordinatesAndEffectiveNormalWordsLiteral=False)
    write_json(out/'report.json',report)
    print(json.dumps(dict(side=args.side,before=proof['initialSelfPairs'],after=proof['finalSelfPairs'],
                          targetBefore=proof.get('initialTargetPairs'),targetAfter=proof.get('finalTargetPairs'),
                          complete=proof['complete'],timedOut=proof['timedOut'],elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
