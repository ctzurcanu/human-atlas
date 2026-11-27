"""Repair only declared source cores; pin original profiles, ports and skin contacts.

Source ARAP supplies bounded shape-restoring proposals. Every proposal is
checked for actual intrinsic crossings and actual bone/skin triangle contacts.
No branch, face or recorded source section is removed or shrunk.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.contact import self_contact_summary
from geometry_tools.spatial import BVH
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.arap_core import arap_core
from geometry_tools.algorithms.source_profile_axis import source_profile_graph
from geometry_tools.commands.check_hand_network import combine
from geometry_tools.commands.filter_hand_components import crossings


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'candidate', 'profiles', 'bones', 'skin-reference', 'part', 'output']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--slide-skin-contacts', action='store_true')
    p.add_argument('--strict-local-uncross', action='store_true')
    p.add_argument('--strict-section-uncross', action='store_true')
    p.add_argument('--strict-body-uncross', action='store_true')
    p.add_argument('--fixed-mouth-corridor',action='store_true',
                   help='Solve the complete source part in the bone/skin corridor with every mouth held literal')
    p.add_argument('--allow-source-curvature-change',action='store_true',
                   help='For a complete annular conduit, allow surface bend signs to change; retain rings, face orientation and actual collision gates')
    p.add_argument('--rigid-section-groups', action='store_true')
    p.add_argument('--declared-core-node', type=int,
                   help='Restore only this explicit unpartitioned source-core node; retain all section/port controls')
    p.add_argument('--maximum-core-motion-mm', type=float, default=3.)
    args = p.parse_args()
    if not 0 < args.maximum_core_motion_mm <= 10:
        p.error('Declared source-core motion bound must be within (0,10] mm')
    started = time.monotonic()
    paths = {k: Path(getattr(args, k.replace('-', '_'))).resolve()
             for k in ['source', 'candidate', 'profiles', 'bones', 'skin-reference']}
    data = {k: json.loads(v.read_text()) for k, v in paths.items()}
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new repair stage')
    source = next(r for r in data['source']['rows'] if r['partId'] == args.part)
    row = next(r for r in data['candidate']['rows'] if r['partId'] == args.part)
    bank = next((r for r in data['profiles']['rows'] if r['partId'] == args.part),None)
    for r in [source, row]+([bank]if bank else[]):
        if sha(r['path']) != r['sha256']:
            raise ValueError('Source binding changed')
    rest, mesh = Mesh.load(source['path']), Mesh.load(row['path'])
    if not np.array_equal(rest.faces, mesh.faces):
        raise ValueError('Core repair requires identical source face ancestry')
    original = Mesh.load(row['sourcePath'])
    if bank:
        with np.load(bank['path'], allow_pickle=False) as z:pairs=z['retained_source_triangle_pairs'].copy()
    else:pairs=np.arange(len(original.faces)//2*2).reshape(-1,2)
    graph = source_profile_graph(original.vertices, original.faces, triangle_pairing=pairs,
                                 extend_annuli=data['candidate'].get('sourceAnnularPropagation', False))
    _, rep, inv = np.unique(original.vertices, axis=0, return_index=True, return_inverse=True)
    protected = np.unique(inv[np.concatenate(graph['profiles']+graph['boundary_raw'])])
    edges = np.unique(np.sort(inv[mesh.faces][:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1), axis=0)
    skin_row = data['skin-reference']
    fixed_rows = [r for r in data['bones']['rows'] if r['partId'].endswith('.l' if row['side'] == 'left' else '.r')]
    for r in [skin_row, *fixed_rows]:
        if sha(r['path']) != r['sha256']:
            raise ValueError('Fixed context changed')
    skin_mesh = Mesh.load(skin_row['path'])
    skin = SurfaceIndex(skin_mesh)
    skin_tree = BVH(skin_mesh)
    bone_mesh = combine([Mesh.load(r['path']) for r in fixed_rows])
    bone_tree = BVH(bone_mesh)
    bone_index = SurfaceIndex(bone_mesh)
    q = mesh.vertices.copy()
    stages = []
    first, _ = self_contact_summary(mesh, example_limit=1000)
    if not first['complete_within_pair_budget']:
        raise ValueError('Intrinsic contact inventory incomplete')
    defects = first['defect_pairs']
    strict_proof = None
    if args.fixed_mouth_corridor:
        from geometry_tools.algorithms.source_control_clearance import clear_source_controls
        protected=np.unique(inv[np.concatenate(graph['boundary_raw'])])
        raw=np.flatnonzero(np.isin(inv,protected))
        if args.allow_source_curvature_change:
            covered=np.concatenate([inv[s]for s in graph['profiles']])
            if (len(graph['boundary_raw'])!=2 or len(np.unique(covered))!=len(rep)
                    or len(covered)!=len(rep)):
                raise ValueError('Curvature changes require a complete source annular conduit')
        proposed,strict_proof=clear_source_controls(mesh.vertices,mesh.faces,graph['profiles'],inv,
                    bone_index,skin,maximum_motion_m=args.maximum_core_motion_mm*.001,
                    wall_seconds=2.,protected_raw=raw,
                    preserve_dihedral_sign=not args.allow_source_curvature_change)
        if strict_proof['complete']:
            summary,_=self_contact_summary(Mesh(proposed,mesh.faces,mesh.frame),example_limit=1000)
            counts=[]
            contexts=[bone_tree,skin_tree]+[BVH(Mesh.load(r['path']))for r in data['candidate']['rows']
                           if r['side']==row['side']and r['partId']!=row['partId']]
            for context in contexts:
                before,after=[crossings(v,context,started)for v in [mesh,Mesh(proposed,mesh.faces,mesh.frame)]]
                counts.append(dict(before=before,after=after))
            accepted=bool(summary['complete_within_pair_budget']and summary['defect_pairs']<=defects
                          and all(c['after']<=c['before']for c in counts)
                          and counts[0]['after']==0 and counts[1]['after']==0)
            strict_proof.update(actualTriangleReplay=counts,wholePartSelfPairs=summary['defect_pairs'],
                                actualTriangleGatePassed=accepted)
            if accepted:q,defects=proposed,summary['defect_pairs']
            else:strict_proof['complete']=False
    elif args.strict_body_uncross:
        from geometry_tools.algorithms.monotone_context_clearance import monotone_context_clearance
        contexts=[skin_tree]
        for neighbour in data['candidate']['rows']:
            if neighbour['side']==row['side'] and neighbour['partId']!=row['partId']:
                if sha(neighbour['path'])!=neighbour['sha256']:raise ValueError('Adjoining artery changed')
                contexts.append(BVH(Mesh.load(neighbour['path'])))
        protected=np.unique(inv[np.concatenate(graph['boundary_raw'])])
        q,strict_proof=monotone_context_clearance(mesh,inv,graph['profiles'],protected,bone_tree,bone_index,contexts)
        defects=strict_proof['finalSelfPairs']
    elif args.declared_core_node is not None:
        if args.declared_core_node < len(graph['profiles']):
            raise ValueError('A declared core cannot be an original section polygon')
        raw = np.flatnonzero(graph['raw_node'] == args.declared_core_node)
        free = np.setdiff1d(np.unique(inv[raw]), protected)
        if not len(free):
            raise ValueError('Declared source core has no unprotected controls')
        boundary = np.setdiff1d(np.arange(len(rep)), free)
        # Initialize the free patch from its valid source shape, then solve
        # its actual prescribed border. Do not fit a new course or radius.
        a, b = rest.vertices[raw], mesh.vertices[raw]
        u, _, vt = np.linalg.svd((a-a.mean(0)).T@(b-b.mean(0)))
        sign = np.linalg.det(vt.T@u.T); vt[-1] *= sign
        R = vt.T@u.T
        initial = mesh.vertices[rep].copy()
        initial[free] = (rest.vertices[rep][free]-a.mean(0))@R.T+b.mean(0)
        target, proof = arap_core(rest.vertices[rep], edges, boundary,
                    mesh.vertices[rep][boundary], initial, maximum_steps=32, wall_seconds=.4)
        trial = mesh.vertices.copy(); trial[:] = target[inv]
        motion = float(np.linalg.norm(trial-mesh.vertices, axis=1).max())
        contexts = [bone_tree, skin_tree]
        for neighbour in data['candidate']['rows']:
            if neighbour['side'] == row['side'] and neighbour['partId'] != row['partId']:
                if sha(neighbour['path']) != neighbour['sha256']:
                    raise ValueError('Adjoining artery changed')
                contexts.append(BVH(Mesh.load(neighbour['path'])))
        candidate = Mesh(trial, mesh.faces, mesh.frame)
        summary, _ = self_contact_summary(candidate, example_limit=1000)
        contacts = []
        motion_bound = args.maximum_core_motion_mm*.001
        if motion <= motion_bound and summary['complete_within_pair_budget'] and summary['defect_pairs'] < defects:
            for context in contexts:
                before, after = [crossings(m, context, started) for m in [mesh, candidate]]
                contacts.append(dict(before=before, after=after))
                if after > before:
                    break
            if len(contacts) == len(contexts) and all(r['after'] <= r['before'] for r in contacts):
                q, defects = trial, summary['defect_pairs']
        strict_proof = dict(complete=defects==0, declaredSourceCore=args.declared_core_node,
                            sourceControls=len(free), originalSectionsAndMouthsLiteral=True,
                            sourceFacesRetained=True, maximumMotion_m=motion, motionBound_m=motion_bound,
                            before=first['defect_pairs'], proposed=summary['defect_pairs'], after=defects,
                            contextCounts=contacts, arap=proof)
    elif args.strict_section_uncross:
        from geometry_tools.algorithms.monotone_sections import repair_sections
        contexts=[('bone',bone_tree),('skin',skin_tree)]
        for neighbour in data['candidate']['rows']:
            if neighbour['side']==row['side']and neighbour['partId']!=row['partId']:
                if sha(neighbour['path'])!=neighbour['sha256']:raise ValueError('Adjoining artery changed')
                contexts.append((neighbour['partId'],BVH(Mesh.load(neighbour['path']))))
        protected=np.unique(inv[np.concatenate(graph['boundary_raw'])])
        q,strict_proof=repair_sections(mesh,rest,inv,graph['profiles'],graph['annularPropagation']['strips'],
                                      protected,contexts,wall_seconds=18.)
        defects=strict_proof['finalSelfPairs']
    elif args.strict_local_uncross:
        from geometry_tools.algorithms.monotone_uncross import monotone_uncross
        contexts=[bone_tree,skin_tree];context_names=['bone','skin']
        for neighbour in data['candidate']['rows']:
            if neighbour['side']==row['side']and neighbour['partId']!=row['partId']:
                if sha(neighbour['path'])!=neighbour['sha256']:raise ValueError('Adjoining artery changed')
                contexts.append(BVH(Mesh.load(neighbour['path'])))
                context_names.append(neighbour['partId'])
        groups=[]
        if args.rigid_section_groups:
            protected=np.unique(inv[np.concatenate(graph['boundary_raw'])])
            groups=[inv[s].tolist()for s in graph['profiles']]
        q, strict_proof = monotone_uncross(mesh, rest, inv, protected, contexts, wall_seconds=18.,
                                          preserve_context_counts=True,rigid_groups=groups)
        strict_proof['contextNames']=context_names
        defects = strict_proof['finalSelfPairs']
    for iteration in range(0 if args.strict_local_uncross or args.strict_section_uncross or args.strict_body_uncross or args.fixed_mouth_corridor or args.declared_core_node is not None else 4):
        current = Mesh(q, mesh.faces, mesh.frame)
        audit, _ = self_contact_summary(current, example_limit=1000)
        if not audit['defect_pairs']:
            break
        fids = np.unique(audit['examples'])
        active = np.unique(inv[mesh.faces[fids]])
        neighbourhood = np.unique(edges[np.isin(edges, active).any(1)])
        hit, _ = skin.closest_points(q[rep])
        tri = skin_mesh.triangles[hit['face_ids']]
        inward = -np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
        inward /= np.linalg.norm(inward, axis=1)[:, None]
        gap = np.einsum('ij,ij->i', q[rep]-hit['points'], inward)
        pinned = protected if args.slide_skin_contacts else np.unique(np.r_[protected, np.flatnonzero(gap < .00022)])
        free = np.setdiff1d(neighbourhood, pinned)
        boundary = np.setdiff1d(np.arange(len(rep)), free)
        target, proof = arap_core(rest.vertices[rep], edges, boundary, q[rep][boundary], q[rep],
                                  maximum_steps=32, wall_seconds=.4)
        delta = target-q[rep]
        maximum = np.linalg.norm(delta, axis=1).max(initial=0)
        delta *= min(1., .0002/max(maximum, 1e-15))
        trials = []
        accepted = False
        for strength in [1., .5, .25]:
            moved = strength*delta.copy()
            if args.slide_skin_contacts:
                bh, _ = bone_index.closest_points(q[rep])
                outward = q[rep]-bh['points']
                outward /= np.linalg.norm(outward, axis=1)[:, None]
                # Contact points may slide tangentially. Numerical halfspaces
                # only propose positions; full triangle contacts below decide.
                for _ in range(3):
                    for normal, distance in [(inward, gap), (outward, bh['distance_m'])]:
                        missing = np.maximum(.00021-distance[free]-np.einsum('ij,ij->i', normal[free], moved[free]), 0)
                        moved[free] += normal[free]*missing[:, None]
            trial = q+moved[inv]
            if np.linalg.norm(moved,axis=1).max(initial=0)>.0002+1e-12:
                trials.append(dict(strength=strength,complete=False,rejected='Core step exceeds 0.2 mm bound'))
                continue
            candidate = Mesh(trial, mesh.faces, mesh.frame)
            summary, _ = self_contact_summary(candidate, example_limit=1000)
            count = summary['defect_pairs']
            record = dict(strength=strength, selfPairs=count, complete=summary['complete_within_pair_budget'])
            if record['complete'] and count < audit['defect_pairs']:
                contacts = [crossings(candidate, tree, started) for tree in [bone_tree, skin_tree]]
                record['boneSkinCrossings'] = contacts
                if contacts == [0, 0]:
                    q, defects, accepted = trial, count, True
            trials.append(record)
            if accepted:
                break
        stages.append(dict(freeControls=len(free), arap=proof, trials=trials, accepted=accepted))
        if not accepted:
            break
        if time.monotonic()-started > 18:
            raise TimeoutError('18-second core repair deadline')
    normals = np.zeros((len(rep), 3))
    if not np.array_equal(q[np.isin(inv,protected)],mesh.vertices[np.isin(inv,protected)]):
        raise ValueError('A protected original profile or port moved during core repair')
    for section in graph['profiles']:
        if np.max(np.abs(np.linalg.norm(q[section][:,None]-q[section][None,:],axis=2)-
                         np.linalg.norm(mesh.vertices[section][:,None]-mesh.vertices[section][None,:],axis=2)))>1e-12:
            raise ValueError('A recorded source section lost its original calibre')
    t = q[mesh.faces]
    fn = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
    np.add.at(normals, inv[mesh.faces.ravel()], np.repeat(fn, 3, axis=0))
    normals = normals[inv]
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    attrs = dict(mesh.attributes, normals=np.rint(normals*32767).astype('<i2'), normalSigns=np.ones(len(q), 'i1'))
    result = Mesh(q, mesh.faces, mesh.frame, attrs, mesh.source)
    out.mkdir()
    file = out/'repaired-part.npz'
    result.save(file)
    report = copy.deepcopy(data['candidate'])
    r = next(r for r in report['rows'] if r['partId'] == args.part)
    r.update(path=str(file), sha256=sha(file), geometrySha256=result.geometry_sha)
    for pair in report.get('pairs',[]):
        if args.part in [pair['firstPart'],pair['secondPart']]:
            pa,pb=pair['firstPart'],pair['secondPart']
            other=pb if args.part==pa else pa
            other_row=next(r for r in report['rows'] if r['partId']==other)
            peer=Mesh.load(other_row['path'])
            own_raw=np.array(pair['firstRaw'] if args.part==pa else pair['secondRaw'],int)
            peer_raw=np.array(pair['secondRaw'] if args.part==pa else pair['firstRaw'],int)
            literal=bool(np.array_equal(result.vertices[own_raw],peer.vertices[peer_raw]))
            own_words=result.attributes['normals'][own_raw].astype(np.int32)*result.attributes['normalSigns'][own_raw,None]
            peer_words=peer.attributes['normals'][peer_raw].astype(np.int32)*peer.attributes['normalSigns'][peer_raw,None]
            normal_literal=bool(np.array_equal(own_words,peer_words))
            pair.update(literalCoordinatesAndEffectiveNormals=literal and normal_literal,coordinatesLiteral=literal,
                        requiresRecoupling=not literal,effectiveNormalWordsRequireSync=not normal_literal,
                        packedCoordinatesAndEffectiveNormalWordsLiteral=False)
    local_complete=defects==0 and(strict_proof is None or strict_proof['complete'])
    report.update(coreRepair=dict(initialSelfPairs=first['defect_pairs'], finalSelfPairs=defects, stages=stages,
                                  recordedSectionsAndPortsLiteral=not(args.strict_section_uncross or args.rigid_section_groups or args.strict_body_uncross),
                                  recordedSectionsCongruent=True,strictLocalRepair=strict_proof),
                  localRepairComplete=local_complete,complete=False, installed=False, installationAllowed=False,
                  inputs={k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()},
                  sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json', report)
    print(json.dumps(dict(before=first['defect_pairs'], after=defects,
                          targetBefore=strict_proof.get('initialTargetPairs') if strict_proof else None,
                          targetAfter=strict_proof.get('finalTargetPairs') if strict_proof else None,
                          localRepairComplete=local_complete, elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
