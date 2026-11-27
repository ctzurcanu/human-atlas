"""Clear declared imported vessel components using their original section bank.

All controls, faces and component ownership come from recorded source data.
This produces a diagnostic candidate, never an installation authorization.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,complete_source_sections
from geometry_tools.algorithms.profile_core_transport import profile_translation_weights
from geometry_tools.commands.clear_hand_profiles import section_proposal
from geometry_tools.commands.check_hand_network import combine
from geometry_tools.algorithms.source_control_clearance import clear_source_controls
from geometry_tools.contact import self_contact_summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['candidate', 'profiles', 'binding', 'bones', 'skin-reference', 'part', 'output']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--maximum-motion-mm', type=float, default=3.)
    p.add_argument('--movable-source-cores', action='store_true')
    p.add_argument('--collision-report')
    p.add_argument('--recipient-context')
    p.add_argument('--recipient-face-guards',help='Hash-bound actual source-face and fixed-recipient separation planes; requires independent final triangle replay')
    p.add_argument('--components',help='Comma-separated complete original physical components; all other source arrays remain literal')
    p.add_argument('--retain-source-ports',action='store_true',help='Pin every recorded mouth of the selected components')
    p.add_argument('--solver-seconds',type=float,default=.4,help='Explicit optimizer cap up to one second; native source-component stage remains capped at two seconds')
    p.add_argument('--harmonic-source-cores',action='store_true',help='Retain every source core point using harmonic section influence, rather than independent core variables')
    p.add_argument('--core-patch-mm',type=float,default=0.,help='Optional source-edge-connected rigid core patches, maximum original diameter up to 3 mm; changes controls, never mesh detail')
    p.add_argument('--recipient-clearance-mm',type=float,default=.2,help='Explicit soft-recipient separation; bone and skin retain their independent 0.2 mm margin')
    p.add_argument('--elastic-recipient-proposal',action='store_true',help='Allow penalized recipient residuals during numerical route search only; final completion still requires every original clearance and actual triangle replay')
    p.add_argument('--rigid-translation',action='store_true')
    p.add_argument('--extend-source-annuli',action='store_true')
    p.add_argument('--capsules-outside-bone',action='store_true')
    p.add_argument('--warm-start',help='Saved complete source-control route proposal; retains source-reference face and fold guards')
    p.add_argument('--allow-source-curvature-change',action='store_true',help='Allow surface curvature sign changes while retaining source section calibre and face orientation; require exact intrinsic-crossing replay')
    args = p.parse_args()
    if not 0 < args.maximum_motion_mm <= 5:
        p.error('Motion bound must be within (0,5] mm')
    started = time.monotonic()
    paths = {k: Path(getattr(args, k.replace('-', '_'))).resolve()
             for k in ['candidate', 'profiles', 'binding', 'bones', 'skin-reference']}
    data = {k: json.loads(v.read_text()) for k, v in paths.items()}
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new clearance evidence stage')
    row = next(r for r in data['candidate']['rows'] if r['partId'] == args.part)
    old = next(r for r in data['profiles']['rows'] if r['partId'] == args.part)
    for r in [row, old]:
        if sha(r['path']) != r['sha256']:
            raise ValueError('Input mesh/profile bank changed')
    m = Mesh.load(row['path'])
    external_guards=[]
    if args.recipient_face_guards:
        path=Path(args.recipient_face_guards).resolve();saved=json.loads(path.read_text())
        if saved['sourceMeshSha256']!=row['sha256'] or saved['partId']!=args.part:
            raise ValueError('Source-face separation guards refer to a different proposal')
        paths['recipient-face-guards']=path
        for guard in saved['guards']:
            if sha(guard['recipientPath'])!=guard['recipientSha256']:
                raise ValueError('Guarded recipient mesh changed')
        external_guards=saved['guards']
    warm=None
    if args.warm_start:
        path=Path(args.warm_start).resolve();warm_report=json.loads(path.read_text());paths['warm-start']=path
        if not (warm_report.get('complete')or warm_report.get('finiteSourceRoutes')and all(p['complete']for p in warm_report['finiteSourceRoutes'])):
            raise ValueError('Warm-start algorithm did not complete its native stage')
        wr=next(r for r in warm_report['rows']if r['partId']==args.part)
        if sha(wr['path'])!=wr['sha256']:raise ValueError('Saved warm-start source changed')
        wm=Mesh.load(wr['path'])
        if not np.array_equal(wm.faces,m.faces):raise ValueError('Warm-start source face ancestry differs')
        warm=wm.vertices
    guard_pairs = []
    if args.collision_report:
        collision_path = Path(args.collision_report).resolve()
        collision = json.loads(collision_path.read_text())
        cr = next(r for r in collision['rows'] if r['partId'] == args.part)
        if sha(cr['path']) != cr['sha256']:
            raise ValueError('Collision example changed')
        witness_mesh = Mesh.load(cr['path'])
        if not np.array_equal(witness_mesh.faces, m.faces):
            raise ValueError('Collision face ancestry differs')
        defects, _ = self_contact_summary(witness_mesh, example_limit=1000)
        if not defects['complete_within_pair_budget'] or defects['defect_pairs'] != len(defects['examples']):
            raise ValueError('Collision guard inventory is incomplete')
        guard_pairs = defects['examples']
        paths['collision-report'] = collision_path
    original = Mesh.load(row['sourcePath'])
    if sha(row['sourcePath']) != row['sourceSha256'] or not np.array_equal(original.faces, m.faces):
        raise ValueError('Source face ancestry changed')
    with np.load(old['path'], allow_pickle=False) as z:
        pairs = z['retained_source_triangle_pairs'].copy()
    graph = source_profile_graph(original.vertices, original.faces, triangle_pairing=pairs,
                                 extend_annuli=args.extend_source_annuli or data['candidate'].get('sourceAnnularPropagation', False))
    physical = np.unique(original.vertices, axis=0, return_inverse=True)[1]
    profiles = complete_source_sections(graph,physical)
    n = len(profiles)
    section_owner={int(physical[raw]):i for i,s in enumerate(profiles)for raw in s}
    links=set()
    for tri in original.faces:
        owners={section_owner.get(int(physical[raw]))for raw in tri}-{None}
        for a in owners:
            for b in owners:
                if a<b:links.add((a,b))
    links=np.array(sorted(links),int).reshape(-1,2)
    groups = [r for r in data['binding']['rows'] if r['partId'] == args.part]
    if not groups or not all(r['anatomicallyReviewed'] for r in groups):
        raise ValueError('Component ownership needs recorded review')
    component_scope=None if args.components is None else set(map(int,args.components.split(',')))
    if component_scope is not None:
        if not component_scope or not component_scope.issubset({r['physicalComponent']for r in groups}):
            raise ValueError('Unknown complete source component in requested scope')
        groups=[r for r in groups if r['physicalComponent']in component_scope]
    port_raw=[]
    if args.retain_source_ports:
        for pair in data['candidate']['pairs']:
            if pair['firstPart']==args.part:port_raw.extend(pair['firstRaw'])
            if pair['secondPart']==args.part:port_raw.extend(pair['secondRaw'])
    port_raw=np.unique(np.asarray(port_raw,int))
    fixed = [r for r in data['bones']['rows'] if r['partId'].endswith('.l' if row['side'] == 'left' else '.r')]
    skin = data['skin-reference']
    for r in [*fixed, skin]:
        if sha(r['path']) != r['sha256']:
            raise ValueError('Fixed context changed')
    if not skin['defaultHiddenExcluded']:
        raise ValueError('Clearance requires the actual active skin')
    bone_index = SurfaceIndex(combine([Mesh.load(r['path']) for r in fixed]))
    skin_index = SurfaceIndex(Mesh.load(skin['path']))
    recipient=[]
    if args.recipient_context:
        if not args.movable_source_cores:
            raise ValueError('Recipient surfaces require the explicit source-control solve')
        path=Path(args.recipient_context).resolve();context=json.loads(path.read_text())
        paths['recipient-context']=path
        for fixed_surface in context['rows']:
            if sha(fixed_surface['path'])!=fixed_surface['sha256']:
                raise ValueError('Fixed recipient surface changed')
            if 'querySource'in fixed_surface:
                evidence=context['queryEnvelopeEvidence'];proof_path=Path(evidence['path'])
                if sha(proof_path)!=evidence['sha256']:raise ValueError('Numerical envelope qualification changed')
                qualification=json.loads(proof_path.read_text())
                query_source=fixed_surface['querySource']
                if (sha(query_source['path'])!=query_source['sha256'] or
                    not qualification['originalFacesLiteral'] or not qualification['queryNeverRenderedOrInstalled'] or
                    qualification['qualifiedArterialMesh']['sha256']!=row['sha256'] or
                    args.maximum_motion_mm*.001>qualification['qualifiedMaximumMotion_m'] or
                    qualification['minimumArterialSampleToQueryCap_m']<=args.maximum_motion_mm*.001+args.recipient_clearance_mm*.001):
                    raise ValueError('Numerical end closures are not remote from this bounded source proposal')
            policy=fixed_surface.get('clearanceSidePolicy','initial-control-side')
            if policy in ['oriented-exterior','reference-positive']:
                if 'orientationWitnesses'not in fixed_surface or not fixed_surface['orientationWitnesses']['sourceWindingConsistent']:
                    raise ValueError('Explicit source surface orientation evidence required')
                if m.vertices[:,1].max()+args.maximum_motion_mm*.001>fixed_surface['qualifiedRegionMaxY_m']:
                    raise ValueError('Source corridor exceeds the qualified distal recipient region')
            if args.capsules_outside_bone and 'capsul' in fixed_surface['partId'].lower() and policy!='oriented-exterior':
                raise ValueError('Capsule exterior requires checked component orientation; the nearest-bone direction can reverse within one capsule face')
            recipient.append((SurfaceIndex(Mesh.load(fixed_surface['path'])),fixed_surface['closedOrientedManifold'],policy,fixed_surface['partId'],fixed_surface.get('escapeDirection')))
    if args.rigid_translation and not args.movable_source_cores:
        raise ValueError('Rigid translation requires the explicit source-control solve')
    physical = np.unique(original.vertices, axis=0, return_inverse=True)[1]
    q = m.vertices.copy()
    covered = np.zeros(len(m.faces), bool)
    evidence = []
    for group in groups:
        fids = np.array(group['sourceFaceIds'], int)
        if covered[fids].any():
            raise ValueError('Component face ownership overlaps')
        covered[fids] = True
        ids = np.unique(m.faces[fids])
        remap = np.full(len(q), -1, int)
        remap[ids] = np.arange(len(ids))
        local_faces = remap[m.faces[fids]]
        chosen = [i for i, section in enumerate(profiles) if np.all(remap[section] >= 0)]
        if not chosen:
            raise ValueError('Component lacks actual source section controls')
        rank = {i: j for j, i in enumerate(chosen)}
        local_profiles = [remap[profiles[i]] for i in chosen]
        local_links = np.array([[rank[a], rank[b]] for a, b in links
                                if a in rank and b in rank], int).reshape(-1, 2)
        if args.movable_source_cores:
            face_rank = {int(k): i for i, k in enumerate(fids)}
            local_guards = [(face_rank[a], face_rank[b]) for a, b in guard_pairs
                            if a in face_rank and b in face_rank]
            moved, proof = clear_source_controls(m.vertices[ids], local_faces, local_profiles,
                                                physical[ids], bone_index, skin_index,
                                                maximum_motion_m=args.maximum_motion_mm*.001,
                                                guard_face_pairs=local_guards,
                                                recipient_surfaces=recipient,
                                                rigid_translation=args.rigid_translation,
                                                initial_positions=None if warm is None else warm[ids],
                                                protected_raw=remap[port_raw][remap[port_raw]>=0],
                                                solver_seconds=args.solver_seconds,
                                                harmonic_cores=args.harmonic_source_cores,
                                                core_patch_diameter_m=args.core_patch_mm*.001,
                                                recipient_clearance_m=args.recipient_clearance_mm*.001,
                                                elastic_recipient_proposal=args.elastic_recipient_proposal,
                                                recipient_plane_guards=[dict(g,sourceRaw=remap[np.unique(m.faces[g['sourceFaceIds']])].tolist())
                                                    for g in external_guards if np.isin(g['sourceFaceIds'],fids).all()],
                                                preserve_dihedral_sign=not args.allow_source_curvature_change)
            if proof['complete']and args.allow_source_curvature_change:
                intrinsic,_=self_contact_summary(Mesh(moved,local_faces,m.frame))
                proof['actualIntrinsicReplay']=dict(pairs=intrinsic['defect_pairs'],complete=intrinsic['complete_within_pair_budget'])
                if intrinsic['defect_pairs']or not intrinsic['complete_within_pair_budget']:
                    proof.update(complete=False,status='IntrinsicCrossingUnresolved')
        else:
            target, proof = section_proposal(m.vertices[ids], local_faces, local_profiles, local_links,
                                            bone_index, skin_index, trust_m=args.maximum_motion_mm*.001,
                                            physical_ids=physical[ids], whole_source=True)
        proof['physicalComponent'] = group['physicalComponent']
        evidence.append(proof)
        if proof['complete']:
            if args.movable_source_cores:
                q[ids] = moved
            else:
                centres = np.array([m.vertices[ids][s].mean(0) for s in local_profiles])
                W = profile_translation_weights(m.vertices[ids], local_faces, local_profiles,
                                                physical_ids=physical[ids])
                q[ids] += W@(target-centres)
        if time.monotonic()-started > 22:
            raise TimeoutError('22-second native section-clearance stage')
    if component_scope is None and not covered.all():
        raise ValueError('Declared components do not cover the whole imported part')
    changed_raw=np.unique(m.faces[covered])
    untouched_raw=np.setdiff1d(np.arange(len(q)),changed_raw)
    if not np.array_equal(q[untouched_raw],m.vertices[untouched_raw]):
        raise ValueError('Unselected source component moved')
    gauge = max(float(np.abs(np.linalg.norm(q[s][:, None]-q[s][None, :], axis=2)-
                             np.linalg.norm(m.vertices[s][:, None]-m.vertices[s][None, :], axis=2)).max())
                for s in profiles)
    if gauge > 1e-12:
        raise ValueError('Original section calibre changed')
    _, inv = np.unique(physical, return_inverse=True)
    normals = np.zeros((inv.max()+1, 3))
    t = q[m.faces]
    fn = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
    np.add.at(normals, inv[m.faces.ravel()], np.repeat(fn, 3, axis=0))
    normals = normals[inv]
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    attrs={k:v.copy()for k,v in m.attributes.items()}
    attrs['normals'][changed_raw]=np.rint(normals[changed_raw]*32767).astype('<i2')
    attrs['normalSigns'][changed_raw]=1
    if len(port_raw):
        attrs['normals'][port_raw]=m.attributes['normals'][port_raw]
        attrs['normalSigns'][port_raw]=m.attributes['normalSigns'][port_raw]
    candidate = Mesh(q, m.faces, m.frame, attrs, m.source)
    out.mkdir()
    file = out/'cleared-part.npz'
    candidate.save(file)
    report = copy.deepcopy(data['candidate'])
    r = next(r for r in report['rows'] if r['partId'] == args.part)
    r.update(path=str(file), sha256=sha(file), geometrySha256=candidate.geometry_sha,
             parentPosePath=row['path'], parentPoseSha256=row['sha256'])
    for pair in report['pairs']:
        if args.part in [pair['firstPart'], pair['secondPart']]:
            pair.update(literalCoordinatesAndEffectiveNormals=False,
                        packedCoordinatesAndEffectiveNormalWordsLiteral=False, requiresRecoupling=True)
    report.update(sectionClearance=evidence, complete=all(x['complete'] for x in evidence),
                  sourceSectionGaugeError_m=gauge, sourceFacesLiteral=True,
                  skeletonChanged=False, portsRequireRecoupling=True,
                  installed=False, installationAllowed=False, sourceCodeSha256=sha(__file__),
                  sourceAnnularPropagation=args.extend_source_annuli or data['candidate'].get('sourceAnnularPropagation',False),
                  componentScope=None if component_scope is None else sorted(component_scope),
                  unselectedSourceArraysLiteral=True,sourceMouthsProtected=args.retain_source_ports,
                  inputs={k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()},
                  elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json', report)
    print(json.dumps(dict(complete=report['complete'], components=len(evidence),
                          proposed=sum(x['complete'] for x in evidence),
                          sections=n, gaugeError_m=gauge, elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
