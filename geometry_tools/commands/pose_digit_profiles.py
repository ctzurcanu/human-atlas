"""Propose an articulated hand course using original rigid vascular sections.

No atlas installation. Fixed bone meshes are read-only. Original branch cores
and every source face remain owned; profile and boundary gauges are audited.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.course_io import sha, write_json
from geometry_tools.model import Mesh
from geometry_tools.algorithms.digit_course_carrier import DigitCourseCarrier
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores


def endpoints(centres, bones):
    a, z = centres[1]-centres[0], centres[-1]-centres[-2]
    a /= np.linalg.norm(a); z /= np.linalg.norm(z)
    first = centres[0]+a*np.min((bones[0]-centres[0])@a)
    last = centres[-1]+z*np.max((bones[-1]-centres[-1])@z)
    return np.r_[first[None], centres, last[None]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['guides', 'reference', 'profiles', 'context', 'output']:
        parser.add_argument('--'+key, required=True)
    parser.add_argument('--side', choices=['left', 'right'], required=True)
    parser.add_argument('--radial-bone-clearance', action='store_true')
    parser.add_argument('--bone-envelopes')
    parser.add_argument('--project-bone-clearance', action='store_true')
    parser.add_argument('--branch-frame-transport', action='store_true')
    parser.add_argument('--round-branch-courses', action='store_true')
    parser.add_argument('--component-ownership')
    parser.add_argument('--system', choices=['arterial', 'venous'])
    args = parser.parse_args()
    start = time.monotonic()
    paths = {k: Path(getattr(args, k)).resolve() for k in ['guides', 'reference', 'profiles', 'context']}
    guides, reference, bank, context = [json.loads(paths[k].read_text()) for k in paths]
    if not bank['complete'] or bank['contextSha256'] != sha(paths['context']):
        raise ValueError('Complete profiles for the actual context are required')
    if reference['contextSha256'] != sha(paths['context']):
        raise ValueError('Source bone frame/context mismatch')
    envelopes = None
    if args.bone_envelopes:
        ep = Path(args.bone_envelopes).resolve()
        envelopes = json.loads(ep.read_text())
        if envelopes['contextSha256'] != sha(paths['context']):
            raise ValueError('Canonical bone envelopes belong to another skeleton context')
        envelopes = {r['sourceObject']: r for r in envelopes['rows']}
        paths['bone_envelopes'] = ep
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new candidate stage')
    out.mkdir()
    courses, radial_surfaces = [], []
    for guide in guides['rows']:
        if guide['side'] != args.side:
            continue
        src, dst, sbones, tbones, smeshes, tmeshes = [], [], [], [], [], []
        for pair in guide['namedBoneCorrespondences']:
            row = next(r for r in reference['rows'] if r['targetOwner'] == pair['targetOwner'])
            if sha(row['sourcePath']) != row['sha256']:
                raise ValueError('Source bone changed')
            m = Mesh.load(row['sourcePath'])
            smeshes.append(m)
            sbones.append(m.vertices[np.unique(m.faces)])
            target = []
            for p in row['targetPieces']:
                if sha(p['path']) != p['sha256']:
                    raise ValueError('Fixed bone context changed')
                q = Mesh.load(p['path'])
                tmeshes.append(q)
                target.append(q.vertices[np.unique(q.faces)])
            tbones.append(np.concatenate(target))
            if envelopes is not None:
                envelope = envelopes[pair['targetOwner']]
                if sha(envelope['path']) != envelope['sha256']:
                    raise ValueError('Canonical fixed bone envelope changed')
                tmeshes = [m for m in tmeshes if m.source.get('canonical_owner') != pair['targetOwner']]
                # Replace this owner's material/cell pieces for radial rays;
                # those pieces remain the literal fixed display geometry.
                del tmeshes[-len(row['targetPieces']):]
                whole = Mesh.load(envelope['path'])
                whole.source['canonical_owner'] = pair['targetOwner']
                tmeshes.append(whole)
            src.append(pair['sourceCentre_m']); dst.append(pair['targetCentre_m'])
        courses.append(dict(digit=guide['digit'], source=endpoints(np.array(src), sbones).tolist(),
                            target=endpoints(np.array(dst), tbones).tolist()))
        if args.radial_bone_clearance:
            from geometry_tools.open3d_backend import SurfaceIndex
            def merged(items):
                tri = np.concatenate([m.triangles for m in items])
                return Mesh(tri.reshape(-1, 3), np.arange(tri.size//3).reshape(-1, 3), items[0].frame)
            radial_surfaces.append((SurfaceIndex(merged(smeshes)), SurfaceIndex(merged(tmeshes))))
    carrier = DigitCourseCarrier(courses, radial_surfaces=radial_surfaces if args.radial_bone_clearance else None)
    clearance_scene = None
    if args.project_bone_clearance:
        if envelopes is None:
            raise ValueError('Clearance projection requires verified whole-bone envelopes')
        from geometry_tools.open3d_backend import SurfaceIndex
        all_bones = [Mesh.load(r['path']) for owner, r in envelopes.items()
                     if owner.endswith('.l' if args.side == 'left' else '.r')]
        t = np.concatenate([m.triangles for m in all_bones])
        clearance_scene = SurfaceIndex(Mesh(t.reshape(-1, 3), np.arange(t.size//3).reshape(-1, 3), all_bones[0].frame))
    by_id = {r['partId']: r for r in context['blood']}
    ownership = None
    if args.component_ownership:
        op = Path(args.component_ownership).resolve()
        declaration = json.loads(op.read_text())
        if declaration['sourceBankSha256'] != sha(paths['profiles']):
            raise ValueError('Component binding belongs to a different source profile bank')
        if not all(r['anatomicallyReviewed'] for r in declaration['rows']):
            raise ValueError('Source component assignments require actual source review')
        ownership = declaration['rows']; paths['component_ownership'] = op
    rows = []
    for row in bank['rows']:
        if '-'+args.side+'_' not in row['partId']:
            continue
        if args.system and by_id[row['partId']]['sourceRecord']['system'] != args.system:
            continue
        if sha(row['path']) != row['sha256'] or sha(row['sourcePath']) != row['sourceSha256']:
            raise ValueError('Source profile bank changed')
        mesh = Mesh.load(row['sourcePath'])
        with np.load(row['path'], allow_pickle=False) as z:
            physical = z['physical_ids'].copy()
            original_profiles = [r.copy() for r in z['profiles']]
            original_links = z['profile_links'].copy() if 'profile_links' in z else None
        mouths = [np.array(m, int) for m in row['originalPhysicalMouths']]
        mouth_sets = [frozenset(physical[m]) for m in mouths]
        profiles, demoted, profile_mapping = [], [], {}
        for i, p in enumerate(original_profiles):
            own = frozenset(physical[p])
            if any(own & other and own != other for other in mouth_sets):
                demoted.append(i)
            else:
                profile_mapping[i] = len(profiles)
                profiles.append(p)
        owned = {frozenset(physical[p]) for p in profiles}
        for mouth, own in zip(mouths, mouth_sets):
            if own not in owned:
                profiles.append(mouth); owned.add(own)
        centres = np.array([mesh.vertices[p].mean(0) for p in profiles])
        support = np.array([np.linalg.norm(mesh.vertices[p]-c, axis=1).max() for p, c in zip(profiles, centres)])
        bindings = [r for r in ownership if r['partId'] == row['partId']] if ownership is not None else []
        if ownership is None:
            target, rotation, carrier_evidence = carrier.poses(centres, support)
        else:
            if not bindings:
                raise ValueError('No reviewed component declaration for '+row['partId'])
            covered_faces = [i for r in bindings for i in r['sourceFaceIds']]
            if sorted(covered_faces) != list(range(len(mesh.faces))):
                raise ValueError('Component declaration must own each original face exactly once')
            target = np.empty_like(centres); rotation = np.empty((len(centres), 3, 3))
            effective = np.zeros((len(centres), 5)); assigned = np.zeros(len(centres), bool); radial = []
            for binding in bindings:
                if binding['originalSha256'] != row['sourceSha256']:
                    raise ValueError('Component source geometry changed')
                raw = set(binding['rawVertexIds'])
                selected = np.array([i for i, p in enumerate(profiles) if int(p[0]) in raw], int)
                if any(not set(map(int, profiles[i])) <= raw for i in selected):
                    raise ValueError('A source profile spans distinct declared components')
                if assigned[selected].any():
                    raise ValueError('A profile has multiple component owners')
                numbers = binding['proposedCourses']
                if not numbers:
                    raise ValueError('Unresolved source component assignment')
                ci = [next(i for i, c in enumerate(courses) if c['digit'] == number) for number in numbers]
                component_carrier = DigitCourseCarrier([courses[i] for i in ci],
                    radial_surfaces=[radial_surfaces[i] for i in ci] if args.radial_bone_clearance else None)
                q, R, evidence = component_carrier.poses(centres[selected], support[selected])
                target[selected] = q; rotation[selected] = R; assigned[selected] = True
                effective[np.ix_(selected, np.array(numbers)-1)] = evidence['effective_weights']
                radial.append(dict(component=binding['physicalComponent'], courses=numbers,
                                   source='reviewed explicit component binding', proposals=evidence['radial_surface_proposals']))
            if not assigned.all():
                raise ValueError('Some source profiles lack an anatomical component owner')
            carrier_evidence = dict(effective_weights=effective, radial_surface_proposals=radial)
        clearance = None
        if clearance_scene is not None:
            from geometry_tools.algorithms.profile_clearance_projection import project_profile_centres
            movable = carrier_evidence['effective_weights'].sum(axis=1) > 0.
            target, clearance = project_profile_centres(clearance_scene, target, support, movable)
        frame_transport = None
        if args.branch_frame_transport:
            if original_links is None:
                raise ValueError('Branch transport needs verified actual source annular links')
            from geometry_tools.algorithms.profile_graph_transport import transport_linked_profiles, round_linked_courses
            links = np.array([[profile_mapping[a], profile_mapping[z]] for a, z in original_links
                              if a in profile_mapping and z in profile_mapping], int).reshape(-1, 2)
            protected = [i for i, p in enumerate(profiles)
                         if frozenset(physical[p]) in mouth_sets or
                         carrier_evidence['effective_weights'][i].sum() == 0.]
            rounded = None
            if args.round_branch_courses:
                if clearance_scene is None:
                    raise ValueError('Constrained rounding needs the actual fixed bone scene')
                target, rounded = round_linked_courses(centres, target, links, support, clearance_scene, protected,
                                                       reference_rotations=rotation)
            rotation, frame_transport = transport_linked_profiles(centres, target, rotation, links, protected)
            frame_transport['roundedCourses'] = rounded
        moved = transport_profile_cores(mesh.vertices, mesh.faces, profiles, target, rotation,
                                        physical_ids=physical)
        vertices = moved['vertices']
        maximum_gauge = 0.
        for ids in profiles:
            a, z = mesh.vertices[ids], vertices[ids]
            maximum_gauge = max(maximum_gauge, float(np.max(abs(
                np.linalg.norm(a[:, None]-a[None], axis=2)-
                np.linalg.norm(z[:, None]-z[None], axis=2)))))
        if maximum_gauge > 5e-8:
            raise ValueError('Actual source profile calibre changed')
        # Diagnostic smooth normals from the actual posed faces, retaining raw
        # source words through the bound source file. Packed install pending.
        face_normals = np.cross(vertices[mesh.faces[:, 1]]-vertices[mesh.faces[:, 0]],
                                vertices[mesh.faces[:, 2]]-vertices[mesh.faces[:, 0]])
        _, inv = np.unique(physical, return_inverse=True)
        normals = np.zeros((int(inv.max())+1, 3))
        np.add.at(normals, inv[mesh.faces.ravel()], np.repeat(face_normals, 3, axis=0))
        normals = normals[inv]
        length = np.linalg.norm(normals, axis=1)
        normals = np.divide(normals, length[:, None], out=np.zeros_like(normals), where=length[:, None] > 0)
        old = mesh.attributes.get('normals', np.zeros_like(normals)).astype(float)
        old *= mesh.attributes.get('normalSigns', np.ones(len(old)))[:, None]
        normals[length == 0] = old[length == 0]
        attrs = dict(mesh.attributes, normals=normals, normalSigns=np.ones(len(vertices)))
        candidate = Mesh(vertices, mesh.faces, mesh.frame, attrs, dict(sourcePartId=row['partId']))
        path = out/f'part-{len(rows):02}.npz'; candidate.save(path)
        rows.append(dict(partId=row['partId'], name=by_id[row['partId']]['name'],
                         system=by_id[row['partId']]['sourceRecord']['system'],
                         path=str(path), sha256=sha(path), sourcePath=row['sourcePath'],
                         sourceSha256=row['sourceSha256'], sourceGeometrySha256=mesh.geometry_sha,
                         geometrySha256=candidate.geometry_sha, rigidProfiles=len(profiles),
                         originalMouths=len(mouths), boundaryConflictingProfilesRetainedAsCore=demoted,
                         maximumProfileGaugeError_m=maximum_gauge,
                         maximumMotion_m=float(np.linalg.norm(vertices-mesh.vertices, axis=1).max()),
                         harmonicResidual_m=moved['harmonic_equation_maximum_residual_m'],
                         radialSurfaceProposals=carrier_evidence['radial_surface_proposals'],
                         centreClearanceProposal=clearance,
                         actualBranchFrameTransport=frame_transport,
                         unanchoredSourceComponents=moved['unanchored_components_retained_literal']))
        if time.monotonic()-start > 20:
            raise TimeoutError('20-second native articulated digit stage')
    report = dict(schema='articulated-digit-profile-candidate/v1', side=args.side, rows=rows,
                  inputs={k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()},
                  sourceCodeSha256=sha(__file__), carrierSha256=sha(Path(__file__).parents[1]/'algorithms/digit_course_carrier.py'),
                  courses=courses, skeletonChanged=False, sourceFacesRetained=True,
                  radialBoneClearanceProposed=args.radial_bone_clearance,
                  profileAndBoundaryGaugeRetained=True, complete=True, installed=False,
                  installationAllowed=False, anatomicalAcceptance=False,
                  pending=['Recipient surface/containment and discrete/self geometry checks.',
                           'Matched junction replay and proximal wrist/palmar connections.',
                           'Packed normal encoding and actual source/candidate rendering.'],
                  elapsedSeconds=time.monotonic()-start)
    write_json(out/'report.json', report)
    print(json.dumps(dict(side=args.side, parts=len(rows), rigidProfiles=sum(r['rigidProfiles'] for r in rows),
                          maximumGaugeError_m=max(r['maximumProfileGaugeError_m'] for r in rows),
                          elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
