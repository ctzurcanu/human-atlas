"""Propagate explicit digital-mouth poses into the existing palm parents.

Only the ninety declared source interfaces may share controls. All source
faces and branches remain; unpaired forearm boundaries retain their pose.
This is a construction stage requiring independent contact/render review.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.source_profile_io import load_source_pairing
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,complete_source_sections
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.commands.couple_hand_arteries import rotation


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'candidate', 'profiles', 'priority-part', 'output']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--section-graph', action='store_true')
    p.add_argument('--open3d-arap', action='store_true',
                   help='Use bounded native surface ARAP and rigid source-profile restoration instead of harmonic interpolation')
    p.add_argument('--fixed-parts',default='',
                   help='Comma-separated explicit interface partners whose coordinates and attributes must remain literal')
    p.add_argument('--geometry-fixed-parts',default='',
                   help='Keep coordinates literal while permitting only declared seam normal-word synchronization')
    p.add_argument('--protected-source-controls',help='Hash-bound actual raw controls that must keep their current coordinates')
    p.add_argument('--prescribed-source-controls',help='Hash-bound explicit source raw controls and targetPositions_m for a coupled assembly')
    p.add_argument('--tangent-parts',default='',
                   help='Explicit child parts whose inlet directions prescribe one existing annular parent collar')
    args = p.parse_args()
    started = time.monotonic()
    paths = {k: Path(getattr(args, k.replace('-', '_'))).resolve() for k in ['source', 'candidate', 'profiles']}
    data = {k: json.loads(v.read_text()) for k, v in paths.items()}
    protected_rows=[]
    if args.protected_source_controls:
        path=Path(args.protected_source_controls).resolve();paths['protected-source-controls']=path
        protected_rows=json.loads(path.read_text())['rows']
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new coupling stage')
    source_rows = {r['partId']: r for r in data['source']['rows']}
    rows = {r['partId']: r for r in data['candidate']['rows']}
    fixed=set(filter(None,args.fixed_parts.split(',')))
    geometry_fixed=set(filter(None,args.geometry_fixed_parts.split(',')))
    immovable=fixed|geometry_fixed
    if not immovable.issubset(rows):raise ValueError('Fixed interface partner missing from the assembly')
    meshes, original, physical, representative, prescribed = {}, {}, {}, {}, {}
    for pid, row in rows.items():
        sr = source_rows[pid]
        for r in [row, sr]:
            if sha(r['path']) != r['sha256']:
                raise ValueError('Source binding changed')
        m, s = Mesh.load(row['path']), Mesh.load(sr['path'])
        if not np.array_equal(m.faces, s.faces):
            raise ValueError('Coupling requires literal source face ancestry')
        meshes[pid], original[pid] = m, s
        _, rep, inv = np.unique(s.vertices, axis=0, return_index=True, return_inverse=True)
        physical[pid], representative[pid] = inv, rep
        prescribed[pid] = {int(k): m.vertices[rep[k]].copy()
                           for loop in boundary_loops(inv[s.faces]) for k in loop}
    pairs = copy.deepcopy(data['source']['pairs'])
    local_fixed={pid:set()for pid in rows}
    for item in protected_rows:
        pid=item['partId']
        if pid not in rows or item['sourceMeshSha256']!=rows[pid]['sha256']:
            raise ValueError('Protected source controls belong to another actual input mesh')
        raw=np.asarray(item['protectedRaw'],int)
        if raw.ndim!=1 or not len(raw)or raw.min()<0 or raw.max()>=len(meshes[pid].vertices):
            raise ValueError('Protected controls must be explicit actual source raw IDs')
        for k in physical[pid][raw]:
            k=int(k);local_fixed[pid].add(k)
            prescribed[pid][k]=meshes[pid].vertices[representative[pid][k]].copy()
    imposed={pid:set()for pid in rows}
    if args.prescribed_source_controls:
        path=Path(args.prescribed_source_controls).resolve();paths['prescribed-source-controls']=path
        spec=json.loads(path.read_text())
        for item in spec['rows']:
            pid=item['partId']
            if pid not in rows or item['sourceMeshSha256']!=rows[pid]['sha256']or spec['frame']!=meshes[pid].frame:
                raise ValueError('Prescribed controls must bind the actual candidate mesh and frame')
            raw=np.asarray(item['sourceRaw'],int);target=np.asarray(item['targetPositions_m'],float)
            if raw.ndim!=1 or not len(raw)or raw.min()<0 or raw.max()>=len(meshes[pid].vertices)or target.shape!=(len(raw),3)or not np.isfinite(target).all():
                raise ValueError('Finite target coordinates for explicit actual source controls required')
            for k,point in zip(physical[pid][raw],target):
                k=int(k)
                if k in local_fixed[pid]and not np.array_equal(prescribed[pid][k],point):
                    raise ValueError('Prescribed control conflicts with a protected source point')
                prescribed[pid][k]=point.copy();local_fixed[pid].add(k);imposed[pid].add(k)
    for pair in pairs:
        pa, pb = pair['firstPart'], pair['secondPart']
        a, b = np.array(pair['firstRaw']), np.array(pair['secondRaw'])
        if not np.array_equal(original[pa].vertices[a], original[pb].vertices[b]):
            raise ValueError('A source interface is not literal in the installed source bank')
        canonical = args.priority_part if args.priority_part in [pa, pb] else pair['canonicalPart']
        fixed_pair=immovable.intersection([pa,pb])
        if fixed_pair:canonical=args.priority_part if args.priority_part in fixed_pair else sorted(fixed_pair)[0]
        else:
            imposed_part=[pid for pid,raw in [(pa,a),(pb,b)]if all(int(k)in imposed[pid]for k in physical[pid][raw])]
            if len(imposed_part)==1:canonical=imposed_part[0]
        ids = a if canonical == pa else b
        target = np.array([prescribed[canonical].get(int(k),meshes[canonical].vertices[i])for i,k in zip(ids,physical[canonical][ids])])
        for pid, raw in [(pa, a), (pb, b)]:
            for k, point in zip(physical[pid][raw], target):
                if int(k)in local_fixed[pid]and not np.array_equal(prescribed[pid][int(k)],point):
                    raise ValueError('A shared mouth conflicts with a protected current source control')
                prescribed[pid][int(k)] = point.copy()
        pair['canonicalPart'] = canonical
    collars=[]
    tangent_parts=set(filter(None,args.tangent_parts.split(',')))
    if not tangent_parts.issubset(rows):raise ValueError('Tangent child is missing from the source assembly')
    if tangent_parts:
        from geometry_tools.algorithms.tube_collar_pose import pose_parent_collar
        def next_source_ring(pid,mouth):
            from geometry_tools.algorithms.source_tube_topology import _annulus
            s=original[pid];inv=physical[pid];mouth_keys=set(map(int,inv[mouth]))
            incident=np.isin(inv[s.faces],list(mouth_keys)).any(1)
            adjacent=set(map(int,inv[s.faces[incident]].ravel()))-mouth_keys
            # A recorded end strip can be non-planar and therefore absent
            # from the regular-profile graph. Its actual half-edge topology
            # must still certify two oriented rims and Euler characteristic 0.
            if len(adjacent)>=3 and _annulus(inv[s.faces[incident]],mouth_keys,adjacent):
                return representative[pid][np.array(sorted(adjacent),int)]
            bank=next(b for b in data['profiles']['rows']if b['partId']==pid)
            pairing=load_source_pairing(bank,s.faces)
            g=source_profile_graph(s.vertices,s.faces,triangle_pairing=pairing,extend_annuli=True)
            candidates=[cycle for cycle in g['profiles']if set(map(int,inv[cycle])).issubset(adjacent)]
            if len(candidates)!=1:
                raise ValueError(f'Expected one actual adjacent source annulus for {pid}, found {len(candidates)}')
            return candidates[0]
        for pair in pairs:
            hits=tangent_parts.intersection([pair['firstPart'],pair['secondPart']])
            if not hits:continue
            if len(hits)!=1:raise ValueError('A tangent interface needs one declared child')
            child=next(iter(hits));parent=pair['secondPart']if child==pair['firstPart']else pair['firstPart']
            ci=np.array(pair['firstRaw']if child==pair['firstPart']else pair['secondRaw'],int)
            pi=np.array(pair['firstRaw']if parent==pair['firstPart']else pair['secondRaw'],int)
            cn=next_source_ring(child,ci);pn=next_source_ring(parent,pi)
            posed,proof=pose_parent_collar(original[parent].vertices[pi],original[parent].vertices[pn],
                                           meshes[child].vertices[ci],meshes[child].vertices[cn])
            for k,point in zip(physical[parent][pn],posed):
                if int(k)in prescribed[parent]and not np.allclose(prescribed[parent][int(k)],point,rtol=0,atol=1e-12):
                    raise ValueError('A tangent collar conflicts with another fixed source interface')
                prescribed[parent][int(k)]=point.copy()
            collars.append(dict(parentPart=parent,childPart=child,sourceCollarRaw=pn.tolist(),**proof))
    result = {}
    moved_parts = []
    for pid, m in meshes.items():
        inv, rep = physical[pid], representative[pid]
        keys = np.array(sorted(prescribed[pid]), int)
        raw = rep[keys]
        target = np.array([prescribed[pid][int(k)] for k in keys])
        if np.array_equal(m.vertices[raw], target):
            result[pid] = m
            continue
        if pid in immovable:raise ValueError('A fixed source interface would need to move')
        bank = next((r for r in data['profiles']['rows'] if r['partId'] == pid), None)
        pairing = load_source_pairing(bank, m.faces)
        g = source_profile_graph(original[pid].vertices, m.faces, triangle_pairing=pairing, extend_annuli=True)
        profiles = [s for s in g['profiles'] if not np.isin(inv[s], keys).any()]
        graph_proof = None
        if args.open3d_arap:
            from geometry_tools.open3d_backend import library
            o3d = library()
            native = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(m.vertices[rep]),
                                               o3d.utility.Vector3iVector(inv[m.faces]))
            first = native.deform_as_rigid_as_possible(o3d.utility.IntVector(keys.astype(np.int32)),
                        o3d.utility.Vector3dVector(target), max_iter=80)
            proposed = np.asarray(first.vertices)[inv].copy()
            prescribed_profiles = {int(k): point.copy() for k, point in zip(keys,target)}
            for section in profiles:
                centre = m.vertices[section].mean(0)
                posed_centre = proposed[section].mean(0)
                R = rotation(m.vertices[section], proposed[section])
                desired = (m.vertices[section]-centre)@R.T+posed_centre
                for k, point in zip(inv[section],desired):prescribed_profiles[int(k)]=point
            ids = np.array(sorted(prescribed_profiles),int)
            values = np.array([prescribed_profiles[int(k)]for k in ids])
            second_native = native.deform_as_rigid_as_possible(o3d.utility.IntVector(ids.astype(np.int32)),
                        o3d.utility.Vector3dVector(values), max_iter=80)
            q = np.asarray(second_native.vertices)[inv].copy()
            for k,point in zip(ids,values):q[inv==k]=point
            for section in profiles:
                if np.max(np.abs(np.linalg.norm(q[section][:,None]-q[section][None],axis=2)-
                                 np.linalg.norm(m.vertices[section][:,None]-m.vertices[section][None],axis=2)))>1e-12:
                    raise ValueError('Native source-profile restoration lost calibre')
            graph_proof = dict(method='Open3D surface ARAP plus rigid original section constraints',
                               nativeIterationsPerPass=80, sourcePhysicalControls=len(rep),
                               fixedSourceMouthControls=len(keys), rigidOriginalProfiles=len(profiles))
        elif args.section_graph:
            from geometry_tools.algorithms.port_section_graph import pose_section_graph
            fixed_positions, fixed_rotations = {}, {}
            for node, loop in zip(g['boundary_nodes'], g['boundary_raw']):
                posed = np.array([prescribed[pid][int(inv[k])] for k in loop])
                fixed_positions[node] = posed.mean(0)
                fixed_rotations[node] = rotation(m.vertices[loop], posed)
            for node,loop in enumerate(g['profiles']):
                if all(int(inv[k])in prescribed[pid]for k in loop):
                    posed=np.array([prescribed[pid][int(inv[k])]for k in loop])
                    fixed_positions[node]=posed.mean(0)
                    fixed_rotations[node]=rotation(m.vertices[loop],posed)
            posed, frames, graph_proof = pose_section_graph(g, fixed_positions, fixed_rotations)
            selected = [i for i, s in enumerate(g['profiles']) if not np.isin(inv[s], keys).any()]
            centres, rotations = posed[selected], frames[selected]
        else:
            first = transport_profile_cores(m.vertices, m.faces, [], np.empty((0, 3)), np.empty((0, 3, 3)),
                                            physical_ids=inv, prescribed_vertices=raw, prescribed_positions=target)
            centres = np.array([first['vertices'][s].mean(0) for s in profiles]).reshape(-1, 3)
            rotations = np.array([rotation(m.vertices[s], first['vertices'][s]) for s in profiles]).reshape(-1, 3, 3)
        if not args.open3d_arap:
            second = transport_profile_cores(m.vertices, m.faces, profiles, centres, rotations,
                                             physical_ids=inv, prescribed_vertices=raw, prescribed_positions=target,
                                             core_method='arap' if args.section_graph else 'harmonic',
                                             core_wall_seconds=.4)
            q = second['vertices']
        for k, point in zip(keys, target):
            q[inv == k] = point
        normal = np.zeros((len(rep), 3))
        t = q[m.faces]
        fn = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
        np.add.at(normal, inv[m.faces.ravel()], np.repeat(fn, 3, axis=0))
        normal = normal[inv]
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        attrs = dict(m.attributes, normals=np.rint(normal*32767).astype('<i2'), normalSigns=np.ones(len(q), 'i1'))
        changed_faces=np.any(q[m.faces]!=m.vertices[m.faces],axis=(1,2))
        changed_star=np.unique(inv[m.faces[changed_faces]])
        untouched=~np.isin(inv,changed_star)
        attrs['normals'][untouched]=m.attributes['normals'][untouched]
        attrs['normalSigns'][untouched]=m.attributes['normalSigns'][untouched]
        result[pid] = Mesh(q, m.faces, m.frame, attrs, m.source)
        actual_sections=complete_source_sections(g,inv)
        gauge=max((float(np.abs(np.linalg.norm(q[s][:,None]-q[s][None,:],axis=2)-
                               np.linalg.norm(m.vertices[s][:,None]-m.vertices[s][None,:],axis=2)).max())
                   for s in actual_sections),default=0.)
        moved_parts.append(dict(partId=pid, interiorSections=len(profiles),
                                sectionGraph=graph_proof,
                                checkedSourceSections=len(actual_sections),
                                sourceProfileGaugeError_m=gauge,sourceProfileGaugePreserved=gauge<=1e-10,
                                maximumMotion_m=float(np.linalg.norm(q-m.vertices, axis=1).max())))
        if time.monotonic()-started > 22:
            raise TimeoutError('22-second native hand coupling stage')
    attributes = {pid: {k: v.copy() for k, v in m.attributes.items()} for pid, m in result.items()}
    for pair in pairs:
        pa, pb = pair['firstPart'], pair['secondPart']
        a, b = np.array(pair['firstRaw']), np.array(pair['secondRaw'])
        if not np.array_equal(result[pa].vertices[a], result[pb].vertices[b]):
            raise ValueError('Declared interface coordinates did not replay literally')
        canonical = pair['canonicalPart']
        raw = a if canonical == pa else b
        ca = attributes[canonical]
        words = ca['normals'][raw].astype(np.int32)*ca.get('normalSigns', np.ones(len(result[canonical].vertices), int))[raw, None]
        for pid, ids in [(pa, a), (pb, b)]:
            for key, word in zip(physical[pid][ids], words):
                mask = physical[pid] == key
                if pid in fixed:
                    actual=attributes[pid]['normals'][mask].astype(np.int32)*attributes[pid]['normalSigns'][mask,None]
                    if not np.all(actual==word):raise ValueError('A fixed partner normal field would change')
                    continue
                attributes[pid]['normals'][mask] = word
                attributes[pid]['normalSigns'][mask] = 1
        pair.update(literalCoordinatesAndEffectiveNormals=True,
                    packedCoordinatesAndEffectiveNormalWordsLiteral=False, requiresRecoupling=False)
    out.mkdir()
    emitted = []
    for pid, m in result.items():
        final = Mesh(m.vertices, m.faces, m.frame, attributes[pid], m.source)
        path = out/f'part-{len(emitted):02}.npz'
        final.save(path)
        emitted.append(dict(rows[pid], path=str(path), sha256=sha(path), geometrySha256=final.geometry_sha,
            historicalSourceComparison=dict(path=rows[pid].get('sourcePath'),sha256=rows[pid].get('sourceSha256')),
            sourcePath=source_rows[pid]['path'],sourceSha256=source_rows[pid]['sha256']))
    report = copy.deepcopy(data['candidate'])
    report.update(rows=emitted, pairs=pairs, palmParentCoupling=moved_parts,
                  constructionComplete=all(r['sourceProfileGaugePreserved']for r in moved_parts),
                  complete=False, portsRequireRecoupling=False,
                  installed=False, installationAllowed=False, skeletonChanged=False,
                  pending=['Whole-network self, body contacts, packed joins and rendered review'],
                  inputs={k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()},
                  sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-started)
    report['comparisonBaseline']=dict(path=str(paths['source']),sha256=sha(paths['source']),
        role='Explicit current source bank whose literal face ancestry was independently checked for every part')
    report['fixedInterfaceParts']=sorted(fixed)
    report['geometryFixedInterfaceParts']=sorted(geometry_fixed)
    report['protectedSourceControls']=[dict(partId=pid,physicalControls=len(ids))for pid,ids in local_fixed.items()if ids]
    report['prescribedSourceControlCounts']=[dict(partId=pid,physicalControls=len(ids))for pid,ids in imposed.items()if ids]
    report['tangentCollars']=collars
    write_json(out/'report.json', report)
    print(json.dumps(dict(pairedPorts=len(pairs), movedPalmParents=len(moved_parts),
                          sourceFacesRetained=True, elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
