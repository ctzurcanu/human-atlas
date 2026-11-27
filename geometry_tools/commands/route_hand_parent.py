"""Route one explicit two-mouth source artery between literal assembly mouths.

This stage preserves every source face and section gauge. It does not accept
anatomy: actual triangles, source interfaces, normals and rendered sections
must be reviewed independently before installing any connected assembly.
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
from geometry_tools.source_profile_io import load_source_pairing
from geometry_tools.commands.route_hand_axis import fixed_context_indices
from geometry_tools.commands.couple_hand_arteries import rotation
from geometry_tools.algorithms.source_profile_axis import source_profile_graph
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores
from geometry_tools.algorithms.layered_tube_corridor import layered_tube_corridor,route_station_correspondence
from geometry_tools.algorithms.section_configuration_field import SectionConfigurationField
from geometry_tools.algorithms.feature_normals import indexed_feature_normals, inward_sheet_clearance
from geometry_tools.algorithms.tube_sections import arclength, rotation_minimizing_frames,path_tangents
from geometry_tools.open3d_backend import library


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ['source', 'candidate', 'profiles', 'part', 'root-peer', 'bones',
              'skin-reference', 'recipient-context', 'output']:
        p.add_argument('--'+k, required=True)
    p.add_argument('--allow-lateral-detours', action='store_true',
                   help='Explicit reversible clearance graph for a previously diagnosed station-DAG barrier')
    p.add_argument('--free-root',action='store_true',
                   help='Permit the wrist attachment to move with an explicit upstream radial-artery recoupling obligation')
    p.add_argument('--guide-from-current-course',action='store_true')
    p.add_argument('--actual-section-outlines',action='store_true',
                   help='Screen unchanged non-circular source ring outlines instead of their bounding spheres')
    p.add_argument('--local-cut-report',help='Diagnostic local refinement of a saved disconnected lattice cut; emits no anatomy')
    p.add_argument('--cut-grid-mm',type=float,help='Explicit prior lattice spacing when historical cut evidence omits it')
    p.add_argument('--cut-half-width-mm',type=float,help='Explicit prior lattice width when historical cut evidence omits it')
    p.add_argument('--grid-mm',type=float,default=.5)
    p.add_argument('--half-width-mm',type=float,default=6.)
    args = p.parse_args(); start = time.monotonic()
    paths = {k: Path(getattr(args, k.replace('-', '_'))).resolve() for k in
             ['source', 'candidate', 'profiles', 'bones', 'skin-reference', 'recipient-context']}
    data = {k: json.loads(v.read_text()) for k, v in paths.items()}
    cut=None
    if args.local_cut_report:
        if not args.actual_section_outlines:raise ValueError('Local refinement requires the actual outline configuration screen')
        paths['local-cut-report']=Path(args.local_cut_report).resolve()
        cut=json.loads(paths['local-cut-report'].read_text())
        for k in paths.keys()-{'local-cut-report'}:
            if cut['inputs'][k]['sha256']!=sha(paths[k]):raise ValueError('Local cut belongs to different source/context inputs')
    out = Path(args.output).resolve()
    if out.exists(): raise ValueError('Use a new bounded parent-route stage')
    src = {r['partId']: r for r in data['source']['rows']}
    rows = {r['partId']: r for r in data['candidate']['rows']}
    selected = rows[args.part]; source = src[args.part]
    pairs = [r for r in data['candidate']['pairs'] if args.part in [r['firstPart'], r['secondPart']]]
    if len(pairs) != 2: raise ValueError('Exactly two explicitly declared source mouths required')
    mouths = []
    for pair in pairs:
        key, other = ('first', 'second') if pair['firstPart'] == args.part else ('second', 'first')
        peer = pair[other+'Part']; row = rows[peer]
        if sha(row['path']) != row['sha256']: raise ValueError('Assembly peer changed')
        m = Mesh.load(row['path']); raw = np.array(pair[key+'Raw'], int)
        target = m.vertices[np.array(pair[other+'Raw'], int)]
        mouths.append((peer, raw, target))
    mouths.sort(key=lambda r: r[0] != args.root_peer)
    if mouths[0][0] != args.root_peer: raise ValueError('Root peer must own an actual declared source mouth')
    if sha(source['path']) != source['sha256']: raise ValueError('Original source changed')
    old = Mesh.load(source['path'])
    if not np.array_equal(old.faces, Mesh.load(selected['path']).faces): raise ValueError('Source face ancestry changed')
    bank = next(r for r in data['profiles']['rows'] if r['partId'] == args.part)
    pairing=load_source_pairing(bank,old.faces)
    g = source_profile_graph(old.vertices, old.faces, triangle_pairing=pairing, extend_annuli=True)
    # A branch must be fitted as a network; this two-port operation cannot
    # discard an unmatched side branch or choose an arbitrary primary path.
    if len(g['boundary_nodes']) != 2 or np.max(np.diff(g['graph'].indptr)) > 2:
        raise ValueError('The declared parent is branched; use the whole assembly')
    nodes = []
    for peer, raw, target in mouths:
        matches = [n for n, ids in zip(g['boundary_nodes'], g['boundary_raw'])
                   if set(ids) == set(raw)]
        if len(matches) != 1: raise ValueError('Declared mouth does not match one source boundary')
        nodes.append(matches[0])
    _, pred = shortest_path(g['graph'], directed=False, indices=nodes[0], return_predecessors=True)
    order = [nodes[1]]
    while order[-1] != nodes[0]:
        k = int(pred[order[-1]])
        if k < 0 or k in order: raise ValueError('Invalid source section course')
        order.append(k)
    order = np.array(order[::-1]); axis = g['centres'][order]
    if len(order) != len(g['centres']): raise ValueError('Source contains unvisited section detail')
    fractions = arclength(axis); fractions /= fractions[-1]
    target = axis.copy()
    # Preserve the source course as the route guide, with a smooth distributed
    # displacement to the explicit downstream mouth; actual context selects
    # the corridor. Neither endpoint nor its source radius may shrink.
    delta = [point.mean(0)-axis[j] for j, (_, raw, point) in zip([0, -1], mouths)]
    u = fractions**2*(3-2*fractions)
    target += (1-u[:, None])*delta[0]+u[:, None]*delta[1]
    if args.guide_from_current_course:
        current=Mesh.load(selected['path'])
        target=np.array([current.vertices[g['profiles'][i]].mean(0)for i in order])
        target[0]=mouths[0][2].mean(0);target[-1]=mouths[-1][2].mean(0)
    s = arclength(target); sample = np.linspace(0, s[-1], int(np.ceil(s[-1]/.001))+1)
    guide = PchipInterpolator(s, target)(sample)
    radii = np.array([np.linalg.norm(old.vertices[g['profiles'][i]]-g['centres'][i], axis=1).max() for i in order])
    radius = np.interp(sample/s[-1], fractions, radii)
    fixed = [r for r in data['bones']['rows'] if r['partId'].endswith('.l' if selected['side'] == 'left' else '.r')]
    for r in [*fixed, data['skin-reference'], *data['recipient-context']['rows']]:
        if sha(r['path']) != r['sha256']: raise ValueError('Fixed recipient changed')
    skin_record = data['skin-reference']
    bone, skin, context, scene = fixed_context_indices(tuple((r['path'], r['sha256']) for r in fixed),
        (skin_record['path'], skin_record['sha256']),
        tuple((r['path'], r['sha256'], r['closedOrientedManifold'], r['partId'])
              for r in data['recipient-context']['rows']), .00002)
    def gap(q):
        result = []
        for first in range(0, len(q), 24000):
            v = q[first:first+24000]; hit, _ = bone.closest_points(v)
            tensor = library().core.Tensor(((v-bone.origin)/bone.scale).astype(np.float32))
            inside = bone.scene.compute_occupancy(tensor, nthreads=2, nsamples=3).numpy() > .5
            value = hit['distance_m']*np.where(inside, -1., 1.)-.0002
            hit, _ = skin.closest_points(v); normal, rim = indexed_feature_normals(skin, hit)
            _, clearance, _ = inward_sheet_clearance(v, hit, normal, rim)
            result.append(np.minimum(np.minimum(value, clearance-.0002), context.clearance(v)))
        return np.concatenate(result)
    outline_field=None
    if args.actual_section_outlines:
        old_frames=rotation_minimizing_frames(axis)
        root_R=rotation(old.vertices[mouths[0][1]],mouths[0][2])
        end_R=rotation(old.vertices[mouths[-1][1]],mouths[-1][2])
        guide_frames=rotation_minimizing_frames(guide,initial_reference=root_R@old_frames[0,:,0],
                                               final_reference=end_R@old_frames[-1,:,0])
        source_station=np.interp(sample,s,np.arange(len(order)))
        nearest=np.rint(source_station).astype(int)
        offsets=[]
        for k,j in enumerate(nearest):
            section=g['profiles'][int(order[j])]
            R=guide_frames[k]@old_frames[j].T
            offsets.append((old.vertices[section]-old.vertices[section].mean(0))@R.T)
        if len({len(x)for x in offsets})!=1:raise ValueError('Outline screen needs explicit equal-count original rings')
        offsets=np.array(offsets)
        if cut:
            from scipy.spatial import cKDTree
            boundary=np.array(cut['routeProof']['closestRegionPoints_m'])
            prior=np.load(paths['local-cut-report'].with_name('route.npz'))['guide']
            if not np.array_equal(prior,guide):raise ValueError('Saved cut guide changed')
            # Exact old lattice ancestry, including its transverse station.
            old_spacing=cut['routeProof'].get('gridSpacing_m')
            old_width=cut['routeProof'].get('halfWidth_m')
            if old_spacing is None:
                if args.cut_grid_mm is None:raise ValueError('Historical cut requires explicit prior grid spacing')
                old_spacing=args.cut_grid_mm*.001
            if old_width is None:
                if args.cut_half_width_mm is None:raise ValueError('Historical cut requires explicit prior lattice width')
                old_width=args.cut_half_width_mm*.001
            if not 0<old_spacing<=.001 or not old_spacing<=old_width<=.01:
                raise ValueError('Invalid explicit historical lattice parameters')
            n=int(np.floor(old_width/old_spacing));ij=np.indices((2*n+1,2*n+1)).reshape(2,-1).T-n
            F=rotation_minimizing_frames(guide)
            lattice=guide[:,None]+np.einsum('sij,kj->ski',F[:,:,:2],ij*old_spacing)
            residual,ids=cKDTree(lattice.reshape(-1,3)).query(boundary)
            if residual.max()>1e-12:raise ValueError('Cut boundary does not replay the actual saved source lattice')
            stations=ids//len(ij);count=max(5,int(np.ceil(np.linalg.norm(boundary[1]-boundary[0])/.00025))+1)
            local_indices=np.rint(np.linspace(*stations,count)).astype(int)
            offsets=offsets[local_indices];radius=radius[local_indices]
            guide=np.linspace(boundary[0],boundary[1],count)
            cut_proof=dict(sourceStations=stations.tolist(),sourceLatticeResidual_m=float(residual.max()),
                           priorGridSpacing_m=old_spacing,priorHalfWidth_m=old_width,
                           actualSectionOffsetsUnscaled=True,refinementOnly=True)
        outline_field=SectionConfigurationField(offsets,gap,wall_seconds=16.,ball_clearance=gap)
    route, proof = layered_tube_corridor(guide, radius, gap, segment_clear=scene.segment_clear,
        fixed_prefix=0 if args.free_root else 1, fixed_distal=True, grid_m=args.grid_mm*.001, half_width_m=args.half_width_mm*.001,
        transverse_neighbour_cells=2, wall_seconds=16.,allow_lateral_detours=args.allow_lateral_detours,
        configuration_clearance=outline_field)
    report = dict(complete=False, installed=False, installationAllowed=False, skeletonChanged=False,
        parentPart=args.part, fixedRootPeer=args.root_peer, sourceFacesRetained=True,
        routeProof=proof, endpointWitnesses=context.witnesses(guide[[0, -1]]),
        upstreamRootRequiresRecoupling=bool(args.free_root),
        endpointRadii_m=radius[[0, -1]].tolist(),
        actualSourceOutlineProposal=args.actual_section_outlines,
        outlineSurfaceQueries=0 if outline_field is None else outline_field.queries,
        outlineContainingBallsAccepted=0 if outline_field is None else outline_field.boundingBallAccepted,
        inputs={k:dict(path=str(v), sha256=sha(v)) for k, v in paths.items()}, sourceCodeSha256=sha(__file__))
    out.mkdir(); np.savez(out/'route.npz', guide=guide, radius=radius, route=np.empty((0,3)) if route is None else route)
    if cut:
        report.update(localCutRefinement=cut_proof,routeFound=route is not None,
                      geometryEmitted=False,wholeParentRouteComplete=False,
                      elapsedSeconds=time.monotonic()-start)
        write_json(out/'local-cut-probe.json',report)
        print(json.dumps(dict(routeProof=proof,localCutRefinement=cut_proof,geometryEmitted=False,
                              elapsedSeconds=report['elapsedSeconds'])))
        return
    if route is None:
        if 'closestRegionPoints_m' in proof:
            points=np.array(proof['closestRegionPoints_m'])
            report['barrierWitnesses']=context.witnesses(np.r_[points,points.mean(0)[None]])
            report['barrierGap_m']=gap(np.r_[points,points.mean(0)[None]]).tolist()
        report.update(elapsedSeconds=time.monotonic()-start)
        write_json(out/'route-failure.json', report)
        print(json.dumps(dict(routeProof={k:v for k,v in proof.items() if k not in ['partialRoute_m','nextFreeStates_m']},
                              endpointRadii_m=report['endpointRadii_m'],elapsedSeconds=report['elapsedSeconds']))); return
    # Transport every original profile by a proper rotation. The two original
    # mouth rings are prescribed literally from their assembly peers.
    if args.free_root:
        peer,raw,points=mouths[0];delta=route[0]-points.mean(0)
        mouths[0]=(peer,raw,points+delta)
        report['upstreamRootDisplacement_m']=delta.tolist()
    rs = arclength(route)
    station_map,station_proof=route_station_correspondence(guide,route,grid_m=args.grid_mm*.001,
                                                       half_width_m=args.half_width_mm*.001)
    posed = PchipInterpolator(rs,route)(np.interp(s,sample,station_map))
    report['sourceStationCorrespondence']=station_proof
    old_t=path_tangents(axis);new_t=path_tangents(posed)
    root_R=rotation(old.vertices[mouths[0][1]],mouths[0][2])
    end_R=rotation(old.vertices[mouths[-1][1]],mouths[-1][2])
    new_t[0]=root_R@old_t[0];new_t[-1]=end_R@old_t[-1]
    before = rotation_minimizing_frames(axis,old_t)
    after = rotation_minimizing_frames(posed,new_t,initial_reference=root_R@before[0,:,0],
                                      final_reference=end_R@before[-1,:,0])
    R = after@np.swapaxes(before, 1, 2)
    for j, (_, raw, point) in zip([0, -1], mouths): R[j] = rotation(old.vertices[raw], point)
    interior = [int(i) for i in order if i not in nodes]
    where = {int(k):i for i,k in enumerate(order)}
    raw = np.concatenate([r for _, r, _ in mouths]); values = np.concatenate([v for _, _, v in mouths])
    transported = transport_profile_cores(old.vertices, old.faces, [g['profiles'][i] for i in interior],
        np.array([posed[where[i]] for i in interior]), np.array([R[where[i]] for i in interior]),
        prescribed_vertices=raw, prescribed_positions=values)
    q = transported['vertices']; _, inv = np.unique(old.vertices, axis=0, return_inverse=True)
    for i, point in zip(raw, values): q[inv == inv[i]] = point
    n = np.zeros((inv.max()+1, 3)); tri = q[old.faces]
    fn = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0]); np.add.at(n, inv[old.faces.ravel()], np.repeat(fn,3,axis=0))
    n = n[inv]; n /= np.linalg.norm(n,axis=1)[:,None]
    attrs = dict(old.attributes, normals=np.rint(n*32767).astype('<i2'), normalSigns=np.ones(len(q),'i1'))
    final = Mesh(q, old.faces, old.frame, attrs, old.source); file = out/'parent.npz'; final.save(file)
    output = copy.deepcopy(data['candidate']); row = next(r for r in output['rows'] if r['partId'] == args.part)
    row.update(path=str(file), sha256=sha(file), geometrySha256=final.geometry_sha,
               sourcePath=source['path'], sourceSha256=source['sha256'])
    output.update(parentRoute=report, complete=False, installed=False, installationAllowed=False,
        portsRequireRecoupling=True, pending=['Packed interface normals, full triangle/containment checks and rendered review'],
        inputs=report['inputs'], sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-start)
    if time.monotonic()-start > 22: raise TimeoutError('22-second native explicit parent-route stage')
    write_json(out/'report.json', output)
    print(json.dumps(dict(route=proof, parentSourceFaces=len(old.faces), sourceSections=len(order),
                         sourceMouthCoordinatesLiteral=not args.free_root,
                         upstreamRootRequiresRecoupling=bool(args.free_root),elapsedSeconds=output['elapsedSeconds'])))


if __name__ == '__main__': main()
