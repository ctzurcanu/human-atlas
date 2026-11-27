"""Place imported palmar branches on matched reference courses, without copying meshes.

Source incidence and original quad pairings determine the imported sections.
Digit and radial/ulnar course matching use recorded source ownership and a
common phalangeal section plane. Original faces and branches remain present.
This is a bounded diagnostic proposal, never an installer.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.algorithms.source_profile_axis import source_profile_graph, primary_axis,complete_source_sections
from geometry_tools.algorithms.tube_axis_map import axis_transport
from geometry_tools.algorithms.tube_sections import arclength, rotation_minimizing_frames, unit


def plane_offset(axis, digit, bones, key):
    bone = bones[digit]
    centre = np.array(bone[-2][key])
    normal = np.array(bone[-1][key])-bone[1][key]
    normal /= np.linalg.norm(normal)
    radial = np.array(bones[1][0][key])-bones[5][0][key]
    radial -= normal*(normal@radial)
    radial /= np.linalg.norm(radial)
    a, delta = axis[:-1], np.diff(axis, axis=0)
    den = delta@normal
    t = np.divide(-(a-centre)@normal, den, out=np.zeros(len(a)), where=abs(den)>1e-12)
    live = (t>=0)&(t<=1)&(abs(den)>1e-12)
    if not live.any():
        raise ValueError('Course does not cross its declared phalangeal reference plane')
    hits = a[live]+t[live, None]*delta[live]
    point = hits[np.linalg.norm(hits-centre, axis=1).argmin()]
    return float((point-centre)@radial)


def smooth_samples(axis, count=192):
    s = arclength(axis)
    curve = PchipInterpolator(s/s[-1], axis)
    t = np.linspace(0, 1, count)
    return curve(t), unit(curve.derivative()(t))


def trim_course(axis, distance):
    """Advance a guide mouth by literal arclength, retaining the distal course."""
    s = arclength(axis)
    if distance < 0 or distance >= s[-1]:
        raise ValueError('Guide mouth trim exceeds its course')
    first = np.array([np.interp(distance, s, axis[:, k]) for k in range(3)])
    return np.r_[first[None], axis[s > distance+1e-12]]


def trim_to_matched_bone_station(axis, source_point, source_bones, target_bones):
    """Place a guide start at the corresponding named-bone station plane.

    The source mouth need not be the beginning of a longer reference artery.
    Preserve its location along the explicit paired bone chain instead of
    stretching the source part over the reference artery's proximal course.
    This determines a proposed mouth station, not a vascular attachment.
    """
    a,b=np.asarray(source_bones,float),np.asarray(target_bones,float)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or len(a)<2:
        raise ValueError('Corresponding named bone chains required')
    da=np.diff(a,axis=0);length2=np.sum(da*da,axis=1)
    if (length2<=1e-20).any():raise ValueError('Distinct bone stations required')
    u=np.clip(np.sum((source_point-a[:-1])*da,axis=1)/length2,0.,1.)
    nearest=a[:-1]+u[:,None]*da
    k=int(np.linalg.norm(nearest-source_point,axis=1).argmin())
    point=b[k]+u[k]*(b[k+1]-b[k]);normal=unit(b[k+1]-b[k])
    d=np.diff(axis,axis=0);den=d@normal
    t=np.divide((point-axis[:-1])@normal,den,out=np.zeros(len(d)),where=abs(den)>1e-12)
    valid=np.flatnonzero((t>=0)&(t<=1)&(abs(den)>1e-12))
    if not len(valid):raise ValueError('Reference route misses the matched bone station')
    hits=axis[valid]+t[valid,None]*d[valid]
    j=int(valid[np.linalg.norm(hits-point,axis=1).argmin()])
    station=float(arclength(axis)[j]+t[j]*np.linalg.norm(d[j]))
    return trim_course(axis,station),dict(sourceBoneSegment=k,segmentFraction=float(u[k]),
        targetStationPoint_m=point.tolist(),guideTrim_m=station)


def retain_proximal_course(old, new, old_t, new_t, fraction=.35):
    """Keep source fork position/tangent and reach the recipient course smoothly."""
    if old.shape != new.shape or not 0 < fraction <= 1:
        raise ValueError('Corresponding course samples and bounded blend fraction required')
    t = np.linspace(0, 1, len(old))
    u = np.clip(t/fraction, 0, 1)
    fade = u**3*(10-15*u+6*u*u)
    posed = old+fade[:, None]*(new-old)
    tangent = unit(np.gradient(posed, t, axis=0))
    tangent[0], tangent[-1] = old_t[0], new_t[-1]
    posed[0] = old[0]
    return posed, tangent


def separate_guide_mouths(matches, *, margin=.0002, maximum=.008, step=.0002):
    """Separate downstream daughter mouths on a zero-length reference fork.

    This sets proposed ports only. The original finite parent fork still needs
    to be transported and coupled; a guide bifurcation is never a replacement
    for that source geometry.
    """
    evidence = []
    for i, a in enumerate(matches):
        for b in matches[i+1:]:
            needed = a['mouthRadius_m']+b['mouthRadius_m']+margin
            if np.linalg.norm(a['target']['axis'][0]-b['target']['axis'][0]) >= needed:
                continue
            for distance in np.arange(0, maximum+step/2, step):
                aa, bb = [trim_course(r['target']['axis'], distance) for r in [a, b]]
                gap = float(np.linalg.norm(aa[0]-bb[0]))
                if gap >= needed:
                    break
            else:
                raise ValueError('Reference fork cannot fit unchanged source mouths within the trim limit')
            for r, axis in [(a, aa), (b, bb)]:
                r['target']['axis'] = axis
                r['target']['proximalTrim_m'] = float(distance)
            evidence.append(dict(sourceComponents=[a['source']['owner']['physicalComponent'], b['source']['owner']['physicalComponent']],
                                 guideAdvance_m=float(distance), mouthCentreDistance_m=gap,
                                 requiredDistance_m=needed, sourceCalibreChanged=False))
    return evidence


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'profiles', 'binding', 'courses', 'bones', 'output']:
        p.add_argument('--'+key, required=True)
    p.add_argument('--side', choices=['left', 'right'], required=True)
    p.add_argument('--retain-source-ports', action='store_true')
    p.add_argument('--pose-strength', type=float, default=1.)
    p.add_argument('--extend-source-annuli', action='store_true')
    p.add_argument('--part-kind', choices=['proper-palmar', 'dorsal-digital', 'proper-thumb'], default='proper-palmar')
    p.add_argument('--target-chain',type=int,
                   help='Explicit reviewed reference chain for the thumb source cohort')
    p.add_argument('--palm-chord-frame',action='store_true',
                   help='Use the thumb mouth-to-tip chord and paired palm planes instead of the artificial source-core zigzag')
    p.add_argument('--bone-radial-frame',action='store_true',
                   help='Match the measured thumb-axis radial side relative to the paired named bone chains')
    p.add_argument('--base-candidate', help='Keep unrelated members of this saved whole-hand assembly')
    p.add_argument('--source-axis-report',help='Hash-bound source-contour axis for one explicitly identified original component')
    p.add_argument('--describe-matches-only',action='store_true',help='Save source/recipient course correspondence only, without constructing any vessel mesh')
    p.add_argument('--target-axis-route',help='Completed hash-bound recipient corridor; its moved mouth must be recoupled before installation')
    p.add_argument('--only-source-components', help='Comma-separated original physical component IDs')
    p.add_argument('--source-bones', help='Verified donor hand bone surfaces for local frame orientation')
    p.add_argument('--recipient-bones', help='Fixed recipient bone envelopes for local frame orientation')
    p.add_argument('--transport-current-course', action='store_true',
                   help='Move only a finite proximal prefix of a saved course; preserve its downstream clearance corrections')
    args = p.parse_args()
    if not 0 < args.pose_strength <= 1:
        p.error('Pose strength must be within (0,1]')
    paths = {k: Path(getattr(args, k)).resolve() for k in ['source', 'profiles', 'binding', 'courses', 'bones']}
    data = {k: json.loads(v.read_text()) for k, v in paths.items()}
    base = None
    if args.base_candidate:
        paths['base-candidate'] = Path(args.base_candidate).resolve()
        base = json.loads(paths['base-candidate'].read_text())
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new course-fitting stage')
    started = time.monotonic()
    bone_indices = None
    if args.source_bones or args.recipient_bones:
        if not args.source_bones or not args.recipient_bones:
            raise ValueError('Both source and recipient bone references required')
        from geometry_tools.open3d_backend import SurfaceIndex
        from geometry_tools.commands.check_hand_network import combine
        from scipy.ndimage import uniform_filter1d
        indices = []
        for key, path, filename in [('source-bones', args.source_bones, 'sourcePath'),
                                     ('recipient-bones', args.recipient_bones, 'path')]:
            paths[key] = Path(path).resolve(); ref = json.loads(paths[key].read_text())
            selected = [r for r in ref['rows'] if r.get('side') == args.side or
                         (r.get('partId', '').endswith('.l' if args.side == 'left' else '.r'))]
            for row in selected:
                if sha(row[filename]) != row['sha256']:
                    raise ValueError('Bone reference changed')
            indices.append(SurfaceIndex(combine([Mesh.load(r[filename]) for r in selected])))
        bone_indices = indices

        def surface_frames(axis, tangent, index):
            hit, _ = index.closest_points(axis); tri = index.mesh.triangles[hit['face_ids']]
            n = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
            radial = axis-hit['points']; radial[np.einsum('ij,ij->i', radial, n) < 0] *= -1
            radial -= tangent*np.einsum('ij,ij->i', radial, tangent)[:, None]
            radial = unit(radial)
            radial = uniform_filter1d(radial, size=9, axis=0, mode='nearest')
            radial -= tangent*np.einsum('ij,ij->i', radial, tangent)[:, None]
            radial = unit(radial)
            return np.stack([radial, np.cross(tangent, radial), tangent], axis=2)
    marker = {'proper-palmar':'proper_palmar_','dorsal-digital':'dorsal_digital_',
              'proper-thumb':'proper_digital_arteries_of_thumb_'}[args.part_kind]
    if args.part_kind=='proper-thumb'and args.target_chain is None:
        raise ValueError('The thumb needs an explicitly matched reference chain')
    if args.palm_chord_frame and args.part_kind!='proper-thumb':
        raise ValueError('The measured chord-frame correction is scoped to the thumb')
    if args.bone_radial_frame and not args.palm_chord_frame:
        raise ValueError('Bone radial matching requires the explicit thumb chord frame')
    source_row = next(r for r in data['source']['rows'] if f'-{args.side}_{marker}' in r['partId'])
    m = Mesh.load(source_row['path'])
    if sha(source_row['path']) != source_row['sha256']:
        raise ValueError('Imported source changed')
    axis_overrides={}
    if args.source_axis_report:
        paths['source-axis-report']=Path(args.source_axis_report).resolve()
        inventory=json.loads(paths['source-axis-report'].read_text())
        if 'rows'in inventory:
            items=[]
            for entry in inventory['rows']:
                if entry['partId']!=source_row['partId']or not entry['complete']:continue
                if sha(entry['path'])!=entry['sha256']:raise ValueError('Source-axis inventory record changed')
                items.append(json.loads(Path(entry['path']).read_text()))
        else:items=[inventory]
        for item in items:
            if (item['partId']!=source_row['partId'] or item['sourceMeshSha256']!=source_row['sha256'] or
                sha(item['axisPath'])!=item['axisSha256'] or item['proof']['sourceGeometryModified']):
                raise ValueError('Source-contour axis does not bind this literal input component')
            axis_overrides[item['component']]=item
    legacy = next(r for r in data['profiles']['rows'] if r['partId'] == source_row['partId'])
    if sha(legacy['path']) != legacy['sha256']:
        raise ValueError('Recorded original triangle pairing changed')
    with np.load(legacy['path'], allow_pickle=False) as z:
        pairing = z['retained_source_triangle_pairs'].copy()
    graph = source_profile_graph(m.vertices, m.faces, triangle_pairing=pairing,
                                 extend_annuli=args.extend_source_annuli)
    owners = [r for r in data['binding']['rows'] if r['partId'] == source_row['partId']]
    if not all(r['anatomicallyReviewed'] for r in owners):
        raise ValueError('Source anatomical ownership was not reviewed')
    bones = {r['digit']: r['namedBoneCorrespondences'] for r in data['bones']['rows'] if r['side'] == args.side}
    tips = np.array([bones[d][-1]['targetCentre_m'] for d in range(1, 6)])
    native_name = 'Dorsal digital' if args.part_kind == 'dorsal-digital' else 'Proper palmar digital'
    native = next(r for r in data['courses']['rows'] if r['side'] == args.side and native_name in r['partId'])
    if sha(native['path']) != native['sha256']:
        raise ValueError('Reference course data changed')
    targets = []
    with np.load(native['path'], allow_pickle=False) as z:
        for i, chain in enumerate(native['chains']):
            if args.part_kind=='proper-thumb'and i!=args.target_chain:continue
            axis = z['centres'][chain]
            root = int(axis[:, 1].argmax())
            branches = [axis[root::-1], axis[root:]] if args.part_kind == 'dorsal-digital' else [axis]
            for branch_index, branch in enumerate(branches):
                if len(branch) < 3:
                    continue
                if branch[0, 1] < branch[-1, 1]:
                    branch = branch[::-1]
                digit = int(np.linalg.norm(tips-branch[-1], axis=1).argmin()+1)
                if (digit==1 if args.part_kind=='proper-thumb'else digit>1 or args.part_kind=='dorsal-digital'):
                    targets.append(dict(chain=i, branch=branch_index, digit=digit, axis=branch,
                                        offset=plane_offset(branch, digit, bones, 'targetCentre_m')))
    sources = []
    for owner in owners:
        digit = int(np.argmin(owner['tipDistanceByDigit_m'])+1)
        if digit not in owner['coveredDigits']:
            raise ValueError('Primary distal course contradicts its reviewed source cohort')
        axis, proof = primary_axis(graph, m.faces[owner['sourceFaceIds']])
        axis_override=axis_overrides.get(owner['physicalComponent'])
        if axis_override:
            with np.load(axis_override['axisPath'],allow_pickle=False)as z:
                if not np.array_equal(np.sort(z['source_face_ids']),np.sort(owner['sourceFaceIds'])):
                    raise ValueError('Source-contour axis owns a different original face cohort')
                refined=z['axis'].copy()
            if np.linalg.norm(refined[0]-axis[0])>1e-10 or np.linalg.norm(refined[-1]-axis[-1])>1e-10:
                raise ValueError('Source-contour axis changed its declared original endpoints')
            axis=refined
            proof=dict(proof,sourceContourRefinement=dict(path=str(paths['source-axis-report']),sha256=sha(paths['source-axis-report'])))
        sources.append(dict(owner=owner, digit=digit, axis=axis, proof=proof,
                            offset=plane_offset(axis, digit, bones, 'sourceCentre_m')))
    matches = []
    digits=[1]if args.part_kind=='proper-thumb'else range(1 if args.part_kind=='dorsal-digital'else 2,6)
    for digit in digits:
        ss = sorted([r for r in sources if r['digit'] == digit], key=lambda r: r['offset'])
        tt = sorted([r for r in targets if r['digit'] == digit], key=lambda r: r['offset'])
        expected = 1 if digit == 1 else 2
        if len(ss) != expected or len(tt) != expected or (expected == 2 and
                min(ss[1]['offset']-ss[0]['offset'], tt[1]['offset']-tt[0]['offset']) < .001):
            raise ValueError('Ambiguous paired digital course assignment')
        for source, target in zip(ss, tt):
            if args.part_kind=='proper-thumb':
                if source['offset']*target['offset']<=0:
                    raise ValueError('Reference chain is on the opposite side of the thumb')
                target=copy.deepcopy(target)
                target['axis'],target['boneStationMatch']=trim_to_matched_bone_station(
                    target['axis'],source['axis'][0],
                    np.array([v['sourceCentre_m']for v in bones[digit]]),
                    np.array([v['targetCentre_m']for v in bones[digit]]))
            mouth = graph['boundary_raw'][graph['boundary_nodes'].index(source['proof']['proximalNode'])]
            centre = m.vertices[mouth].mean(0)
            radius = float(np.linalg.norm(m.vertices[mouth]-centre, axis=1).max())
            matches.append(dict(source=source, target=copy.deepcopy(target), mouthRadius_m=radius))
    fork_evidence = separate_guide_mouths(matches)
    if args.describe_matches_only:
        records=[dict(sourcePartId=source_row['partId'],sourceComponent=r['source']['owner']['physicalComponent'],
            digit=r['source']['digit'],targetChain=r['target']['chain'],sourceAxis=r['source']['axis'].tolist(),
            sourceAxisProof=r['source']['proof'],targetAxis=r['target']['axis'].tolist(),
            sourceContourAxisResolved=r['source']['owner']['physicalComponent']in axis_overrides,
            proximalGuideTrim_m=r['target'].get('proximalTrim_m',0))for r in matches]
        out.mkdir();report=dict(complete=True,courseMatches=records,guideForkSeparation=fork_evidence,
            sourceGeometryChanged=False,installed=False,installationAllowed=False,
            inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},
            sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
        write_json(out/'report.json',report)
        print(json.dumps(dict(courseMatches=len(records),curvedSourceAxes=sum(r['sourceContourAxisResolved']for r in records),sourceGeometryChanged=False)))
        return
    fixed_route_ports=False
    if args.target_axis_route:
        if args.transport_current_course:
            raise ValueError('Saved course transport cannot replace its guide with a different corridor')
        paths['target-axis-route']=Path(args.target_axis_route).resolve()
        route_inventory=json.loads(paths['target-axis-route'].read_text())
        if 'rows'in route_inventory:
            if not args.only_source_components:raise ValueError('A route batch needs an explicit completed component scope')
            requested=set(map(int,args.only_source_components.split(',')));routes=[]
            entries=[r for r in route_inventory['rows']if r['partId']==source_row['partId']and r['component']in requested]
            if {r['component']for r in entries}!=requested or not all(r['complete']for r in entries):
                raise ValueError('At least one selected route did not complete its independent stage')
            for entry in entries:
                if sha(entry['path'])!=entry['sha256']:raise ValueError('Completed route record changed')
                routes.append((json.loads(Path(entry['path']).read_text()),Path(entry['path'])))
        else:routes=[(route_inventory,paths['target-axis-route'])]
        for routed,route_record in routes:
            if (not routed['proof']['complete'] or routed['partId']!=source_row['partId'] or
                sha(routed['routePath'])!=routed['routeSha256']):
                raise ValueError('Recipient corridor is incomplete or bound to another source')
            if args.retain_source_ports and routed['parentsRequireRecoupling']:
                raise ValueError('A free-mouth route cannot retain its source ports')
            fixed_route_ports=args.retain_source_ports and not routed['parentsRequireRecoupling']
            matched=[r for r in matches if r['source']['owner']['physicalComponent']==routed['component']]
            if len(matched)!=1:raise ValueError('Recipient route needs one explicit source component owner')
            with np.load(routed['routePath'],allow_pickle=False)as z:route=z['route'].copy()
            if route.ndim!=2 or route.shape[1]!=3 or len(route)<3:raise ValueError('Finite recipient route required')
            matched[0]['target']['axis']=route
            matched[0]['target']['routedAxisEvidence']=dict(path=str(route_record),sha256=sha(route_record))
    base_mesh = Mesh.load(next(r['path'] for r in base['rows'] if r['partId'] == source_row['partId'])) if base else m
    geometry_input = base_mesh if args.transport_current_course else m
    if args.transport_current_course and not base:
        raise ValueError('Saved course transport requires an explicit base candidate')
    if not np.array_equal(base_mesh.faces, m.faces):
        raise ValueError('Base proposal changed imported source face ancestry')
    only = set(map(int, args.only_source_components.split(','))) if args.only_source_components else None
    if only is not None and not base:
        raise ValueError('A scoped correction requires a retained whole-hand candidate')
    q = base_mesh.vertices.copy()
    R = np.repeat(np.eye(3)[None], len(q), axis=0)
    covered = np.zeros(len(q), bool)
    raw_owner = np.full(len(q), -1, int)
    poses = []
    records = []
    for match in matches:
        source, target = match['source'], match['target']
        if only is not None and source['owner']['physicalComponent'] not in only:
            continue
        digit = source['digit']
        ids = np.unique(m.faces[source['owner']['sourceFaceIds']])
        if covered[ids].any():
            raise ValueError('Independently posed source components share controls')
        transport_axis=source['axis'][[0,-1]]if args.palm_chord_frame else source['axis']
        old, old_t = smooth_samples(transport_axis)
        new, new_t = smooth_samples(target['axis'])
        if args.transport_current_course:
            from geometry_tools.algorithms.port_section_graph import shortest_rotation
            prior = next(r for r in base['courseMatches'] if r['sourceComponent'] == source['owner']['physicalComponent']
                         and (r.get('sourcePartId')==source_row['partId']or
                              ('sourcePartId'not in r and np.array_equal(r['sourceAxis'],source['axis']))))
            old, old_t = smooth_samples(np.array(prior['targetAxis']))
            advance = target.get('proximalTrim_m', 0)-prior.get('proximalGuideTrim_m', 0)
            if advance <= 0:
                continue
            s = arclength(old)
            transition = min(s[-1]*.45, max(.012, 2.5*advance))
            if transition <= 1.875*advance:
                raise ValueError('Saved course is too short for a monotone finite-mouth advance')
            t = np.clip(s/transition, 0, 1); fade=t**3*(10-15*t+6*t*t)
            station = s+advance*(1-fade)
            new = np.column_stack([np.interp(station,s,old[:,k]) for k in range(3)])
            new_t = unit(np.gradient(new, axis=0))
            tail = s >= transition; new[tail] = old[tail]; new_t[tail] = old_t[tail]
        if args.retain_source_ports and not fixed_route_ports:
            new, new_t = retain_proximal_course(old, new, old_t, new_t)
        if args.pose_strength < 1:
            new = old+args.pose_strength*(new-old)
            new_t = unit(np.gradient(new, axis=0))
            if args.retain_source_ports:
                new_t[0] = old_t[0]
        old_ref,new_ref=np.array([0.,0.,1.]),np.array([0.,0.,1.])
        if args.palm_chord_frame:
            def palm_normal(key):
                centres=np.array([bones[d][0][key]for d in range(1,6)])
                _,_,vh=np.linalg.svd(centres-centres.mean(0));n=vh[-1]
                return n if n[2]>=0 else -n
            old_ref,new_ref=[palm_normal(k)for k in ['sourceCentre_m','targetCentre_m']]
        if args.bone_radial_frame:
            def bone_radial(axis,key):
                chain=np.array([v[key]for v in bones[digit]])
                centre=axis[len(axis)//2];delta=np.diff(chain,axis=0)
                u=np.clip(np.sum((centre-chain[:-1])*delta,axis=1)/np.sum(delta*delta,axis=1),0,1)
                nearest=chain[:-1]+u[:,None]*delta
                radial=centre-nearest[np.linalg.norm(nearest-centre,axis=1).argmin()]
                if np.linalg.norm(radial)<1e-6:raise ValueError('Bone chain does not determine the vessel radial side')
                return radial/np.linalg.norm(radial)
            old_ref=bone_radial(old,'sourceCentre_m');new_ref=bone_radial(new,'targetCentre_m')
        before = rotation_minimizing_frames(old, old_t, initial_reference=old_ref)
        after = rotation_minimizing_frames(new, new_t, initial_reference=new_ref)
        if bone_indices:
            before = surface_frames(old, old_t, bone_indices[0])
            after = surface_frames(new, new_t, bone_indices[1])
        rotations = after@before.transpose(0, 2, 1)
        if args.transport_current_course:
            rotations = np.array([shortest_rotation(a,b) for a,b in zip(old_t,new_t)])
            rotations[tail] = np.eye(3)
        raw_owner[ids] = len(poses)
        poses.append((old, new, rotations))
        q[ids], R[ids] = axis_transport(geometry_input.vertices[ids], old, new, rotations)
        records.append(dict(sourcePartId=source_row['partId'],sourceComponent=source['owner']['physicalComponent'], digit=digit,
                            targetChain=target['chain'], sourceOffset_m=source['offset'], targetOffset_m=target['offset'],
                            sourceAxis=source['axis'].tolist(), targetAxis=(new if args.transport_current_course else target['axis']).tolist(),
                            proximalGuideTrim_m=target.get('proximalTrim_m', 0),
                            boneStationMatch=target.get('boneStationMatch'),
                            transportAxis=transport_axis.tolist(),
                            frameConstruction='paired-bone-radial-and-source-chord'if args.bone_radial_frame else
                                              'paired-palm-plane-and-source-chord'if args.palm_chord_frame else 'source-axis',
                            sourceFrameReference=old_ref.tolist(),targetFrameReference=new_ref.tolist(),
                            sourceAxisProof=source['proof']))
        if 'routedAxisEvidence'in target:records[-1]['recipientRoute']=target['routedAxisEvidence']
        covered[ids] = True
    if only is None and not covered[np.unique(m.faces)].all():
        raise ValueError('Course matches do not own the complete imported proper-palmar surface')
    # Restore each literal source section by its own proper pose. All exact
    # attribute aliases move with it; no proximity weld or radius shrink occurs.
    _, physical = np.unique(m.vertices, axis=0, return_inverse=True)
    actual_sections=complete_source_sections(graph,physical)
    for section in actual_sections:
        if not covered[section].any():
            continue
        centre = geometry_input.vertices[section].mean(0)
        owner = np.unique(raw_owner[section])
        if len(owner) != 1 or owner[0] < 0:
            raise ValueError('A literal source profile crosses component ownership')
        old, new, rotations = poses[owner[0]]
        point, rotation = axis_transport(centre[None], old, new, rotations)
        for raw in section:
            aliases = physical == physical[raw]
            q[aliases] = (geometry_input.vertices[aliases]-centre)@rotation[0].T+point[0]
            R[aliases] = rotation[0]
    if args.retain_source_ports:
        port_ids = []
        for pair in data['source'].get('pairs', []):
            if pair['firstPart'] == source_row['partId']:
                port_ids.extend(pair['firstRaw'])
            if pair['secondPart'] == source_row['partId']:
                port_ids.extend(pair['secondRaw'])
        anchored = np.isin(physical, physical[np.array(port_ids, int)])
        q[anchored] = m.vertices[anchored]
        R[anchored] = np.eye(3)
    attrs = {k: v.copy() for k, v in base_mesh.attributes.items()}
    original_normals = geometry_input.attributes['normals'].astype(float)*geometry_input.attributes.get('normalSigns', np.ones(len(q)))[:, None]
    normals = np.einsum('nij,nj->ni', R[covered], original_normals[covered])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    attrs['normals'][covered]=np.rint(normals*32767).astype('<i2')
    if 'normalSigns'not in attrs:attrs['normalSigns']=np.ones(len(q),'i1')
    attrs['normalSigns'][covered]=1
    if args.retain_source_ports:
        attrs['normals'][anchored]=m.attributes['normals'][anchored]
        attrs['normalSigns'][anchored]=m.attributes.get('normalSigns',np.ones(len(q),'i1'))[anchored]
    candidate = Mesh(q, m.faces, m.frame, attrs, m.source)
    out.mkdir()
    file = out/f'{args.part_kind}.npz'
    candidate.save(file)
    report = copy.deepcopy(base or data['source'])
    reference_rows={r['partId']:r for r in data['source']['rows']}
    for row in report['rows']:
        if row['partId'] == source_row['partId']:
            row.update(path=str(file), sha256=sha(file), geometrySha256=candidate.geometry_sha,
                       sourceTrianglePairing=dict(path=legacy['path'], sha256=legacy['sha256'], array='retained_source_triangle_pairs'))
        reference=reference_rows[row['partId']]
        if not np.array_equal(Mesh.load(row['path']).faces,Mesh.load(reference['path']).faces):
            raise ValueError('An unrelated base member changed the explicit current-source face ancestry')
        row.update(historicalSourceComparison=dict(path=row.get('sourcePath'),sha256=row.get('sourceSha256')),
                   sourcePath=reference['path'],sourceSha256=reference['sha256'])
    for key in ['wholeAssembly', 'trials', 'constructedJoins']:
        report.pop(key, None)
    for key in ['catalogueSha256', 'bankSha256', 'sourcePackedCandidateSha256']:
        if key in report:
            report['input'+key[0].upper()+key[1:]] = report.pop(key)
    for pair in report.get('pairs', []):
        if not args.retain_source_ports and source_row['partId'] in [pair['firstPart'], pair['secondPart']]:
            pair.update(literalCoordinatesAndEffectiveNormals=False,
                        packedCoordinatesAndEffectiveNormalWordsLiteral=False,
                        requiresRecoupling=True)
    report.update(courseMatches=records, retainedImportedProfiles=len(actual_sections),
                  guideForkSeparation=fork_evidence,
                  sourceOnly=False, installed=False, installationAllowed=False, complete=True,
                  sourceCodeSha256=sha(__file__), inputs={k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()},
                  portsRequireRecoupling=not args.retain_source_ports,
                  sourceForkPortsRetained=args.retain_source_ports,
                  poseStrength=args.pose_strength,
                  sourceAnnularPropagation=args.extend_source_annuli,
                  elapsedSeconds=time.monotonic()-started)
    if base:
        updated={(r['sourcePartId'],r['sourceComponent']):r for r in records}
        report['courseMatches']=[updated.pop((r.get('sourcePartId'),r['sourceComponent']),r) for r in base.get('courseMatches',[])]
        report['courseMatches'].extend(updated.values())
    report['onlySourceComponents']=sorted(only) if only is not None else None
    report['savedDownstreamCorrectionsRetained']=args.transport_current_course
    if report['elapsedSeconds'] > 18:
        raise TimeoutError('18-second digital course proposal')
    write_json(out/'report.json', report)
    print(json.dumps(dict(matchedImportedCourses=len(records), rigidSourceSections=len(actual_sections),
                          oldGeometryCopied=False, elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
