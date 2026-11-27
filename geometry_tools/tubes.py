"""Rounded courses, section-preserving transport and actual open-mouth joins."""
from pathlib import Path
import json
import numpy as np
from .common import atomic_json,atomic_npz,digest_file
from .model import Mesh
from .audit import audit
from .contact import pair_contacts
from .distance import triangle_pairs


def load_fields(path,frame=None):
    with np.load(path,allow_pickle=False)as z:
        actual_frame=frame or(str(z['frame'])if 'frame'in z else None)
        if not actual_frame:raise ValueError('An explicit frame is required for every tube bank')
        fields=z['fields'].astype(float)+(z['geometryOrigin']if 'geometryOrigin'in z else 0);faces=z['faces'].copy()
        if fields.ndim!=3 or fields.shape[2]!=3 or len(fields)<2 or not np.isfinite(fields).all():raise ValueError('At least two finite corresponding wall interfaces required')
    return fields,faces,actual_frame


def audit_fields(fields,faces,frame,*,mouths=2,max_pairs=250000,centres=None,minimum_lumen_radius_m=0.):
    from geometry_tools.algorithms.tube_continuity import continuity_audit,finite_wall_cells
    sheets={str(k):audit(Mesh(v,faces,frame),max_pairs=max_pairs)for k,v in enumerate(fields)}
    carriers=[continuity_audit(v,faces,expected_mouths=mouths)for v in fields]
    finite=finite_wall_cells(fields,faces);cross=[]
    for i in range(len(fields)):
        for j in range(i+1,len(fields)):
            arrays,summary=pair_contacts(Mesh(fields[i],faces,frame),Mesh(fields[j],faces,frame),0.,max_pairs=max_pairs)
            cross.append(dict(first=i,second=j,**summary,passed=not summary['truncated']and len(arrays['pairs'])==0))
    lumen=None
    if centres is not None:
        from .spatial import BVH
        center=np.asarray(centres,float);sheet=Mesh(fields[-1],faces,frame);bvh=BVH(sheet)
        if center.ndim!=2 or center.shape[1]!=3 or len(center)<2:raise ValueError('A lumen guide needs at least two metre points')
        # Segments are represented as degenerate triangles only for distance;
        # edge/face tests still detect a segment piercing a sheet interior.
        gaps=[bvh.closest_point(q)['distance_m']for q in center];hits=[]
        for k,(a,b)in enumerate(zip(center[:-1],center[1:])):
            line=Mesh(np.array([a,b,a]),np.array([[0,1,2]]),frame);arr,r=pair_contacts(line,sheet,0.,max_pairs=max_pairs)
            if len(arr['pairs'])or r['truncated']:hits.append(k)
        lumen=dict(guide_segments=len(center)-1,guide_wall_crossing_segments=hits,minimum_sampled_radius_m=float(min(gaps)),
                   required_minimum_radius_m=minimum_lumen_radius_m,passed=not hits and min(gaps)>=minimum_lumen_radius_m,
                   scope='Entire guide segments tested against inner sheet; radius bound measured at guide samples, not a full inscribed-radius proof')
    passed=all(x['passed']for x in sheets.values())and all(x['topological_continuity_pass']for x in carriers)and finite['finite_wall_cells_pass']and all(x['passed']for x in cross)and(lumen is None or lumen['passed'])
    return dict(passed=passed,sheets=sheets,carriers=carriers,finite_wall=finite,interface_separation=cross,lumen_guide=lumen,anatomical_acceptance=False)


def path(spec,output):
    from geometry_tools.algorithms.tube_sections import filleted_path,centripetal_curve
    if not spec.get('frame'):raise ValueError('An explicit frame is required')
    if spec.get('method','fillet')=='fillet':q,t=filleted_path(spec['waypoints'],spec['bend_radius_m'],spec.get('spacing_m',.002))
    else:q,t=centripetal_curve(spec['waypoints'],spec.get('chord_tolerance_m',.0001),spec.get('spacing_m',.002))
    atomic_npz(output,centres=q,tangents=t,frame=np.array(spec['frame']))
    return dict(passed=True,path=str(Path(output).resolve()),sha256=digest_file(output),stations=len(q),length_m=float(np.linalg.norm(np.diff(q,axis=0),axis=1).sum()),anatomical_acceptance=False)


def sweep(spec,output):
    from geometry_tools.algorithms.tube_sections import TubeTolerances,sweep_sections,resample_path
    fields,faces,frame=load_fields(spec['source'],spec.get('frame'));ring=int(spec['ring_size'])
    if ring<3 or fields.shape[1]%ring:raise ValueError('Tube controls must be explicitly ordered into complete corresponding rings')
    profiles=fields.reshape(len(fields),-1,ring,3)
    with np.load(spec['path'],allow_pickle=False)as z:
        if str(z['frame'])!=frame:raise ValueError('Path and tube frames differ')
        centers=z['centres'].copy();tangents=z['tangents'].copy()if 'tangents'in z else None
    if len(centers)!=profiles.shape[1]:centers,tangents=resample_path(centers,profiles.shape[1],tangents)
    tol=TubeTolerances(**spec.get('tolerances',{}));keep=spec.get('keep_source_roll',False)
    result,frames,transport=sweep_sections(profiles,centers,None if keep else tangents,tol,preserve_source_roll=keep)
    result=result.reshape(fields.shape)
    # Replay the exact Float32 coordinates intended for the viewer.
    origin=np.asarray(spec.get('geometry_origin',[0,0,0]),float);encoded=(result-origin).astype('<f4');actual=encoded.astype(float)+origin
    checks=audit_fields(actual,faces,frame,max_pairs=spec.get('max_pairs',250000),centres=centers,minimum_lumen_radius_m=spec.get('minimum_lumen_radius_m',0.))
    from geometry_tools.algorithms.tube_sections import section_metrics
    before=section_metrics(profiles);after=section_metrics(actual.reshape(profiles.shape));areaerr=float(np.max(abs(after['area']-before['area'])/before['area']))
    checks['encoded_area_change_relative']=areaerr;checks['passed']=checks['passed']and areaerr<=tol.area_relative and transport['section_transport_pass']
    atomic_npz(output,fields=encoded,faces=faces,centres=centers,frames=frames,frame=np.array(frame),geometryOrigin=origin)
    out=dict(schema='tube-transport/v1',passed=checks['passed'],output=str(Path(output).resolve()),sha256=digest_file(output),source_sha256=digest_file(spec['source']),path_sha256=digest_file(spec['path']),transport=transport,checks=checks,atlas_modified=False,anatomical_acceptance=False)
    atomic_json(str(output)+'.json',out);return out


def join(spec,output):
    from geometry_tools.algorithms.tube_continuity import join_mouths
    a,fa,af=load_fields(spec['first'],spec.get('frame'));b,fb,bf=load_fields(spec['second'],spec.get('frame'))
    if af!=bf:raise ValueError('Adjoining tubes must use one declared frame')
    fields,faces,mapping,proof=join_mouths(a,fa,spec['first_mouth'],b,fb,spec['second_mouth'],spec.get('position_tolerance_m',0.),spec.get('area_tolerance_relative',.01),spec.get('tangent_tolerance_rad',.2617993878))
    checks=audit_fields(fields,faces,af,max_pairs=spec.get('max_pairs',250000))
    atomic_npz(output,fields=fields,faces=faces,second_control_mapping=mapping,frame=np.array(af))
    out=dict(passed=checks['passed'],output=str(Path(output).resolve()),sha256=digest_file(output),join=proof,checks=checks,central_lumen_capped=False,atlas_modified=False,anatomical_acceptance=False)
    atomic_json(str(output)+'.json',out);return out


def add_parser(sub):
    p=sub.add_parser('tube');s=p.add_subparsers(dest='action',required=True)
    for action in ['path','sweep','join']:
        q=s.add_parser(action);q.add_argument('spec');q.add_argument('--output',required=True)
    q=s.add_parser('audit');q.add_argument('source');q.add_argument('--frame');q.add_argument('--mouths',type=int,default=2);q.add_argument('--max-pairs',type=int,default=250000)
    p.set_defaults(handler=run)


def run(args,db):
    if args.action=='audit':
        fields,faces,frame=load_fields(args.source,args.frame);return audit_fields(fields,faces,frame,mouths=args.mouths,max_pairs=args.max_pairs)
    return {'path':path,'sweep':sweep,'join':join}[args.action](json.loads(Path(args.spec).read_text()),args.output)
