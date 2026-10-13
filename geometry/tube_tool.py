"""Command-line tools; all NPZ coordinates are metres in one declared frame.

Examples:
  python3 -m geometry.tube_tool bend --waypoints turns.json --radius-mm 30 --output path.npz
  python3 -m geometry.tube_tool sweep --source wall.npz --path path.npz --ring-size 12 --output bent.npz
  python3 -m geometry.tube_tool sweep --source wall.npz --path edited-centres.npz --ring-size 12 --keep-source-roll --output local-bend.npz
  python3 -m geometry.tube_tool audit --source bent.npz --report checks.json
  python3 -m geometry.tube_tool join --first a.npz --second b.npz --mouth-a a.json --mouth-b b.json --output joined.npz

source NPZ: fields[interfaces,vertices,3], faces[n,3], optional geometryOrigin[3].
path NPZ: centres[stations,3], optional tangents[stations,3], optional geometryOrigin[3].
Mouth JSON files contain explicitly ordered vertex indices. No atlas is modified.
Exit 2 means a geometry gate failed; a diagnostic report remains available.
"""
from pathlib import Path
import argparse, json, sys, importlib.util
import numpy as np
from .tube_sections import (TubeGeometryError,TubeTolerances,filleted_path,centripetal_curve,
    resample_path,sweep_sections,bend_diagnostics,path_tangents)
from .tube_continuity import (continuity_audit,finite_wall_cells,join_mouths,loft_faces)


def _write_report(path, report):
    if path:
        p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2,allow_nan=False))


def _full_mesh_checks(fields, faces, validator):
    """Use an explicit existing narrow-phase validator, never silently skip it."""
    if not validator:return dict(whole_mesh_intersections_checked=False)
    path=Path(validator);spec=importlib.util.spec_from_file_location('tube_narrow_phase',path)
    if spec is None or spec.loader is None:raise TubeGeometryError('cannot load the specified mesh validator')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    if not callable(getattr(module,'sweep_hits',None)):raise TubeGeometryError('validator must expose sweep_hits')
    checks=[]
    for k in range(len(fields)):
        for j in range(k,len(fields)):
            hits=module.sweep_hits(fields[k][faces],fields[j][faces],selfcheck=k==j,max_examples=4)
            checks.append(dict(first=k,second=j,intersections=hits['properOrCoplanarIntersections'],examples=hits['examples']))
    return dict(whole_mesh_intersections_checked=True,interface_pair_checks=checks,
        interface_intersection_pass=not any(r['intersections'] for r in checks))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    sub=parser.add_subparsers(dest='operation',required=True)
    bend=sub.add_parser('bend',help='rounded guide, with controllable bend radius')
    bend.add_argument('--waypoints',required=True);bend.add_argument('--radius-mm',type=float,default=30)
    bend.add_argument('--method',choices=['fillet','centripetal'],default='fillet')
    bend.add_argument('--spacing-mm',type=float,default=2);bend.add_argument('--chord-tolerance-mm',type=float,default=.1)
    bend.add_argument('--output',required=True);bend.add_argument('--report')
    sweep=sub.add_parser('sweep',help='transport all sections/wall layers together')
    sweep.add_argument('--source',required=True);sweep.add_argument('--path',required=True)
    sweep.add_argument('--output',required=True);sweep.add_argument('--report');sweep.add_argument('--ring-size',type=int,required=True)
    sweep.add_argument('--first-control',type=int,default=0);sweep.add_argument('--area-tolerance',type=float,default=.01)
    sweep.add_argument('--diameter-tolerance',type=float,default=.01);sweep.add_argument('--thickness-tolerance',type=float,default=.01)
    sweep.add_argument('--shape-tolerance-mm',type=float,default=.0001);sweep.add_argument('--planarity-tolerance-mm',type=float,default=.1)
    sweep.add_argument('--curvature-support-limit',type=float,default=.85);sweep.add_argument('--turn-step-deg',type=float,default=15)
    sweep.add_argument('--length-tolerance',type=float);sweep.add_argument('--mesh-validator')
    sweep.add_argument('--keep-source-roll',action='store_true',help='local course adjustment; unchanged sections stay literal, all layers share each rigid transform')
    audit=sub.add_parser('audit',help='real carrier topology and every finite wall interval')
    audit.add_argument('--source',required=True);audit.add_argument('--report');audit.add_argument('--mesh-validator')
    audit.add_argument('--expected-mouths',type=int,default=2)
    audit.add_argument('--expected-collars',help='JSON list of expected 3D mouth collars in the source frame')
    audit.add_argument('--collar-tolerance-mm',type=float,default=0)
    join=sub.add_parser('join',help='join explicit matching open collars, never fill the lumen')
    join.add_argument('--first',required=True);join.add_argument('--second',required=True)
    join.add_argument('--mouth-a',required=True);join.add_argument('--mouth-b',required=True)
    join.add_argument('--position-tolerance-mm',type=float,default=.001);join.add_argument('--area-tolerance',type=float,default=.01)
    join.add_argument('--tangent-tolerance-deg',type=float,default=15)
    join.add_argument('--output',required=True);join.add_argument('--report');join.add_argument('--mesh-validator')
    args=parser.parse_args(argv)
    # A diagnostic tool must never overwrite its input bank or path in place.
    inputs=[Path(getattr(args,k)).resolve() for k in ('source','first','second','path','waypoints','mouth_a','mouth_b','expected_collars') if getattr(args,k,None)]
    outputs=[Path(getattr(args,k)).resolve() for k in ('output','report') if getattr(args,k,None)]
    if any(p in inputs for p in outputs):
        print(json.dumps(dict(construction_checks_pass=False,error='output/report must differ from every input',atlas_modified=False)))
        return 2
    if getattr(args,'output',None):
        if Path(args.output).suffix!='.npz':parser.error('--output must use the .npz extension')
        Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    try:
        if args.operation=='bend':
            waypoints=json.loads(Path(args.waypoints).read_text())
            if isinstance(waypoints,dict):waypoints=waypoints['waypoints']
            if args.method=='fillet':p,t=filleted_path(waypoints,args.radius_mm/1000,args.spacing_mm/1000)
            else:p,t=centripetal_curve(waypoints,args.chord_tolerance_mm/1000,args.spacing_mm/1000)
            np.savez_compressed(args.output,centres=p,tangents=t)
            report=dict(constructed=True,method=args.method,stations=len(p),radius_m=args.radius_mm/1000 if args.method=='fillet' else None,
                source_units='metres',tube_section_checked=False,anatomical_acceptance=False)
            _write_report(args.report,report);return 0
        if args.operation=='sweep':
            z=np.load(args.source);fields=z['fields'];origin=z['geometryOrigin'] if 'geometryOrigin' in z else np.zeros(3)
            start=args.first_control
            if start<0 or start>=fields.shape[1] or args.ring_size<3 or (fields.shape[1]-start)%args.ring_size:
                raise TubeGeometryError('controls do not form the requested ordered section bank')
            profiles=fields[:,start:].reshape(len(fields),-1,args.ring_size,3)
            path=np.load(args.path);centres=path['centres'];tangents=path['tangents'] if 'tangents' in path else None
            pathorigin=path['geometryOrigin'] if 'geometryOrigin' in path else origin
            centres=centres+pathorigin-origin
            if len(centres)!=profiles.shape[1]:centres,tangents=resample_path(centres,profiles.shape[1],tangents)
            tol=TubeTolerances(area_relative=args.area_tolerance,diameter_relative=args.diameter_tolerance,
                thickness_relative=args.thickness_tolerance,absolute_m=args.shape_tolerance_mm/1000,
                planarity_m=args.planarity_tolerance_mm/1000,curvature_support_limit=args.curvature_support_limit,
                max_station_turn_rad=np.deg2rad(args.turn_step_deg),maximum_length_change_relative=args.length_tolerance)
            if args.keep_source_roll:
                # Only the centres define a local redirection. A stored analytic
                # tangent belongs to the complete-guide sweep, not this mode.
                tangents=None
            result,frames,report=sweep_sections(profiles,centres,tangents,tol,preserve_source_roll=args.keep_source_roll)
            bank=result.reshape(len(fields),-1,3)
            if 'faces' in z:
                faces=z['faces'];mask=np.all(faces>=start,axis=1);faces=faces[mask]-start
                if not len(faces):raise TubeGeometryError('no carrier faces after the control crop')
            else:faces=loft_faces(profiles.shape[1],args.ring_size);mask=None
            report['finite_wall']=finite_wall_cells(bank,faces)
            report['carrier']=continuity_audit(bank[0],faces)
            report.update(_full_mesh_checks(bank,faces,args.mesh_validator))
            data=dict(fields=bank,faces=faces,centres=centres,frames=frames,tangents=path_tangents(centres) if tangents is None else tangents,geometryOrigin=origin)
            if 'owners' in z and mask is not None:data['owners']=z['owners'][mask]
            np.savez_compressed(args.output,**data)
            passed=report['section_transport_pass'] and report['finite_wall']['finite_wall_cells_pass'] and report['carrier']['topological_continuity_pass'] and report.get('interface_intersection_pass',True)
        elif args.operation=='join':
            a,b=np.load(args.first),np.load(args.second)
            origin=a['geometryOrigin'] if 'geometryOrigin' in a else np.zeros(3)
            borigin=b['geometryOrigin'] if 'geometryOrigin' in b else np.zeros(3)
            am=json.loads(Path(args.mouth_a).read_text());bm=json.loads(Path(args.mouth_b).read_text())
            fields,faces,mapping,report=join_mouths(a['fields'],a['faces'],am,
                b['fields'].astype(float)+borigin-origin,b['faces'],bm,args.position_tolerance_mm/1000,args.area_tolerance,np.deg2rad(args.tangent_tolerance_deg))
            report['finite_wall']=finite_wall_cells(fields,faces)
            report.update(_full_mesh_checks(fields,faces,args.mesh_validator))
            np.savez_compressed(args.output,fields=fields,faces=faces,geometryOrigin=origin,secondControlMapping=mapping)
            passed=report['finite_wall']['finite_wall_cells_pass'] and report.get('interface_intersection_pass',True)
        else:
            z=np.load(args.source);fields=z['fields'];faces=z['faces']
            collars=json.loads(Path(args.expected_collars).read_text()) if args.expected_collars else None
            report=dict(carrier=continuity_audit(fields[0],faces,args.expected_mouths,expected_collars=collars,collar_tolerance_m=args.collar_tolerance_mm/1000),finite_wall=finite_wall_cells(fields,faces))
            report.update(_full_mesh_checks(fields,faces,args.mesh_validator))
            passed=report['carrier']['topological_continuity_pass'] and report['finite_wall']['finite_wall_cells_pass'] and report.get('interface_intersection_pass',True)
        report.update(construction_checks_pass=bool(passed),anatomical_acceptance=False,atlas_modified=False)
        _write_report(args.report,report);return 0 if passed else 2
    except (TubeGeometryError,OSError,KeyError,ValueError) as err:
        _write_report(getattr(args,'report',None),dict(construction_checks_pass=False,error=str(err),atlas_modified=False))
        return 2


if __name__=='__main__':sys.exit(main())
