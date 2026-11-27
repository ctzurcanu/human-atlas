"""Fit a hash-bound source-cell tube course; save every phase, never install.

Usage: python -m geometry_tools.commands.fit_tube_course SPEC.json --output NEW_DIRECTORY
Run through geometry_tools.runtime with a finite deadline and memory limit.
Y is fixed; X/Z fit in the same declared source/context coordinate frame.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from geometry_tools.course_io import sha,load_course_spec,write_json
from geometry_tools.algorithms.constrained_course_fit import fit_course
from geometry_tools.algorithms.surface_deformation_guards import surface_deformation_guards


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('spec');parser.add_argument('--output',required=True);parser.add_argument('--resume',help='Continue a terminal saved fit with the exact same immutable inputs');parser.add_argument('--reference-fit',help='Use a source-identical terminal field as a mathematical reference for added source-surface guards');args=parser.parse_args();start=time.monotonic();out=Path(args.output).resolve()
    if args.resume and args.reference_fit:raise ValueError('Choose exact resume or a new guarded reference, not both')
    if out.exists():raise ValueError('Use a new output directory; prior evidence is immutable')
    spec,g,source,plan,course,active,fixed,curve_edges,reference=load_course_spec(args.spec);out.mkdir();bindings=dict(specSha256=sha(args.spec),sourceGraphSha256=spec['source_graph_sha256'],graphArraysSha256=spec['graph_arrays_sha256'],constraintsSha256=spec['constraints_sha256'],catalogueSha256=spec['catalogue_sha256'],sourceCodeSha256=sha(__file__),solverCodeSha256=sha(ROOT/'geometry_tools/algorithms/bounded_cell_fairing.py'),orchestratorCodeSha256=sha(ROOT/'geometry_tools/algorithms/constrained_course_fit.py'));saved=[]
    start_strict=False;surface_guards=None
    if spec.get('maximum_source_surface_gradient')is not None:
        surface_guards=surface_deformation_guards(g['vertices'],g['faces'],g['cells'],active-fixed,maximum_norm=spec['maximum_source_surface_gradient']);bindings['surfaceGuardCodeSha256']=sha(ROOT/'geometry_tools/algorithms/surface_deformation_guards.py');bindings['surfaceGuardFaces']=len(surface_guards['cell_ids'])
    if args.resume or args.reference_fit:
        parent=Path(args.resume or args.reference_fit).resolve();previous=json.loads((parent/'run.json').read_text())
        if previous['state']!='finished'or previous.get('displacementsSha256')!=sha(parent/previous['displacements']):raise ValueError('Reference requires a terminal saved field')
        if args.resume and previous['bindings']['specSha256']!=bindings['specSha256']:raise ValueError('Resume changed immutable inputs')
        if args.reference_fit:
            if any(previous['bindings'][key]!=bindings[key]for key in ['sourceGraphSha256','graphArraysSha256','constraintsSha256','catalogueSha256']):raise ValueError('Guarded reference changed geometry or tissue constraints')
            previous_spec=json.loads(Path(previous['specPath']).read_text())
            if sorted(previous_spec['moving_parts'])!=sorted(spec['moving_parts']):raise ValueError('Guarded reference changed the moving closure')
        with np.load(parent/previous['displacements'],allow_pickle=False)as arrays:reference=g['centres']+arrays['delta']
        selected=previous.get('selectedStep',len(previous['steps'])-1);selected=len(previous['steps'])-1 if selected is None else selected;start_strict=bool(args.resume and previous['steps'][selected]['proof'].get('curvaturePhaseSlack_per_m',1.)<=.01);bindings['resumedFitSha256'if args.resume else'referenceFitSha256']=sha(parent/'run.json')
    write_json(out/'run.json',dict(state='running',bindings=bindings,anatomicalAcceptance=False,installed=False))
    def save_step(step,delta):
        row=dict(step)
        if delta is not None:
            name=f'phase-{step["index"]:02}.npz';np.savez_compressed(out/name,delta=delta);row.update(displacements=name,displacementsSha256=sha(out/name))
        write_json(out/f'phase-{step["index"]:02}.json',row);saved.append(row)
        p=row['proof'];print(json.dumps(dict(phase=row['index'],status=p['status'],qualified=row['qualifiedCourse'],hardViolation_m=p.get('maximumHalfspaceViolation_m'),geometricViolation_per_m=p.get('maximumGeometricCurvatureViolation_per_m'),curvatureSlack_per_m=p.get('curvaturePhaseSlack_per_m'),elapsedSeconds=row['elapsedSeconds'])),flush=True)
    delta,report=fit_course(g['centres'],g['edges'],active,fixed,{int(k):v for k,v in plan['targets'].items()},profile_cells=g['profile_cells'],halfspaces=plan['halfspaces'],plane_count=plan['planeCount'],curvature_edges=curve_edges,curvature_limits_per_m=plan['calibreCurvaturePlan']['limits_per_m'],reference_centres=reference,maximum_motion_m=spec['maximum_motion_m'],maximum_steps=spec['maximum_steps'],step_seconds=spec['step_seconds'],wall_seconds=spec['wall_seconds'],on_step=save_step,start_strict=start_strict,surface_guards=surface_guards)
    if sha(spec['catalogue'])!=spec['catalogue_sha256']:report.update(qualifiedCourse=False,complete=False,reason='Catalogue changed during fitting')
    report.update(bindings=bindings,steps=saved,allProtectedCells=sorted(fixed),specPath=str(Path(args.spec).resolve()),state='finished',elapsedSecondsIncludingInputs=time.monotonic()-start)
    if delta is not None:
        name='candidate-displacements.npz'if report['qualifiedCourse']else'diagnostic-displacements.npz';np.savez_compressed(out/name,delta=delta);report.update(displacements=name,displacementsSha256=sha(out/name))
    write_json(out/'run.json',report);print(json.dumps(dict(output=str(out),qualifiedCourse=report['qualifiedCourse'],reason=report['reason'],elapsedSeconds=report['elapsedSecondsIncludingInputs'])))
    if not report['qualifiedCourse']:raise SystemExit(1)


if __name__=='__main__':main()
