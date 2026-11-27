"""Hash-bound, once-per-job array loading for reusable course tools."""
from pathlib import Path
import hashlib,json
import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_arrays(path):
    with np.load(path,allow_pickle=False)as archive:return {k:archive[k]for k in archive.files}


def write_json(path,value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp');temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def load_course_spec(path):
    path=Path(path).resolve();spec=json.loads(path.read_text())
    if spec.get('schema')!='constrained-tube-course/v1':raise ValueError('Expected constrained-tube-course/v1 spec')
    for key in ['graph_arrays','source_graph','constraints','course','catalogue','reference','source_curvature_arrays']:
        if sha(spec[key])!=spec[key+'_sha256']:raise ValueError('Stale source: '+key)
    graph=load_arrays(spec['graph_arrays']);source=json.loads(Path(spec['source_graph']).read_text());plan=json.loads(Path(spec['constraints']).read_text());course=json.loads(Path(spec['course']).read_text())
    if not(plan['complete']and plan.get('completedRequestedSections',False)and not plan['unsupportedSupportPlanes']):raise ValueError('Incomplete physical constraint proposal')
    if set(plan['coveredSourceIndices'])!={r['index']for r in course['rows']}:raise ValueError('Constraint stations do not cover the whole proposed course')
    if not plan.get('calibreCurvaturePlan'):raise ValueError('Source-calibre geometric quality limits required')
    if spec.get('fixed_axis')!=1:raise ValueError('This fitter preserves Y; transform the complete source/context frame together before fitting another axis')
    if plan['graphArraysSha256']!=spec['graph_arrays_sha256']or plan['sourceGraphSha256']!=spec['source_graph_sha256']or plan['courseReportSha256']!=spec['course_sha256']or source['sourceProfileSha256']!=spec['catalogue_sha256']:
        raise ValueError('Geometry, physical context or catalogue belongs to another revision')
    moving=set(spec['moving_parts']);records={r['partId']:r for r in source['records']}
    if not moving or not moving<=set(records):raise ValueError('Explicit source-owned moving parts required')
    allowed=set()
    for key in moving:
        record=records[key];allowed.update(map(int,graph['cells'][record['offset']:record['offset']+record['vertices']]))
    active=allowed&set(map(int,graph['moving_cells']));fixed=set(map(int,graph['fixed_cells']))|(set(range(len(graph['centres'])))-active)
    old=load_arrays(spec['source_curvature_arrays'])
    if not np.array_equal(graph['centres'][:len(old['centres'])],old['centres']):raise ValueError('Refinement changed original source centre identities')
    reference=load_arrays(spec['reference'])['centres']
    if reference.shape!=graph['centres'].shape or not np.isfinite(reference).all():raise ValueError('Invalid mathematical reference')
    if spec['maximum_motion_m']!=plan['maximumMotion_m']:raise ValueError('Motion domain differs from physical constraints')
    return spec,graph,source,plan,course,active,fixed,old['edges'],reference
