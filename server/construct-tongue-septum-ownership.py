"""Give a contained MLS prototype one reciprocal T/V tissue interface.

No anatomy is fitted and no atlas is modified. When two simple closed source
boundaries do not intersect and the inner one is inside the outer one, their
set difference is represented exactly by retaining the outer chart and using
the inner chart with opposite orientation. No Boolean remeshing is necessary.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
import igl
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.tissue import pack,normal_words
from tongue_whole_input import WORK,fixed_bone_sources
_spec=importlib.util.spec_from_file_location('tongue_septum_proposal',Path(__file__).with_name('propose-tongue-septum.py'))
_module=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_module)
inputs=_module.inputs


def main(prototype,output):
    began=time.monotonic()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new local stage')
    report=json.loads((prototype/'candidate.json').read_text())
    if not report['passed_geometry']or report['installed']:raise ValueError('Use passed uninstalled prototype')
    bindings,meshes=inputs()
    if bindings!=report['inputs']:raise ValueError('Selected geometry changed')
    source=Path(report['mesh']['path'])
    if digest_file(source)!=report['mesh']['sha256']:raise ValueError('Changed prototype')
    septum=Mesh.load(source);tv=meshes['transverse-vertical-shared']
    origin=np.array([0.,1.5,0.])
    if not np.array_equal((tv.vertices-origin).astype('<f4').astype(float)+origin,tv.vertices):
        raise ValueError('Do not reencode source T/V controls')
    positions=np.r_[tv.vertices,septum.vertices]
    actual=(positions-origin).astype('<f4').astype(float)+origin
    corrected=Mesh(actual[len(tv.vertices):],septum.faces,septum.frame)
    # Recheck the actual common-bank Float32 position field before sharing it.
    _,summary=pair_contacts(tv,corrected,0.,max_pairs=100000,max_tests=1000000)
    q=np.r_[corrected.vertices,corrected.triangles.mean(1)]
    winding=np.abs(igl.winding_number(tv.vertices,tv.faces,q))
    if not summary['complete_within_pair_budget']or summary['crossing_face_pairs']or summary['coplanar_overlap_face_pairs']or not np.all(winding>.99):
        raise ValueError('The actual shared-bank septum is not strictly contained')
    surfaces={'transverse-vertical:source-outer':tv.faces,
              'septum:reciprocal':septum.faces+len(tv.vertices)}
    cells={'transverse-vertical-shared':[
        dict(surface='transverse-vertical:source-outer',side=1,neighbor=None),
        dict(surface='septum:reciprocal',side=-1,neighbor='septum')],
        'septum':[dict(surface='septum:reciprocal',side=1,neighbor='transverse-vertical-shared')]}
    words=normal_words(actual-origin,surfaces)
    with np.load(bindings['transverse-vertical-shared']['path'],allow_pickle=False)as raw:
        if 'corner_normals'not in raw:raise ValueError('Source T/V artist chart absent')
        words['transverse-vertical:source-outer']=raw['corner_normals'].copy()
    output.mkdir(parents=True)
    result=pack(positions,surfaces,cells,
        {'transverse-vertical-shared':'TA98TONGUE:male:intrinsic:transverse-vertical-shared',
         'septum':'TA98TONGUE:male:septum'},tv.frame,output/'packed',origin=origin,corner_words=words,
        construction=dict(kind='strict-contained-source-boundary-complement',
            finite_wall_cells_pass=True,prototype_sha256=digest_file(prototype/'candidate.json'),
            source_TV_outer_literal=True,partial_anterior_extent=True,anatomical_acceptance=False))
    if not result['passed']:
        atomic_json(output/'ownership.json',dict(schema='tongue-septum-shared-ownership/v1',
            passed=False,anatomical_acceptance=False,installed=False,installation_ready=False,
            prototype=dict(path=str(prototype/'candidate.json'),sha256=digest_file(prototype/'candidate.json')),
            actual_tissue_manifest=dict(path=str(output/'packed/tissue.json'),sha256=digest_file(output/'packed/tissue.json')),
            state='stopped-at-occupied-volume-gate',checks=result['checks'],
            no_geometry_parameters_or_thresholds_relaxed=True,
            source_files_modified=False,elapsed_seconds=time.monotonic()-began))
        raise ValueError('The actual common tissue bank failed: '+str(result['checks']['failures']))
    actual_septum=Mesh.load(output/'packed/TA98TONGUE_male_septum.npz')
    actual_tv=Mesh.load(output/'packed/TA98TONGUE_male_intrinsic_transverse-vertical-shared.npz')
    # Exact source face coordinates are an explicit invariant of this complement.
    outer_literal=np.array_equal(actual_tv.triangles[:len(tv.faces)],tv.triangles)
    reciprocal=np.array_equal(actual_tv.triangles[len(tv.faces):, [0,2,1]],actual_septum.triangles)
    checks=[]
    for name,mesh in meshes.items():
        if name in ['whole','transverse-vertical-shared']:continue
        _,s=pair_contacts(actual_septum,mesh,0.,max_pairs=100000,max_tests=1000000)
        checks.append(dict(second=name,summary=s))
    for row in fixed_bone_sources():
        _,s=pair_contacts(actual_septum,Mesh.load(row['meshSavePath']),0.,max_pairs=100000,max_tests=1000000)
        checks.append(dict(second=row['partId'],summary=s))
    clear=all(r['summary']['complete_within_pair_budget']and not(
        r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs'])for r in checks)
    passed=bool(result['passed']and outer_literal and reciprocal and clear)
    atomic_json(output/'ownership.json',dict(schema='tongue-septum-shared-ownership/v1',passed=passed,
        prototype=dict(path=str(prototype/'candidate.json'),sha256=digest_file(prototype/'candidate.json')),
        actual_tissue_manifest=dict(path=str(output/'packed/tissue.json'),sha256=digest_file(output/'packed/tissue.json')),
        actual_boundary_containment=summary,minimum_septum_absolute_winding=float(winding.min()),
        source_TV_outer_coordinates_and_indices_literal=outer_literal,
        reciprocal_744_facets_literal=reciprocal,
        source_TV_outer_raw_corner_words_literal=True,
        actual_packed_tissue_checks=result['checks'],fixed_and_tissue_contacts=checks,
        source_files_modified=False,source_fit_or_remeshing_performed=False,
        total_packed_bytes=result['total_bytes'],image_textures=0,
        anatomical_acceptance=False,installed=False,installation_ready=False,
        remaining=['Partial posterior/body septum extent; anterior extent requires anatomical reconstruction',
                   'Epimysial fiber attachment architecture remains an illustrative finite interface',
                   'Actual Blender anatomical surface and section review required'],
        elapsed_seconds=time.monotonic()-began))
    print(json.dumps(dict(passed=passed,source_TV_outer_literal=outer_literal,
        reciprocal_facets=len(septum.faces),total_packed_bytes=result['total_bytes'],
        elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--prototype',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    main(args.prototype.resolve(),args.output.resolve())
