"""Independent literal-interface and fixed-neighbor replay of a papillary field."""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.model import Mesh
from tongue_whole_input import WORK,fixed_bone_sources


def main(candidate,output):
    began=time.monotonic()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new audit stage')
    report=json.loads((candidate/'construction.json').read_text())
    if not report['passed']:raise ValueError('Failed papillary construction')
    old=json.loads(Path(report['source']['path']).read_text())
    if digest_file(report['source']['path'])!=report['source']['sha256']:raise ValueError('Source manifest changed')
    actual=json.loads((candidate/'packed/tissue.json').read_text())
    spec=importlib.util.spec_from_file_location('septum_proposal',ROOT/'server/propose-tongue-septum.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    selected,meshes=mod.inputs()
    muco=Mesh.load(candidate/'packed/TA98TONGUE_male_mucosa.npz')
    old_m=Mesh.load(WORK/'shared-coverings-exposed-apo-sl-v3/packed/TA98TONGUE_male_mucosa.npz')
    changed=report['old_source_face'];count=report['replacement_faces']
    old_tri=old_m.triangles;tri=muco.triangles
    literal=np.array_equal(tri[:changed],old_tri[:changed])and np.array_equal(tri[changed+count:],old_tri[changed+1:])
    def corner_rows(manifest,cell):
        data=Path(manifest['bank']).read_bytes();row=manifest['parts'][cell];n=row['vertexCount']
        index=np.ndarray((row['indexCount']//3,3),'<u4',buffer=data,offset=row['indices'])
        words=np.ndarray((n,3),'<i2',buffer=data,offset=row['normals'])[index]
        signs=np.ndarray((n,),'i1',buffer=data,offset=row['normalSigns'])[index]
        return words.astype(np.int32)*signs[:,:,None]
    ow=corner_rows(old,'mucosa');nw=corner_rows(actual,'mucosa')
    words_literal=np.array_equal(nw[:changed],ow[:changed])and np.array_equal(nw[changed+count:],ow[changed+1:])
    shared=[]
    for key in ['aponeurosis','superior-longitudinal']:
        name='TA98TONGUE_male_'+('intrinsic_superior-longitudinal'if key=='superior-longitudinal'else key)+'.npz'
        a=Mesh.load(candidate/'packed'/name);b=Mesh.load(WORK/'shared-coverings-exposed-apo-sl-v3/packed'/name)
        shared.append(dict(cell=key,oriented_triangle_coordinates_literal=np.array_equal(a.triangles,b.triangles),
            actual_corner_words_literal=np.array_equal(corner_rows(actual,key),corner_rows(old,key))))
    pairs=[]
    for key,other in meshes.items():
        if key in ['whole','mucosa','aponeurosis','superior-longitudinal']:continue
        _,s=pair_contacts(muco,other,0.,max_pairs=100000,max_tests=1000000)
        pairs.append(dict(second=key,summary=s,geometry_sha256=other.geometry_sha))
    stylo=json.loads((WORK/'styloglossus-v2/candidate.json').read_text())
    # Named source manifests differ; use the actual review-spec mesh records,
    # not guessed positions or an unreviewed older styloglossus version.
    review=json.loads((WORK/'shared-coverings-exposed-apo-sl-v3/review-spec.json').read_text())
    for row in review['meshes']:
        if 'styloglossus'in row['name'].lower():
            other=Mesh.load(row['mesh']);_,s=pair_contacts(muco,other,0.,max_pairs=100000,max_tests=1000000)
            pairs.append(dict(second=row['name'],summary=s,mesh=row['mesh'],sha256=digest_file(row['mesh'])))
    for row in fixed_bone_sources():
        other=Mesh.load(row['meshSavePath']);_,s=pair_contacts(muco,other,0.,max_pairs=100000,max_tests=1000000)
        pairs.append(dict(second=row['partId'],summary=s,geometry_sha256=other.geometry_sha))
    clear=all(r['summary']['complete_within_pair_budget']and not(
        r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs'])for r in pairs)
    passed=bool(actual['passed']and literal and words_literal and clear and all(
        r['oriented_triangle_coordinates_literal']and r['actual_corner_words_literal']for r in shared))
    output.mkdir(parents=True)
    atomic_json(output/'audit.json',dict(schema='tongue-papillary-independent-replay/v1',passed=passed,
        construction=dict(path=str(candidate/'construction.json'),sha256=digest_file(candidate/'construction.json')),
        packed_manifest=dict(path=str(candidate/'packed/tissue.json'),sha256=digest_file(candidate/'packed/tissue.json')),
        unchanged_mucosal_faces_literal=literal,unchanged_mucosal_raw_normal_words_literal=words_literal,
        all_inner_and_shared_cells_literal=shared,fixed_and_adjoining_pairs=pairs,
        complete_pair_count=len(pairs),actual_finite_tissue_checks_passed=actual['passed'],
        field_diameter_mm=3.6,filiform_units=13,fungiform_units=1,
        filiform_mean_diameter_um=355,fungiform_diameter_um=878,
        actual_maximum_relief_um=report['maximum_relief_m']*1e6,
        added_physical_meshes=0,image_textures=0,extra_common_bank_bytes=actual['total_bytes']-old['total_bytes'],
        full_dorsum_coverage=False,anatomical_acceptance=False,installed=False,
        elapsed_seconds=time.monotonic()-began))
    print(json.dumps(dict(passed=passed,complete_pairs=len(pairs),unchanged_geometry_and_words_literal=literal and words_literal,
        extra_common_bank_bytes=actual['total_bytes']-old['total_bytes'],elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--candidate',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();main(a.candidate.resolve(),a.output.resolve())
