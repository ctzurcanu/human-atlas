"""Bounded, read-only audit of explicit tongue ROI NPZs against fixed native bones."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
FRAME='model:local-ta98:assembled-metres'


def main(candidate,output):
    began=time.monotonic();deadline=began+50
    candidate,output=candidate.resolve(),output.resolve()
    if candidate.is_dir():candidate=candidate/'construction.json'
    assert candidate.is_relative_to(WORK) and output.is_relative_to(WORK) and not output.exists()
    construction=json.loads(candidate.read_text())
    source_path=WORK/'source-v1/source-audit.json';source=json.loads(source_path.read_text())
    catalogue_path=Path(source['sourceCatalogues']['male']['path']);catalogue_sha=digest_file(catalogue_path)
    assert catalogue_sha==source['sourceCatalogues']['male']['sha256']==construction['source_catalogue_sha256']
    atlas=json.loads(catalogue_path.read_text());by_id={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    bones=[]
    for row in source['contextParts']:
        if row['sourcePart'].get('system')!='skeletal':continue
        assert by_id[row['partId']]==row['sourcePart'] and digest_file(row['bankPath'])==row['bankSha256']
        mesh=reader.read_part(atlas,by_id[row['partId']],frame=FRAME,catalogue_path=catalogue_path,catalogue_sha=catalogue_sha)
        assert mesh.geometry_sha==row['geometrySha256']
        bones.append((row,mesh,BVH(mesh)))
    assert len(bones)==8 and construction['meshes']
    assert len({row['key'] for row in construction['meshes']})==len(construction['meshes'])
    output.mkdir(parents=True);results=[]
    for row in construction['meshes']:
        assert digest_file(row['mesh'])==row['sha256']
        roi=Mesh.load(row['mesh']);assert roi.geometry_sha==row['geometry_sha256'] and roi.frame==FRAME
        roi_bvh=BVH(roi)
        for bone_row,bone,bone_bvh in bones:
            if time.monotonic()>=deadline:raise TimeoutError('Tongue bone audit reached its 50 s internal budget; saved pair stages remain incomplete')
            stage=time.monotonic()
            arrays,summary=pair_contacts(roi,bone,0.,first_bvh=roi_bvh,second_bvh=bone_bvh,max_pairs=50000,max_tests=1000000)
            bad=np.isin(arrays['relation'],['crossing','coplanar-overlap']);pairs=arrays['pairs'][bad]
            stem=row['key']+'--'+bone_row['partId'].replace(':','_').replace('/','_').replace(' ','_')
            witnesses=output/(stem+'.npz')
            atomic_npz(witnesses,**arrays,defect_pairs=pairs,defect_first_triangles=roi.triangles[pairs[:,0]],defect_second_triangles=bone.triangles[pairs[:,1]])
            def bounds(triangles):
                return None if not len(triangles) else [triangles.reshape(-1,3).min(0).tolist(),triangles.reshape(-1,3).max(0).tolist()]
            examples=[dict(roi_face=int(a),bone_face=int(b),relation=str(relation),roi_triangle=roi.triangles[a].tolist(),bone_triangle=bone.triangles[b].tolist()) for (a,b),relation in zip(pairs[:20],arrays['relation'][bad][:20])]
            result=dict(roi=row['key'],roi_path=row['mesh'],roi_sha256=row['sha256'],roi_geometry_sha256=roi.geometry_sha,bone=bone_row['partId'],bone_geometry_sha256=bone.geometry_sha,bone_bank_path=bone.source['bank_path'],bone_bank_sha256=bone.source['bank_sha256'],summary=summary,
                all_witnesses_npz=str(witnesses),all_witnesses_sha256=digest_file(witnesses),defect_pairs=len(pairs),roi_patch_bounds_m=bounds(roi.triangles[pairs[:,0]]),bone_patch_bounds_m=bounds(bone.triangles[pairs[:,1]]),first_twenty_witnesses=examples,elapsed_seconds=time.monotonic()-stage)
            atomic_json(output/(stem+'.json'),result);results.append(result)
            print(json.dumps(dict(roi=row['key'],bone=bone_row['partId'],crossing=summary['crossing_face_pairs'],coplanar=summary['coplanar_overlap_face_pairs'],tests=summary['candidate_triangle_pairs_tested'],complete=summary['complete_within_pair_budget'],seconds=result['elapsed_seconds'])),flush=True)
    expected_pairs=len(construction['meshes'])*len(bones)
    complete=len(results)==expected_pairs and all(r['summary']['complete_within_pair_budget'] for r in results)
    report=dict(schema='tongue-roi-native-bone-audit/v1',construction_path=str(candidate),construction_sha256=digest_file(candidate),source_audit_path=str(source_path),source_audit_sha256=digest_file(source_path),source_catalogue_path=str(catalogue_path),source_catalogue_sha256=catalogue_sha,frame=FRAME,pair_budget=dict(max_pairs=50000,max_tests=1000000),expected_pairs=expected_pairs,pairs=results,complete=complete,
        crossing_pairs=sum(r['summary']['crossing_face_pairs'] for r in results),coplanar_pairs=sum(r['summary']['coplanar_overlap_face_pairs'] for r in results),passed=complete and not any(r['defect_pairs'] for r in results),gate_relaxed=False,geometry_modified=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    assert digest_file(catalogue_path)==catalogue_sha
    atomic_json(output/'summary.json',report)
    print(json.dumps({k:report[k] for k in ['expected_pairs','complete','crossing_pairs','coplanar_pairs','passed','elapsed_seconds']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--candidate',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();main(args.candidate,args.output)
