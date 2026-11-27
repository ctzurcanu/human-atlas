"""Independent saved checks of the actual emitted six-part female revision."""
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'server')]
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH
from tongue_basal_contract import require
FRAME='model:local-female-ta98:assembled-metres'
BONES=['LOCAL:female:mandible:960','LOCAL:female:hyoid:1002','LOCAL:female:skull_cranium:999','LOCAL:female:skull_cranium:1000','LOCAL:female:skull_cranium:1001']
SOFT=['LOCAL:female:mylohyoid:238','LOCAL:female:l_geniohyoid:241','LOCAL:female:r_geniohyoid:254','LOCAL:female:uvula:1088','LOCAL:female:pharynx:1087','LOCAL:female:nasopharynx:1047']

def main(candidate,stage,output):
    began=time.monotonic();candidate=candidate.resolve();require(candidate.is_file() and not output.exists(),'Explicit candidate and new audit stage required')
    atlas=json.loads(candidate.read_text());evidence_path=candidate.with_name('transfer-evidence.json');evidence=json.loads(evidence_path.read_text())
    require(evidence['candidate']['sha256']==digest_file(candidate) and evidence['schema']=='tongue-female-core-native-partial/v1','Actual partial candidate binding changed')
    new=evidence['new_part_ids'];require(len(set(new))==6,'Only six new carriers may enter this audit')
    parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();meshes={pid:reader.read_part(atlas,parts[pid],frame=FRAME) for pid in new}
    rows=[];output.mkdir(parents=True)
    if stage=='self':
        for pid in new:
            started=time.monotonic();m=meshes[pid];result=audit(m,max_pairs=250000,max_tests=2000000,require_closed=True)
            t=m.triangles;raw=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);norm=m.attributes['normals'][m.faces].astype(float)*m.attributes['normalSigns'][m.faces,None]
            bad=int(np.count_nonzero(np.einsum('fci,fi->fc',norm,raw)<=0))
            rows.append(dict(partId=pid,geometry_sha256=m.geometry_sha,result=result,nonpositive_corner_dots=bad,
                             passed=result['passed'] and bad==0,elapsed_seconds=time.monotonic()-started))
        passed=all(r['passed'] for r in rows)
        report=dict(schema='female-partial-tongue-actual-self-audit/v1',checks=rows,passed=passed)
    else:
        context=BONES if stage=='bones' else SOFT
        fixed={pid:reader.read_part(atlas,parts[pid],frame=FRAME) for pid in context}
        first_trees={pid:BVH(m) for pid,m in meshes.items()};second_trees={pid:BVH(m) for pid,m in fixed.items()}
        for i,first in enumerate(new):
            for j,second in enumerate(context):
                data,summary=pair_contacts(meshes[first],fixed[second],0.,first_bvh=first_trees[first],second_bvh=second_trees[second],max_pairs=250000,max_tests=2000000)
                path=output/(str(i).zfill(2)+'--'+str(j).zfill(2)+'.npz');atomic_npz(path,**data)
                rows.append(dict(partId=first,targetPartId=second,first_geometry_sha256=meshes[first].geometry_sha,
                    second_geometry_sha256=fixed[second].geometry_sha,
                    constraint=dict(kind='disjoint',parameters=dict(distance_m=0.,max_pairs=250000,max_tests=2000000,
                        scope='Actual complete finite triangle crossings/coplanar area; solid containment independently rechecked by install job')),
                    summary=summary,witnesses=str(path.resolve()),witnesses_sha256=digest_file(path)))
        complete=all(r['summary']['complete_within_pair_budget'] and not r['summary']['truncated'] for r in rows)
        crossed=sum(r['summary']['crossing_face_pairs'] for r in rows);coplanar=sum(r['summary']['coplanar_overlap_face_pairs'] for r in rows)
        failed=[dict(partId=r['partId'],targetPartId=r['targetPartId'],crossings=r['summary']['crossing_face_pairs'],coplanar=r['summary']['coplanar_overlap_face_pairs']) for r in rows if r['summary']['crossing_face_pairs'] or r['summary']['coplanar_overlap_face_pairs'] or r['summary']['degenerate_pair_count']]
        passed=complete and not failed
        report=dict(schema='female-partial-tongue-context-audit/v1',group=stage,pairs=rows,all_pairs_complete=complete,
                    proper_crossing_pairs=crossed,coplanar_overlap_pairs=coplanar,failed_pairs=failed,passed=passed)
    report.update(candidate=str(candidate),candidate_sha256=digest_file(candidate),transfer_evidence_sha256=digest_file(evidence_path),
        fixed_skeleton=True,source_geometry_modified=False,retained_native_muscles_audited_as_unchanged_baseline=False,
        geometry_fit_accepted=False,qualifying_candidate=False,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    path=output/'evidence.json';atomic_json(path,report)
    print(json.dumps(dict(stage=stage,report=str(path.resolve()),passed=passed,count=len(rows),
                         crossings=report.get('proper_crossing_pairs'),failed=report.get('failed_pairs'),elapsed_seconds=report['elapsed_seconds'])))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--stage',choices=['self','bones','soft'],required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.candidate,a.stage,a.output.resolve())
