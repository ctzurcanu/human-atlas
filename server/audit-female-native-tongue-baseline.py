"""Literal native female tongue baseline in the exact transfer context scope.

Open anatomy is never capped or interpreted as an occupied solid. Saved
triangle witnesses distinguish pre-existing source crossings from new fits.
"""
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import audit
from tongue_whole_input import WORK

PARTS=['LOCAL:female:tongue:1089','LOCAL:female:genioglossus:239',
       'LOCAL:female:l_hyoglossus:242','LOCAL:female:r_hyoglossus:255',
       'LOCAL:female:l_styloglossus:243','LOCAL:female:r_styloglossus:256']
BONES=['LOCAL:female:mandible:960','LOCAL:female:hyoid:1002','LOCAL:female:skull_cranium:999',
       'LOCAL:female:skull_cranium:1000','LOCAL:female:skull_cranium:1001']
SOFT=['LOCAL:female:mylohyoid:238','LOCAL:female:l_geniohyoid:241','LOCAL:female:r_geniohyoid:254',
      'LOCAL:female:uvula:1088','LOCAL:female:pharynx:1087','LOCAL:female:nasopharynx:1047']
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
EXPECTED='770f47f802c79cf6f3e4336bf2fa7c42cf710998069221f533a1ea3f8322e470'
FRAME='model:local-female-ta98:assembled-metres'

def bound(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=digest_file(p))

def main(output):
    began=time.monotonic()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new source-bound stage')
    if digest_file(CAT)!=EXPECTED:raise ValueError('Native source changed; do not audit a transferred catalogue as baseline')
    atlas=json.loads(CAT.read_text());parts={r['id']:r for r in atlas['parts']};reader=AtlasReader();output.mkdir()
    meshes={};trees={};sources=[]
    for i,pid in enumerate(PARTS+BONES+SOFT):
        mesh=reader.read_part(atlas,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED)
        path=output/f'source-{i:02d}.npz';mesh.save(path);meshes[pid]=mesh;trees[pid]=BVH(mesh)
        sources.append(dict(part_id=pid,name=parts[pid]['name'],mesh=bound(path),geometry_sha256=mesh.geometry_sha,
            bank_sha256=mesh.source['bank_sha256'],triangles=len(mesh.faces),fixed=True,
            role='tongue-assembly'if pid in PARTS else'bone'if pid in BONES else'floor'if pid in SOFT[:3]else'context'))
    audits=[]
    for pid in PARTS:
        m=meshes[pid];result=audit(m,require_closed=False,max_pairs=250000,max_tests=2000000)
        audits.append(dict(part_id=pid,checks=result,source_geometry_unmodified=True,
            solid_occupancy_established=False,no_cap_added=True))
    rows=[]
    for i,pid in enumerate(PARTS):
        for j,target in enumerate(BONES+SOFT):
            result,summary=pair_contacts(meshes[pid],meshes[target],0.,first_bvh=trees[pid],second_bvh=trees[target],max_pairs=250000,max_tests=2000000)
            path=output/f'pair-{i:02d}-{j:02d}.npz';atomic_npz(path,**result)
            rows.append(dict(part_id=pid,target_part_id=target,context_group='bones'if target in BONES else'floor'if target in SOFT[:3]else'other-soft',
                summary=summary,witnesses=bound(path),pre_existing=True,
                first_geometry_sha256=meshes[pid].geometry_sha,second_geometry_sha256=meshes[target].geometry_sha))
        atomic_json(output/f'progress-{i+1:02d}.json',dict(source_catalogue=bound(CAT),completed_pairs=len(rows),all_pairs=66,elapsed_seconds=time.monotonic()-began))
    def summarize(selected):return dict(pairs=len(selected),crossing_pairs=sum(r['summary']['crossing_face_pairs']for r in selected),
        coplanar_overlap_pairs=sum(r['summary']['coplanar_overlap_face_pairs']for r in selected),
        failed_pair_objects=sum(bool(r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs'])for r in selected),
        all_complete=all(r['summary']['complete_within_pair_budget']for r in selected))
    out=dict(schema='female-native-tongue-context-baseline/v1',source_catalogue=bound(CAT),frame=FRAME,
        sources=sources,source_part_audits=audits,pairs=rows,
        totals=summarize(rows),by_assembly_object={pid:summarize([r for r in rows if r['part_id']==pid])for pid in PARTS},
        by_context={pid:summarize([r for r in rows if r['target_part_id']==pid])for pid in BONES+SOFT},
        by_context_group={group:summarize([r for r in rows if r['context_group']==group])for group in ['bones','floor','other-soft']},
        scope='Same exact 5-bone/6-soft-object proper triangle-crossing scope as the 132-pair twelve-part transfer check; native has six parts and therefore 66 pairs. Counts do not measure overlap volume or establish attachment.',
        no_cap_added=True,no_fit_or_mutation=True,pre_existing_native_source_only=True,installed=False,anatomical_acceptance=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'baseline.json',out)
    render=[]
    for row in sources:
        if row['part_id']not in PARTS+BONES[:2]:continue
        color=[.76,.38,.40,.5]if row['part_id']==PARTS[0]else[.66,.27,.28,1]if row['part_id']in PARTS else[.78,.74,.64,.22]
        render.append(dict(mesh=row['mesh']['path'],name=row['name'],color=color))
    atomic_json(output/'review-spec.json',dict(frame=FRAME,resolution=768,wall_seconds=45,meshes=render,
        views=[dict(name='native-female-tongue-superior',direction=[0,1,0]),dict(name='native-female-extrinsic-right',direction=[1,.1,0]),
            dict(name='native-female-extrinsic-sagittal',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.547,.05],clear_outer=True))],
        qualification='Unchanged native source course guides. Open Tongue is not capped and not interpreted as an occupied solid; sparse native muscle envelopes are not accepted attachment reconstructions.'))
    print(json.dumps(dict(report=str(output/'baseline.json'),totals=out['totals'],by_context_group=out['by_context_group'],
        nonzero_pairs=[dict(part=r['part_id'],target=r['target_part_id'],crossings=r['summary']['crossing_face_pairs'],coplanar=r['summary']['coplanar_overlap_face_pairs'])for r in rows if r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs']],elapsed_seconds=time.monotonic()-began)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
