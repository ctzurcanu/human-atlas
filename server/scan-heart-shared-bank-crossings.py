"""Three disjoint, saved BVH stages cover the complete physical boundary bank."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import self_contact_summary
from geometry_tools.common import atomic_json
ap=argparse.ArgumentParser();ap.add_argument('--bank',required=True)
ap.add_argument('--prepare',action='store_true');ap.add_argument('--piece',type=int,choices=[0,1,2])
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.bank
e=json.loads((src/'emission.json').read_text());z=np.load(src/'arrangement.npz')
f=np.concatenate([z[key]for key in e['surfaces'].values()]);mesh=Mesh(z['positions'],f,e['frame'])
cache=src/'physical-bvh.npz'
if args.prepare:
    tree=BVH(mesh);tree.save(cache);left,right=map(int,tree.children[0])
    atomic_json(src/'physical-scan-plan.json',dict(geometry_sha256=mesh.geometry_sha,
        pairs=[[left,left],[left,right],[right,right]],
        complete_disjoint_bvh_cover=True,triangles=len(f)))
    print('Prepared complete physical boundary scan',len(f),'triangles',flush=True)
else:
    plan=json.loads((src/'physical-scan-plan.json').read_text());assert plan['geometry_sha256']==mesh.geometry_sha
    tree=BVH.load(mesh,cache);pair=plan['pairs'][args.piece]
    report,_=self_contact_summary(mesh,bvh=tree,root_pair=pair,max_pairs=4000000,
        max_tests=22000000,example_limit=4000)
    assert report['defect_pairs']==len(report['examples']), 'Complete defect list required'
    atomic_json(src/f'physical-scan-{args.piece}.json',dict(report,
        geometry_sha256=mesh.geometry_sha,root_pair=pair,elapsed_seconds=time.monotonic()-start))
    print(json.dumps(dict(piece=args.piece,complete=report['complete_within_pair_budget'],
        defects=report['defect_pairs'],tests=report['candidate_triangle_pairs_tested'],
        elapsed_seconds=time.monotonic()-start)),flush=True)
