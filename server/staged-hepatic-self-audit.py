"""Save and replay disjoint, bounded BVH self-audit stages for a large wall."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.contact import self_contact_summary
from geometry_tools.spatial import BVH
ap=argparse.ArgumentParser();ap.add_argument('--version',required=True);ap.add_argument('--role',default='hepatic_vein');ap.add_argument('--stage',choices=['plan','check','collect','diagnose'],required=True);ap.add_argument('--piece',type=int);ap.add_argument('--partition-dir',default='');ap.add_argument('--depth',type=int,default=1);args=ap.parse_args();began=time.monotonic()
source_folder=ROOT/'.local-models/geometry-workbench/liver-completion'/('hepatic-contained-banks-'+args.version);record=json.loads((source_folder/(args.role+'-pending.json')).read_text());layout=record['layout'];raw=Path(record['bank']).read_bytes();assert digest_file(record['bank'])==record['bankSha256'];n=layout['vertexCount'];folder=source_folder/args.partition_dir;folder.mkdir(exist_ok=True)
v=np.ndarray((n,3),'<f4',buffer=raw,offset=layout['positions']).astype(float)+layout['geometryOrigin'];f=np.ndarray((layout['indexCount']//3,3),'<u4',buffer=raw,offset=layout['indices']);attributes={k:np.ndarray((n,3)if k=='normals'else(n,),dtype,buffer=raw,offset=layout[k])for k,dtype in [('normals','<i2'),('normalSigns','i1')]};mesh=Mesh(v,f,'model:local-ta98:assembled-metres',attributes)
planpath=folder/(args.role+'-audit-plan.json');bvhpath=folder/(args.role+'-audit-bvh.npz');topopath=folder/(args.role+'-topology.json')
implementations={name:digest_file(ROOT/'geometry_tools'/name)for name in ['audit.py','contact.py','distance.py','spatial.py','model.py']}
if args.stage=='plan':
    topo=audit(mesh,self_intersections=False,require_closed=True);atomic_json(topopath,topo)
    if not topo['passed']:raise ValueError('Actual encoded topology failed')
    tree=BVH(mesh);tree.save(bvhpath);frontier=[0];assert 1<=args.depth<=4
    for depth in range(args.depth):frontier=[int(c)for node in frontier for c in tree.children[node]];assert min(frontier)>=0
    roots=[[x,y]for i,x in enumerate(frontier)for y in frontier[i:]];atomic_json(planpath,dict(bankSha256=record['bankSha256'],geometrySha256=mesh.geometry_sha,attributesSha256=mesh.attributes_sha,implementations=implementations,bvhSha256=digest_file(bvhpath),rootPairs=roots,allPairsCoveredByDisjointBinaryRootPartition=True));print(json.dumps(dict(stage='plan',pieces=len(roots),triangles=len(f),seconds=time.monotonic()-began)));sys.exit(0)
plan=json.loads(planpath.read_text());assert plan['bankSha256']==record['bankSha256']and plan['geometrySha256']==mesh.geometry_sha and plan['attributesSha256']==mesh.attributes_sha
if args.stage!='diagnose':assert plan['implementations']==implementations
if args.stage in ['check','diagnose']:
    assert args.piece in range(len(plan['rootPairs']))and digest_file(bvhpath)==plan['bvhSha256'];tree=BVH.load(mesh,bvhpath);summary,_=self_contact_summary(mesh,bvh=tree,root_pair=plan['rootPairs'][args.piece],max_pairs=3000000,max_tests=24000000,stop_on_defect=args.stage=='check',example_limit=10000 if args.stage=='diagnose'else 20)
    result=dict(planSha256=digest_file(planpath),piece=args.piece,rootPair=plan['rootPairs'][args.piece],summary=summary,seconds=time.monotonic()-began,passed=summary['complete_within_pair_budget']and summary['defect_pairs']==0)
    atomic_json(folder/(args.role+'-self-'+('all-'if args.stage=='diagnose'else'piece-')+str(args.piece)+'.json'),result);print(json.dumps({**result,'summary':{k:v for k,v in summary.items()if k!='examples'}}),flush=True)
    if not result['passed']:raise ValueError('A staged self-audit found a defect or exhausted its budget')
    sys.exit(0)
pieces=[json.loads((folder/(args.role+'-self-piece-'+str(i)+'.json')).read_text())for i in range(len(plan['rootPairs']))]
assert all(p['passed']and p['planSha256']==digest_file(planpath)and p['piece']==i for i,p in enumerate(pieces))
summary=dict(distance_threshold_m=0.,truncated=False,complete_within_pair_budget=True,defect_pairs=0,examples=[],partitionPlan=str(planpath),partitionPlanSha256=digest_file(planpath))
for key in ['candidate_triangle_pairs_tested','returned_triangle_pairs','crossing_face_pairs','coplanar_overlap_face_pairs','touching_face_pairs','degenerate_pair_count']:summary[key]=sum(p['summary'][key]for p in pieces)
assert summary['candidate_triangle_pairs_tested']<=48000000 and summary['returned_triangle_pairs']<=3000000
result=json.loads(topopath.read_text());result.update(self_intersections={**summary,'triangle_test_budget':48000000,'pair_result_budget':3000000},passed=True,failures=[]);atomic_json(source_folder/(args.role+'-audit.json'),result)
record.update(passed=True,geometrySha256=mesh.geometry_sha,partitionedAudit=str(planpath));atomic_json(source_folder/(args.role+'-evidence.json'),record);print(json.dumps(dict(passed=True,triangles=len(f),tests=summary['candidate_triangle_pairs_tested'])))
