"""Resolve only the five saved adaptive-depth uncertainties by finite segment/triangle distance."""
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.distance import triangle_pairs

def outside(a,b,center,r):
    u=b-a;A=float(u@u);v=a-center;B=float(2*v@u);C=float(v@v-r*r);disc=B*B-4*A*C
    if disc<=0:return [(0.,1.)]if np.linalg.norm((a+b)/2-center)>=r else []
    roots=sorted([(-B-np.sqrt(disc))/(2*A),(-B+np.sqrt(disc))/(2*A)]);lo=max(0.,roots[0]);hi=min(1.,roots[1])
    if lo>=hi:return [(0.,1.)]if np.linalg.norm((a+b)/2-center)>=r else []
    return ([(0.,lo)]if lo>0 else [])+([(hi,1.)]if hi<1 else [])

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'global-palatine-actual-radius-v39';E=O/'edge-search-evidence-v40.json';e=json.loads(E.read_bytes());c=json.loads((O/'construction-contract.json').read_bytes());q=np.load(D/'global-palatine-graph-v33/clearance-domain.npz')['points'];z=np.load(O/'literal-surface-query-bank-v40.npz');assert digest_file(O/'literal-surface-query-bank-v40.npz')==e['query_bank_sha256'];start=np.array(c['unchanged_start_m']);end=np.array(c['unchanged_recipient']['point_m']);threshold=c['radius_m']+c['margin_m'];frame='model:local-ta98:assembled-metres'
    parent=Mesh.load(D.parent/'male-left-facial-cervical-routes-anatomical-v5/facial-left.npz');pp=Mesh.load(D.parent/'male-left-facial-cervical-routes-anatomical-v5/palatopharyngeus-left.npz');rows=[]
    uncertain=[r for r in e['nontrivial_edge_certificates']if r['reason']=='subinterval-not-certified-at-depth-limit'];assert len(uncertain)==5
    for r in uncertain:
        a,b=q[r['nodes']];checks=[]
        for name,v,f,intervals,own,orig in [
          ('composite',z['vertices'],z['faces'],[(0.,1.)],z['source_index'],z['native_face']),
          ('ZA:Facial artery.l',parent.vertices,parent.faces,outside(a,b,start,.003),None,None),
          ('ZA:Palatopharyngeus muscle.l',pp.vertices,pp.faces,outside(a,b,end,.002),None,None)]:
            tr=v[f]
            for lo,hi in intervals:
                first=a+lo*(b-a);last=a+hi*(b-a);low=np.minimum(first,last)-threshold;high=np.maximum(first,last)+threshold;selected=np.flatnonzero(np.all(tr.max(1)>=low,axis=1)&np.all(tr.min(1)<=high,axis=1))
                if not len(selected):checks.append(dict(surface=name,interval=[lo,hi],passed=True,no_literal_triangle_AABB_within_full_clearance=True));continue
                seg=np.array([first,last,first]);result=triangle_pairs(np.broadcast_to(seg,(len(selected),3,3)),tr[selected]);k=int(np.argmin(result['distance_m']));fi=int(selected[k]);bound=float(result['distance_m'][k]-result['numerical_roundoff_m'][k]);pid=e['source_bindings'][int(own[fi])]['part_id']if own is not None else name;sourceface=int(orig[fi])if orig is not None else fi
                checks.append(dict(surface=name,interval=[lo,hi],passed=bound>=threshold,source_id=pid,source_face=sourceface,segment_point_m=result['first'][k].tolist(),source_point_m=result['second'][k].tolist(),distance_m=float(result['distance_m'][k]),distance_lower_bound_m=bound,required_clearance_m=threshold,triangle_pairs_tested=len(selected)))
        rows.append(dict(nodes=r['nodes'],first_m=a.tolist(),last_m=b.tolist(),passed=all(x['passed']for x in checks),checks=checks,source_occupancy_basis='Both endpoints satisfy the179-source valid closed-source field. Full radius-clear surface segments cannot cross into these closed sources. Parent/PP outside-interface portions explicitly checked. Open-source occupied volumes remain unresolved.'))
    result=dict(schema='global-palatine-uncertain-edge-analytic-resolution/v41',source_catalogue_sha256=c['source_catalogue_sha256'],prior_search_sha256=digest_file(E),query_bank_sha256=e['query_bank_sha256'],rows=rows,passed_edges=sum(r['passed']for r in rows),failed_edges=sum(not r['passed']for r in rows),no_new_route_search_or_geometry=True,no_anatomical_acceptance_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'uncertain-edge-resolutions-v41.json',result);print(json.dumps(dict(passed_edges=result['passed_edges'],failed_edges=result['failed_edges'],rows=[dict(nodes=r['nodes'],passed=r['passed'],limiting=[dict(source=x.get('source_id'),distance_m=x.get('distance_m'))for x in r['checks']if not x['passed']])for r in rows],elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
