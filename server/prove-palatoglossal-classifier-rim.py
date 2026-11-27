"""Exact finite triangle proof for the declared classifier reflection rim.

No proximity/coplanarity tolerance is used. For every reported cap/M contact,
the complete exact triangle intersection is confined to literal declared rim
vertices or edges. Generic geometry predicates remain unchanged.
"""
import argparse,json,sys,time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature
from geometry_tools.algorithms.exact_source_intersections import contains,edge_plane,_sub,_dot,_cross
from palatoglossal_incidence_fixture import qvec,plane_coordinates,clip_polygon,oriented_area
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def read(p):return json.loads(Path(p).read_text())
def require(x,message):
    if not x:raise ValueError(message)
def in_segment(p,a,b):
    d=_sub(b,a);w=_sub(p,a)
    return not any(_cross(d,w))and 0<=_dot(w,d)<=_dot(d,d)
def intersection_vertices(a,b):
    na=_cross(_sub(a[1],a[0]),_sub(a[2],a[0]));nb=_cross(_sub(b[1],b[0]),_sub(b[2],b[0]))
    if any(_cross(na,nb)):
        points=set(a)&set(b)
        for t,other in [(a,b),(b,a)]:
            for i in range(3):
                hit=edge_plane(t[i],t[(i+1)%3],other)
                if hit is not None and contains(hit,other)and contains(hit,t):points.add(hit)
        return points,'nonparallel-exact-edge-plane'
    if _dot(na,_sub(b[0],a[0])):return set(),'parallel-disjoint'
    axes=plane_coordinates(a);drop=next(i for i in range(3)if i not in axes);ac=[tuple(p[i]for i in axes)for p in a];bc=[tuple(p[i]for i in axes)for p in b]
    common=clip_polygon(ac,bc);require(oriented_area(common)==0,'A cap/M contact has exact coplanar area')
    points=set()
    for p in common:
        q=[Q(0)]*3
        for i,x in zip(axes,p):q[i]=x
        q[drop]=(_dot(na,a[0])-sum(na[i]*q[i]for i in axes))/na[drop];points.add(tuple(q))
    return points,'coplanar-zero-area-exact-clip'
def main(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use a new exact rim proof directory')
    source_path=W/'palatoglossal-direct-inner-domain-v1/evidence.json';source=read(source_path);require(source['classifier_geometry_passed'],'Direct classifier geometry failed')
    plan=read(W/'palatoglossal-original-parent-v2/plan.json');mu=Mesh.load(plan['source_files']['mucosa']);cap_path=W/'palatoglossal-direct-inner-domain-v1/basal-classifier-only.npz';cap=Mesh.load(cap_path)
    outer=mu.faces[:348];edge=boundary_signature(outer);rim_edges={tuple(sorted((qvec(mu.vertices[a]),qvec(mu.vertices[b]))))for a,b in edge};rim_controls={p for e in rim_edges for p in e};require(len(rim_edges)==64,'Declared original reflection loop changed')
    row=next(r for r in source['basal_pair_checks']if r['host']=='mucosa');z=np.load(row['witnesses']);require(len(z['pairs'])==770,'Actual retained contact witness count changed')
    aq={};bq={};proof=[]
    for pair,relation in zip(z['pairs'],z['relation']):
        require(time.monotonic()-began<17,'Exact rim proof approached its finite deadline')
        i,j=map(int,pair);require(str(relation)=='touching','A cap/M pair is crossing or has coplanar area')
        a=aq.setdefault(i,tuple(qvec(p)for p in cap.triangles[i]));b=bq.setdefault(j,tuple(qvec(p)for p in mu.triangles[j]));common=set(a)&set(b)
        require(common and common<=rim_controls,'Pair has no literal declared reflection control')
        exact,method=intersection_vertices(a,b);require(exact,'Reported contact lost its literal exact intersection')
        if len(common)==1:require(exact==common,'Exact cap intersection extends beyond the literal reflection vertex')
        elif len(common)==2:
            ca,cb=sorted(common);require((ca,cb)in rim_edges,'Shared controls are not one declared reflection edge');require(all(in_segment(p,ca,cb)for p in exact),'Exact intersection leaves its declared reflection edge')
        else:raise ValueError('A cap/M contact has at least3common controls')
        proof.append(dict(cap_face=i,mucosa_face=j,method=method,declared_shared_controls=len(common),exact_intersection_vertices=len(exact),entire_intersection_on_declared_rim=True))
    arch_rows=[r for r in source['basal_pair_checks']if r['host'].startswith('arch')];require(all(r['summary']['complete_within_pair_budget']and r['summary']['returned_triangle_pairs']==0 for r in arch_rows),'Basal classifier touches an arch host')
    output.mkdir(parents=True);report=dict(schema='PG-exact-classifier-reflection-rim-proof/v1',passed=True,all770pairs_proved=True,actual_touching_pairs=len(proof),declared_reflection_edges=len(rim_edges),declared_reflection_controls=len(rim_controls),
        proofs=proof,off_rim_points_or_segments=0,proper_crossings=0,coplanar_overlap_area=0,basal_arch_contacts=0,generic_predicates_changed=False,generic_tolerances_changed=False,
        intentional_domain_ancestry_not_tissue_overlap_exception=True,source_bindings={str(p):digest_file(p)for p in [source_path,Path(row['witnesses']),cap_path,Path(plan['source_files']['mucosa']),Path(__file__).resolve()]},
        anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'evidence.json',report);print(json.dumps({k:v for k,v in report.items()if k not in ['proofs','source_bindings']}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
