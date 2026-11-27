"""Exact source intersection curves before any numeric anterior-circle fit."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
import sys,time,json,heapq
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import winding_contains
from geometry_tools.algorithms.exact_source_intersections import edge_plane,contains

CAT=ROOT/'public/models/atlas.json';BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1';OUT=BASE/'anterior-donor-counterpart-curves-v10'
FRAME='source:BodyParts3D-4.0:packed-native-metres-y-up'

def exact_segment(a,b):
    aa,bb=[tuple(tuple(Q.from_float(float(x))for x in p)for p in t)for t in [a,b]];points=set()
    for first,second in [(aa,bb),(bb,aa)]:
        for p in first:
            if contains(p,second):points.add(p)
        for i in range(3):
            q=edge_plane(first[i],first[(i+1)%3],second)
            if q is not None and contains(q,first)and contains(q,second):points.add(q)
    if len(points)<2:return None
    pp=list(points);axis=max(range(3),key=lambda k:max(p[k]for p in pp)-min(p[k]for p in pp));ends=min(pp,key=lambda q:q[axis]),max(pp,key=lambda q:q[axis])
    return ends if ends[0]!=ends[1]else None

def curves(first,second):
    d,summary=pair_contacts(first,second,0.,max_pairs=250000,max_tests=2000000);assert summary['complete_within_pair_budget'];assert not summary['coplanar_overlap_face_pairs']
    edges={};nodefaces=defaultdict(set)
    for af,bf in d['pairs']:
        seg=exact_segment(first.triangles[af],second.triangles[bf])
        if seg is None:continue
        key=tuple(sorted(seg));edges.setdefault(key,[]).append([int(af),int(bf)])
        for p in seg:nodefaces[p].add(int(af))
    graph=defaultdict(set)
    for p,q in edges:graph[p].add(q);graph[q].add(p)
    unseen=set(graph);components=[]
    while unseen:
        start=min(unseen);stack=[start];group=set()
        while stack:
            p=stack.pop()
            if p in group:continue
            group.add(p);unseen.discard(p);stack.extend(graph[p]-group)
        closed=all(len(graph[p])==2 for p in group);ordered=[]
        if closed:
            p=start;previous=None
            while not ordered or p!=start:
                ordered.append(p);options=graph[p]-({previous}if previous is not None else set());q=min(options);previous,p=p,q
                assert len(ordered)<=len(group)
            assert len(ordered)==len(group)
        else:ordered=sorted(group)
        xyz=np.array([[float(x)for x in p]for p in ordered]);lengths=np.linalg.norm(np.roll(xyz,-1,axis=0)-xyz,axis=1)if closed else np.array([])
        centroid=(lengths[:,None]*(xyz+np.roll(xyz,-1,axis=0))/2).sum(0)/lengths.sum()if closed else xyz.mean(0)
        _,_,vt=np.linalg.svd(xyz-centroid,full_matrices=False);normal=vt[-1];plane=(xyz-centroid)@normal
        components.append(dict(closed=closed,nodes=len(group),sourceCurvePoints_m=xyz.tolist(),exactRationalKeys=[[str(x)for x in p]for p in ordered],arcLengthCentroid_m=centroid.tolist(),perimeter_m=float(lengths.sum())if closed else None,diagnosticBestFitPlaneNormal=normal.tolist(),planeMaximumResidual_m=float(abs(plane).max()),planeRmsResidual_m=float(np.sqrt(np.mean(plane*plane))),actualSourceFacePairs=[pair for edge,owners in edges.items()if edge[0]in group for pair in owners],distanceWelding=False))
    return dict(summary=summary,components=components,degreeHistogram={str(k):sum(len(graph[p])==k for p in graph)for k in sorted(set(map(len,graph.values())))},rationalEndpointConnectivity=True,sourceMouthOrLumenInferred=False)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();loaded={}
    def get(pid):
        if pid not in loaded:loaded[pid]=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT)
        return loaded[pid]
    rows=[]
    for side,aca,ica,calloso in [('l','FJ1654M','FJ1682M','FJ1676M'),('r','FJ1654','FJ1682','FJ1676')]:
        root=curves(get(aca),get(ica));distal=curves(get(aca),get(calloso));m=get(aca)
        rootfaces=np.unique(np.array([p for c in root['components']for p in c['actualSourceFacePairs']],int)[:,0]);anchors=np.unique(m.faces[rootfaces]);graph=defaultdict(dict)
        for tri in m.faces:
            for x,y in zip(tri,np.roll(tri,-1)):
                w=float(np.linalg.norm(m.vertices[x]-m.vertices[y]));graph[int(x)][int(y)]=w;graph[int(y)][int(x)]=w
        distances=np.full(len(m.vertices),np.inf);distances[anchors]=0.;pending=[(0.,int(x))for x in anchors]
        while pending:
            d,x=heapq.heappop(pending)
            if d!=distances[x]:continue
            for y,w in graph[x].items():
                n=d+w
                if n<distances[y]:distances[y]=n;heapq.heappush(pending,(n,y))
        far=np.flatnonzero(distances>=np.quantile(distances,.99));where=winding_contains(get(calloso),m.vertices[far]);inside=all(where['inside'])and not any(where['boundary_ambiguous'])
        ready=len(root['components'])==1 and root['components'][0]['closed']and len(distal['components'])==1 and distal['components'][0]['closed']and inside
        rows.append(dict(side=side,donorAca=aca,donorIca=ica,donorDistalContinuation=calloso,sourceRootIntersection=root,sourceDistalContinuationIntersection=distal,geodesicFarEndpointSourceVertexIds=far.tolist(),geodesicFarEndpointSourcePoints_m=m.vertices[far].tolist(),geodesicFarEndpointInsideActualSourceCallosomarginalCast=where,farEndpointWithinDistalContinuation=inside,counterpartCurvesEstablished=ready,remainingLimitation='The source curves bound intersecting outer casts, not paired wall annuli or an established fluid mouth. A distal junction curve is accepted for a numeric trial only when it encloses the actual far end of this named ACA source, using root-relative source graph geodesics rather than axis or nearest-face guesses.'))
    source_records={pid:dict(record=parts[pid],geometrySha256=m.geometry_sha,bankSha256=m.source['bank_sha256'],attributeSha256=m.attributes_sha)for pid,m in loaded.items()}
    proof=dict(schema='anterior-source-counterpart-curves/v10',sourceCatalogueSha256=digest_file(CAT),sourceFrame=FRAME,rows=rows,sourceRecords=source_records,allFourCounterpartCurvesEstablished=all(r['counterpartCurvesEstablished']for r in rows),numericFitRun=False,noGeometryOrLiveMetadataChanged=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',proof)
    print(json.dumps(dict(output=str(OUT/'proof.json'),curvesReady=proof['allFourCounterpartCurvesEstablished'],sides=[dict(side=r['side'],rootCurves=[dict(closed=c['closed'],nodes=c['nodes'])for c in r['sourceRootIntersection']['components']],distalCurves=[dict(closed=c['closed'],nodes=c['nodes'])for c in r['sourceDistalContinuationIntersection']['components']],farEndpointInsideActualDistalCast=r['farEndpointWithinDistalContinuation'])for r in rows],elapsedSeconds=proof['elapsedSeconds'])))
if __name__=='__main__':main()
