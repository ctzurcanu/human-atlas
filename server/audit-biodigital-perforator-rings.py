"""Follow exact source cylindrical-strip topology; compare matched wall rings."""
import json
import time
from pathlib import Path
import numpy as np

W = Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')


def ordered_cycle(vertices, edges):
    graph = {int(v): [] for v in vertices}
    for a, b in edges:
        if a in graph and b in graph:
            graph[int(a)].append(int(b)); graph[int(b)].append(int(a))
    assert all(len(x) == 2 for x in graph.values()), 'Ring must be a literal closed edge cycle'
    first = min(graph); order = [first]; previous, current = first, min(graph[first])
    while current != first:
        order.append(current)
        previous, current = current, next(v for v in graph[current] if v != previous)
        assert len(order) <= len(graph)
    assert len(order) == len(graph)
    return np.asarray(order)


def rings(vertices, faces):
    parent = np.arange(len(vertices))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a, b, c in faces:
        parent[find(a)] = find(b); parent[find(c)] = find(b)
    labels = np.asarray([find(i) for i in range(len(vertices))])
    all_edges = np.sort(np.concatenate([faces[:, [0,1]], faces[:, [1,2]], faces[:, [2,0]]]), axis=1)
    unique, counts = np.unique(all_edges, axis=0, return_counts=True)
    adjacency = {i:set() for i in range(len(vertices))}
    for a, b in unique:
        adjacency[int(a)].add(int(b)); adjacency[int(b)].add(int(a))
    result=[]
    for label in np.unique(labels):
        members=set(np.flatnonzero(labels == label).tolist())
        boundary=unique[(counts == 1)&np.isin(unique[:,0],list(members))]
        assert len(boundary) > 0
        bg={}
        for a,b in boundary:
            bg.setdefault(int(a),set()).add(int(b));bg.setdefault(int(b),set()).add(int(a))
        assert all(len(v)==2 for v in bg.values())
        unseen=set(bg);loops=[]
        while unseen:
            first=min(unseen); connected={first}; todo=[first]
            while todo:
                a=todo.pop()
                for b in bg[a]-connected:connected.add(b);todo.append(b)
            unseen-=connected;loops.append(connected)
        assert len(loops)==2, 'Two actual mouths required for a strip'
        current=min(loops,key=min);level=[];used=set()
        while current:
            order=ordered_cycle(current,unique)
            level.append(order);used |= current
            current=set.union(*(adjacency[i] for i in current))-used
            if current:assert len(current)==len(order), 'Branching or non-strip layer'
        assert used == members
        lookup={int(v):i for i,ring in enumerate(level)for v in ring}
        tf=faces[np.isin(faces[:,0],list(members))]
        layers=np.array([[lookup[int(v)]for v in face]for face in tf])
        assert np.all(np.ptp(layers,axis=1)==1), 'All triangles must join adjacent rings'
        result.append(level)
    return result


def metric(p):
    c=p.mean(0); vector=np.cross(p-c,np.roll(p,-1,axis=0)-c).sum(0)/2
    area=np.linalg.norm(vector); normal=vector/area
    return dict(area_m2=float(area), centroid_m=c.tolist(),
        maximumPlaneResidual_m=float(np.max(np.abs((p-c)@normal))),
        diameter_m=float(np.linalg.norm(p[:,None]-p[None],axis=2).max()))


def main():
    started=time.perf_counter();out=W/'saphenous-transport-v2/actual-perforator-ring-audit-v1';out.mkdir(exist_ok=False)
    ownership=json.loads((W/'lower-limb-joints-v4/candidate-v1/instance-ownership.json').read_text())['instances']
    source=[]
    for stage in ['lower-limb-joints-v4','saphenous-transport-v2']:
        with np.load(W/stage/'candidate-v1/venous.npz')as z:source.append((z['vertices'],z['faces']))
    rows=[]
    for part in ownership:
        if part['name'] not in ['Left femoral canal perforator veins','Right femoral canal perforator veins']:continue
        start=part['candidateVertexStart'];stop=start+part['candidateVertexCount']
        bv,af=source[0];av,bf=source[1]
        f=af[part['candidateFaceStart']:part['candidateFaceStart']+part['candidateFaceCount']]-start
        assert np.array_equal(f,bf[part['candidateFaceStart']:part['candidateFaceStart']+part['candidateFaceCount']]-start)
        groups=rings(bv[start:stop],f)
        for ci,group in enumerate(groups):
            profile=[]
            for ri,ring in enumerate(group):
                a=metric(bv[start+ring]);b=metric(av[start+ring])
                profile.append(dict(ring=ri,sourceVertexIds=ring.tolist(),before=a,after=b,
                    areaRatio=b['area_m2']/a['area_m2'],diameterRatio=b['diameter_m']/a['diameter_m']))
            rows.append(dict(sourceId=part['sourceInstanceId'],component=ci,rings=profile,
                minimumAreaRatio=min(p['areaRatio']for p in profile),minimumDiameterRatio=min(p['diameterRatio']for p in profile),
                topology='Two literal source mouths; all complete edge rings and every strip triangle accounted for. No loft, smoothing or synthetic section.'))
    report=dict(schema='exact-matched-source-strip-ring-audit/v1',rows=rows,geometryModified=False,anatomicalAcceptance=False,
        minimumRingAreaRatio=min(r['minimumAreaRatio']for r in rows),elapsedSeconds=time.perf_counter()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['minimumRingAreaRatio','elapsedSeconds']}))
    print(json.dumps([dict(source=r['sourceId'],component=r['component'],areaRatio=r['minimumAreaRatio'],diameterRatio=r['minimumDiameterRatio'])for r in rows]))


if __name__=='__main__':main()
