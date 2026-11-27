"""Actual source transverse-handle cycles and qualified foramen-magnum contours.

Probe planes locate source holes; they are not anatomical region boundaries.
No source faces or positions are edited, repaired, capped, or installed.
"""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.algorithms.passage_cycles import handle_cycle
from geometry_tools.algorithms.simple_compartment_section import polygon_contains,polygon_clearance

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'

def literal_section(v,f,y):
    value=v[:,1]-y;nodes={};point=[];segments=[];coplanar=[]
    def node(a,b):
        if value[a]==0:key=('v',int(a));p=v[a]
        elif value[b]==0:key=('v',int(b));p=v[b]
        else:
            a,b=sorted((int(a),int(b)));key=('e',a,b);u=-value[a]/(value[b]-value[a]);p=v[a]+u*(v[b]-v[a])
        if key not in nodes:nodes[key]=len(point);point.append(p)
        return nodes[key]
    for i,ff in enumerate(f):
        if np.all(value[ff]==0):coplanar.append(i);continue
        nn=[]
        for a,b in zip(ff,np.roll(ff,-1)):
            if value[a]*value[b]<0 or(value[a]==0)!=(value[b]==0):nn.append(node(a,b))
        nn=list(dict.fromkeys(nn))
        if len(nn)==2:segments.append((nn[0],nn[1],i))
    point=np.array(point);adj={}
    for k,(a,b,_)in enumerate(segments):adj.setdefault(a,[]).append((b,k));adj.setdefault(b,[]).append((a,k))
    visited=set();rows=[];bad=[]
    for start in adj:
        if start in visited:continue
        todo=[start];component=[]
        while todo:
            q=todo.pop()
            if q in visited:continue
            visited.add(q);component.append(q);todo.extend(b for b,_ in adj[q])
        if any(len(adj[q])!=2 for q in component):bad.append(dict(nodes=component,degrees=[len(adj[q])for q in component]));continue
        ring=[];faces=[];last=None;cur=start
        while not ring or cur!=start:
            ring.append(cur);choices=[z for z in adj[cur]if z[0]!=last];n,k=choices[0];faces.append(segments[k][2]);last,cur=cur,n
            if len(ring)>len(component):raise ValueError('Invalid section walk')
        xyz=point[ring];xy=xyz[:,[0,2]];cross=xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1);a=.5*cross.sum()
        if abs(a)<1e-15:bad.append(dict(nodes=ring,reason='zero area'));continue
        center=((xy+np.roll(xy,-1,axis=0))*cross[:,None]).sum(0)/(6*a)
        rows.append(dict(worldPoints_m=xyz,planePoints_m=xy,sourceFaceIds=faces,area_m2=abs(a),center_m=[center[0],y,center[1]],centroidInside=polygon_contains(xy,center),centroidClearance_m=polygon_clearance(xy,center)))
    for i,r in enumerate(rows):
        r['containingContourIds']=[j for j,q in enumerate(rows)if i!=j and polygon_contains(q['planePoints_m'],r['planePoints_m'][0])]
        r['holeByNestedContourParity']=len(r['containingContourIds'])%2==1
    return rows,dict(noncycleSectionComponents=bad,coplanarSourceFaces=coplanar)

def main():
    began=time.monotonic();level=sys.argv[1];assert level in ['C1','C2','C4','C6','FM']
    out=BASE/('actual-aperture-discovery-independent-v1-'+level);out.mkdir(exist_ok=False)
    inv=json.loads((BASE/'source-inventory.json').read_bytes());bindings=inv['sourceBindings']
    name={'C1':'ZA:Atlas (C1)','C2':'ZA:Axis (C2)','C4':'ZA:Vertebra C4','C6':'ZA:Vertebra C6','FM':'ZA:Occipital bone'}[level]
    selected=[(pid,row)for pid,row in bindings.items()if pid.startswith(name+':')]
    arrays=[];sourceFaces=[]
    for pid,row in selected:
        assert digest_file(row['mesh'])==row['sha256'];z=np.load(row['mesh']);base=sum(len(v)for v,_ in arrays);arrays.append((z['vertices'],z['faces']+base));sourceFaces += [dict(partId=pid,sourceFace=i)for i in range(len(z['faces']))]
    vertices=np.vstack([v for v,_ in arrays]);faces=np.vstack([f for _,f in arrays]);v,ix=np.unique(vertices,axis=0,return_inverse=True);f=ix[faces]
    edges=np.sort(np.vstack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1);_,ct=np.unique(edges,axis=0,return_counts=True)
    topology=dict(physicalVertices=len(v),triangles=len(f),boundaryEdges=int(sum(ct==1)),nonmanifoldEdges=int(sum(ct>2)),materialPieces=[p for p,_ in selected],exactCoordinateWeldOnly=True)
    st=np.load(BASE/'l-literal-source-stations.npz');centers=st['centers'];probe={'C1':90,'C2':70,'C4':50,'C6':25,'FM':188}[level];y=float(centers[probe,1])
    contours,diagnostic=literal_section(v,f,y);holes=[(i,c)for i,c in enumerate(contours)if c['holeByNestedContourParity']and c['centroidInside']]
    result=[]
    for i,c in holes:
        if level!='FM'and abs(c['center_m'][0])<.008:continue
        if level=='FM'and abs(c['center_m'][0])>.008:continue
        row=dict(contourId=i,**c,sourceSectionFaceProvenance=[sourceFaces[k]for k in c['sourceFaceIds']],nonseparatingHandleCycle=None)
        if topology['boundaryEdges']==0 and topology['nonmanifoldEdges']==0:
            try:
                cycle,proof=handle_cycle(v,f,c['center_m'],axis=(0,1,0),half_length=.025,radius=.018 if level!='FM'else .04,
                    excluded_centers=[q['center_m']for j,q in holes if j!=i],maximum_starts=12,wall_seconds=4.)
                row['nonseparatingHandleCycle']=dict(physicalSourceVertexIds=cycle,worldPoints_m=v[cycle],proof=proof)
            except (ValueError,TimeoutError)as ex:row['cycleFailure']=dict(kind=type(ex).__name__,reason=str(ex))
        else:row['cycleNotRun']='Whole source has open or nonmanifold edges; global handle proof is not licensed.'
        result.append(row)
    atomic_json(out/'evidence.json',dict(schema='native-vertebral-source-aperture-discovery/v1',level=level,sourceCatalogueSha256=inv['sourceCatalogueSha256'],sourceInventorySha256=digest_file(BASE/'source-inventory.json'),
        analysisAssemblyTopology=topology,probePlane=dict(normal=[0,1,0],Y_m=y,sourceArteryStationUsedOnlyAsSearchProbe=probe,notARegionBoundary=True),
        allSourceContours=contours,sectionDiagnostics=diagnostic,selectedApertures=result,sourcePositionsAndAllFacesUnchanged=True,noCandidateOrPartition=True,anatomicalAcceptance=False,
        limitation='A topological throat cycle and a literal empty section contour establish a source passage reference, not an anatomical ostium plane or dural entry. Global source defects and local noncycle section components remain explicit.',elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(evidence=str(out/'evidence.json'),topology=topology,contours=len(contours),sectionBadComponents=len(diagnostic['noncycleSectionComponents']),apertures=[dict(contourId=r['contourId'],center_mm=[x*1000 for x in r['center_m']],clearance_mm=r['centroidClearance_m']*1000,cycle=r['nonseparatingHandleCycle']['proof']if r['nonseparatingHandleCycle']else r.get('cycleFailure',r.get('cycleNotRun')))for r in result],elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
