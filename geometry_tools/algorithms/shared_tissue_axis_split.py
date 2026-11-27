"""Split explicitly named finite tissue cells using shared plane intersections.

Each existing interface is clipped once and reused by both owners. New tissue
cut faces support nested loops (tube lumens); they never fill a hole with an
organ-centre fan. The caller supplies anatomical ownership and must still check
stored topology, crossings, caps, references and actual sections after export.
"""
from collections import Counter
import numpy as np
import manifold3d as md


def _faces(surfaces,refs):
    return np.concatenate([surfaces[q['surface']]if q['side']==1 else surfaces[q['surface']][:,[0,2,1]]for q in refs])


def _boundary(faces):
    uses=Counter(tuple(sorted((int(a),int(b))))for t in faces for a,b in zip(t,np.roll(t,-1)))
    if any(v>2 for v in uses.values()):raise ValueError('Nonmanifold tissue boundary')
    return [(int(a),int(b))for t in faces for a,b in zip(t,np.roll(t,-1))if uses[tuple(sorted((int(a),int(b))))]==1]


def _volume(v,f,origin):
    p=v[f]-origin
    terms=np.einsum('ti,ti->t',p[:,0],np.cross(p[:,1],p[:,2]))/6
    return float(terms.sum()),float(abs(terms).sum())


def split_shared_cells_at_axis(positions,surfaces,cells,split_names,*,axis=None,coordinate=0.,normal=None):
    source=np.asarray(positions,float);original={k:np.asarray(f,int)for k,f in surfaces.items()}
    if not np.isfinite(coordinate)or(normal is None and axis not in [0,1,2]):raise ValueError('Explicit finite plane required')
    if normal is None:n=np.eye(3)[axis]
    else:
        n=np.asarray(normal,float);length=np.linalg.norm(n)
        if n.shape!=(3,)or not np.isfinite(n).all()or length==0:raise ValueError('Finite nonzero plane normal required')
        n=n/length;coordinate=coordinate/length
    distance=source@n-coordinate
    if not set(split_names)<=set(cells):raise ValueError('Unknown tissue to split')
    new_names=[q for pair in split_names.values()for q in pair]
    if len(set(new_names))!=len(new_names)or set(new_names)&set(cells):raise ValueError('Distinct new tissue names required')
    points=source.tolist();cache={};output={};ancestry={};clipped={};area_checks=[]
    def cross(a,b):
        key=tuple(sorted((int(a),int(b))))
        if key not in cache:
            aa,bb=source[list(key)];da,db=distance[list(key)];fraction=da/(da-db)
            if not 0<fraction<1:raise ValueError('Cut intersection is not inside its source edge')
            p=aa+(bb-aa)*fraction
            if normal is None:p[axis]=coordinate
            else:p+=n*(coordinate-p@n)
            cache[key]=len(points);points.append(p.tolist())
        return cache[key]
    for sid,faces in original.items():
        byhalf={-1:[],1:[]};parents={-1:[],1:[]}
        for ordinal,face in enumerate(faces):
            for half in [-1,1]:
                polygon=[]
                for a,b in zip(face,np.roll(face,-1)):
                    da,db=distance[a],distance[b]
                    if half*da>=0:polygon.append(int(a))
                    if da*db<0:polygon.append(cross(a,b))
                clean=[]
                for k in polygon:
                    if not clean or clean[-1]!=k:clean.append(k)
                if len(clean)>1 and clean[0]==clean[-1]:clean.pop()
                for j in range(1,len(clean)-1):
                    tri=[clean[0],clean[j],clean[j+1]]
                    if len(set(tri))==3:byhalf[half].append(tri);parents[half].append(ordinal)
        clipped[sid]={}
        for half,word in [(-1,'negative'),(1,'positive')]:
            fs=np.array(byhalf[half],int).reshape(-1,3)
            if len(fs):
                key=sid+':'+word;output[key]=fs;ancestry[key]=np.array(parents[half],int);clipped[sid][half]=key
    v=np.array(points,float)
    for sid,faces in original.items():
        tri=source[faces];area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2;replayed=np.zeros(len(faces));perimeter=np.linalg.norm(tri-np.roll(tri,-1,axis=1),axis=2).sum(1);new_perimeter=np.zeros(len(faces))
        for key in clipped[sid].values():
            q=v[output[key]];np.add.at(replayed,ancestry[key],np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)/2);np.add.at(new_perimeter,ancestry[key],np.linalg.norm(q-np.roll(q,-1,axis=1),axis=2).sum(1))
        coordinate_scale=max(float(abs(source[faces]).max()),abs(float(coordinate)),float(np.linalg.norm(tri-np.roll(tri,-1,axis=1),axis=2).max()));guard=128*np.finfo(float).eps*coordinate_scale*(perimeter+new_perimeter)+128*np.finfo(float).eps*area
        if np.any(abs(replayed-area)>guard):raise ValueError('Source tissue area lost or duplicated at '+sid)
        area_checks.append(dict(surface=sid,maximumAreaResidualMetersSquared=float(abs(replayed-area).max()),maximumAreaGuardMetersSquared=float(guard.max())))
    result={};cuts=[]
    for name,refs in cells.items():
        source_faces=_faces(original,refs);y=distance[np.unique(source_faces)];crosses=y.min()<0<y.max()
        if _boundary(source_faces):raise ValueError('Input tissue is open before the plane cut: '+name)
        if name in split_names and not crosses:raise ValueError('Named tissue is not divided by the plane: '+name)
        if name not in split_names:
            # Retain both halves of an adjoining tissue as one cell. It uses
            # the same subdivided interface, so neither a T-junction nor an
            # unnecessary regional tissue partition is introduced.
            uses=[]
            for q in refs:
                for half,key in clipped[q['surface']].items():
                    use=dict(surface=key,side=q['side']);peer=q.get('neighbor')
                    if peer is not None:use['neighbor']=split_names[peer][0 if half==-1 else 1]if peer in split_names else peer
                    uses.append(use)
            result[name]=uses;continue
        for half in [-1,1]:
            target=split_names[name][0 if half==-1 else 1];uses=[]
            for q in refs:
                key=clipped[q['surface']].get(half)
                if key is None:continue
                use=dict(surface=key,side=q['side']);peer=q.get('neighbor')
                if peer is not None:use['neighbor']=split_names[peer][0 if half==-1 else 1]if peer in split_names else peer
                uses.append(use)
            result[target]=uses
        negative,positive=split_names[name];fs=_faces(output,result[negative]);boundary=_boundary(fs);follow={a:b for a,b in boundary}
        if not boundary or len(follow)!=len(boundary)or len(set(follow.values()))!=len(follow):raise ValueError('Cut is not simple directed loops')
        todo=set(follow);loops=[]
        while todo:
            first=min(todo);loop=[];j=first
            while j in todo:todo.remove(j);loop.append(j);j=follow[j]
            if j!=first or len(loop)<3:raise ValueError('Incomplete tissue cut loop')
            loops.append(np.array(loop[::-1],int))
        guard=512*np.finfo(float).eps*max(float(abs(v).max()),abs(coordinate),1.)
        if not all(np.all(abs(v[q]@n-coordinate)<=guard)for q in loops):raise ValueError('An unaccounted boundary survives the cut')
        u=np.cross(n,np.eye(3)[np.argmin(abs(n))]);u/=np.linalg.norm(u);w=np.cross(n,u);basis=np.c_[u,w]
        anchor=v[loops[0][0]];polygons=[(v[q]-anchor)@basis for q in loops];flat=np.concatenate(loops);indices=np.asarray(md.triangulate(polygons,epsilon=0.,allow_convex=False),int);caps=flat[indices];cap_normal=np.cross(v[caps[:,1]]-v[caps[:,0]],v[caps[:,2]]-v[caps[:,0]])
        if not len(caps)or not np.all(cap_normal@n>0):raise ValueError('Invalid annular tissue cap orientation')
        sid='partition:'+name;output[sid]=caps;result[negative].append(dict(surface=sid,side=1,neighbor=positive));result[positive].append(dict(surface=sid,side=-1,neighbor=negative))
        for child in [negative,positive]:
            unclosed=_boundary(_faces(output,result[child]))
            if unclosed:
                error=ValueError('A split tissue remains open: '+child+'; '+str(len(unclosed))+' boundary edges')
                error.diagnostics=dict(cell=child,sourceCell=name,sourceCutBoundary=boundary,capFaces=caps.tolist(),remainingEdges=unclosed,loops=[q.tolist()for q in loops],positions=v.tolist())
                raise error
        origin=source[np.unique(source_faces)].mean(0);before,absolute=_volume(source,source_faces,origin);low,absl=_volume(v,_faces(output,result[negative]),origin);high,absh=_volume(v,_faces(output,result[positive]),origin);guard=512*np.finfo(float).eps*(absolute+absl+absh)
        if min(low,high)<=0 or abs(low+high-before)>guard:raise ValueError('Finite tissue volume is not conserved')
        cuts.append(dict(cell=name,children=[negative,positive],directedCutLoops=len(loops),capTriangles=len(caps),nestedHolesRetained=True,sourceVolumeMetersCubed=before,negativeVolumeMetersCubed=low,positiveVolumeMetersCubed=high,volumeResidualMetersCubed=low+high-before,numericalSummationGuardMetersCubed=guard))
    owners={}
    for name,refs in result.items():
        for use in refs:owners.setdefault(use['surface'],[]).append((name,use))
    for sid,uses in owners.items():
        if len(uses)==2:
            (a,x),(b,y)=uses
            if x.get('neighbor')!=b or y.get('neighbor')!=a or x['side']!=-y['side']:raise ValueError('Nonreciprocal tissue interface')
        elif len(uses)!=1 or uses[0][1].get('neighbor'):raise ValueError('Unmatched tissue ownership')
    output={sid:fs for sid,fs in output.items()if sid in owners}
    # Retain exact source ancestry for transporting authored corner fields.
    # A regional plane must not replace an organ's existing smooth normals.
    edge_intersections=[dict(control=idx,edge=list(edge),fraction=float(distance[edge[0]]/(distance[edge[0]]-distance[edge[1]])))for edge,idx in cache.items()]
    return dict(positions=v,surfaces=output,cells=result,
        source_surface_ancestry={sid:ancestry[sid]for sid in output if sid in ancestry},
        edge_intersections=edge_intersections,
        evidence=dict(axis=axis,normal=n.tolist(),coordinateMeters=float(coordinate),sourceControlsLiteral=bool(np.array_equal(v[:len(source)],source)),newIntersectionControls=len(cache),sourceFaceAreaChecks=area_checks,cuts=cuts,sharedSurfacesStoredOnce=True,noLumenDisk=True,anatomicallyAccepted=False))
