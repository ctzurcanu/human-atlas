"""Literal intersection witnesses for flexible profile transport.

Matching the two barycentric images preserves an existing contact endpoint
without making all incident profiles one rigid body. No tissue union or
anatomical validity follows from these equations.
"""
import numpy as np
from geometry_tools.distance import triangle_pairs,point_triangle,exact_intersection_segment


def barycentric_contact_transport(meshes,contacts):
    endpoints=[];coplanar=[];maximum_error=0.;crossing_count=0
    for contact in contacts:
        first,second=contact['first'],contact['second'];a,b=meshes[first],meshes[second];pairs=np.asarray(contact['pairs'],int);result=triangle_pairs(a.triangles[pairs[:,0]],b.triangles[pairs[:,1]]);flat=[]
        for pair,relation in zip(pairs,result['relation']):
            i,j=map(int,pair)
            if relation=='coplanar-overlap':flat.append([i,j]);continue
            if relation!='crossing':raise ValueError('A declared source contact is not a literal crossing or coplanar overlap')
            segment=exact_intersection_segment(a.triangles[i],b.triangles[j]);assert segment['relation']=='crossing';crossing_count+=1
            for point in segment['endpoints']:
                point=np.asarray(point,float);weights=[]
                for mesh,face in [(a,i),(b,j)]:
                    _,replay,w=point_triangle(point[None],mesh.triangles[face][None]);error=float(np.linalg.norm(replay[0]-point));maximum_error=max(maximum_error,error);assert error<=1e-10;assert np.all(w>=-1e-10)and np.all(w<=1+1e-10);weights.append(w[0].tolist())
                endpoints.append(dict(first=first,second=second,firstFace=i,secondFace=j,firstBarycentric=weights[0],secondBarycentric=weights[1],point_m=point.tolist()))
        if flat:coplanar.append(dict(first=first,second=second,pairs=flat))
    return dict(complete=True,endpoints=endpoints,coplanarContacts=coplanar,crossingPairs=crossing_count,coplanarPairs=sum(len(x['pairs'])for x in coplanar),maximumWitnessReplayError_m=maximum_error,sourceMeshesMerged=False,anatomicalInterfaceEstablished=False)


def endpoint_weight_rows(meshes,vertex_weights,endpoints):
    rows=[]
    for e in endpoints:
        a,b=e['first'],e['second'];wa=np.array(e['firstBarycentric'])@vertex_weights[a][meshes[a].faces[e['firstFace']]];wb=np.array(e['secondBarycentric'])@vertex_weights[b][meshes[b].faces[e['secondFace']]];rows.append(wa-wb)
    return np.asarray(rows).reshape(-1,vertex_weights[0].shape[1])


def endpoint_replay(meshes,endpoints):
    maximum=0.
    for e in endpoints:
        a,b=meshes[e['first']],meshes[e['second']];p=np.array(e['firstBarycentric'])@a.triangles[e['firstFace']];q=np.array(e['secondBarycentric'])@b.triangles[e['secondFace']];maximum=max(maximum,float(np.linalg.norm(p-q)))
    return maximum


def contact_translation_kernel(weights,*,maximum_motion_m,maximum_endpoint_error_m):
    a=np.asarray(weights,float);_,singular,v=np.linalg.svd(a,full_matrices=True);rank=int(np.count_nonzero(singular>max(a.shape,default=1)*np.finfo(float).eps*max(float(singular.max(initial=0.)),1.)));basis=v[rank:].T
    if not basis.shape[1]:raise ValueError('Contact equations leave no movable profile direction')
    residual=float(np.linalg.norm(a@basis,ord=2))if len(a)else 0.;bound=residual*np.sqrt(a.shape[1])*maximum_motion_m
    if bound>maximum_endpoint_error_m:raise ValueError('Contact kernel does not certify the original endpoint tolerance')
    return basis,dict(complete=True,originalStations=a.shape[1],parameters=basis.shape[1],rank=rank,maximumEndpointErrorBound_m=bound,originalTolerance_m=maximum_endpoint_error_m,coordinatesWelded=False)
