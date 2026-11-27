"""Finite dynamic program over qualified configurations and explicit links."""
import numpy as np
import time


def configuration_path(configurations, data_costs, links, *, transition_weight=10.,
                       maximum_primary_step_m=.012):
    if len(configurations)<2 or len(data_costs)!=len(configurations)or len(links)!=len(configurations)-1:
        raise ValueError('Ordered configurations, costs and explicit adjacent links required')
    arrays=[np.asarray(x,float)for x in configurations]
    costs=[np.asarray(x,float)for x in data_costs]
    if any(x.ndim!=3 or x.shape[2]!=2 or not np.isfinite(x).all()or not len(x)for x in arrays):
        raise ValueError('Finite nonempty planar source-profile configurations required')
    if any(c.shape!=(len(x),)or not np.isfinite(c).all()for c,x in zip(costs,arrays)):
        raise ValueError('One finite data cost per configuration required')
    if not np.isfinite([transition_weight,maximum_primary_step_m]).all()or min(transition_weight,maximum_primary_step_m)<=0:
        raise ValueError('Positive finite transition weight and primary step bound required')
    running=costs[0].copy();parents=[]
    for index,(first,second,relations)in enumerate(zip(arrays,arrays[1:],links)):
        distance=np.zeros((len(first),len(second)))
        if not relations:raise ValueError('Source ownership must define each transition')
        for a,b,weight in relations:
            if not(0<=a<first.shape[1]and 0<=b<second.shape[1])or weight<=0:raise ValueError('Invalid source-profile correspondence')
            distance+=weight*np.sum((first[:,None,a]-second[None,:,b])**2,axis=2)
        primary=np.linalg.norm(first[:,None,0]-second[None,:,0],axis=2)
        proposal=running[:,None]+transition_weight*distance+costs[index+1][None]
        proposal[primary>maximum_primary_step_m]=np.inf
        parent=np.argmin(proposal,axis=0);running=proposal[parent,np.arange(len(second))]
        if not np.isfinite(running).any():return None,dict(feasible=False,failedTransition=index)
        parents.append(parent)
    pick=int(np.argmin(running));path=[pick]
    for parent in reversed(parents):pick=int(parent[pick]);path.append(pick)
    path.reverse();return path,dict(feasible=True,cost=float(running[path[-1]]),states=[len(x)for x in arrays],maximumPrimaryStep_m=max(float(np.linalg.norm(arrays[i][path[i],0]-arrays[i+1][path[i+1],0]))for i in range(len(path)-1)))


def curved_primary_configuration_path(configurations, data_costs, links, station_y_m, *,
                                      maximum_curvature_per_m=80.,
                                      maximum_primary_step_m=.00475,
                                      transition_weight=20.,time_limit_seconds=10.,
                                      measure_minimum_maximum_curvature=False):
    """Second-order finite search, with literal 3D primary-centre curvature.

    Source links still price every constituent. The primary geometric bound
    applies to station triplets, not to clinical curvature or all mesh faces.
    """
    start=time.monotonic();arrays=[np.asarray(x,float)for x in configurations];costs=[np.asarray(x,float)for x in data_costs];ys=np.asarray(station_y_m,float)
    if len(arrays)<3 or len(costs)!=len(arrays)or len(links)!=len(arrays)-1 or ys.shape!=(len(arrays),)or np.any(np.diff(ys)<=0)or not np.isfinite(ys).all():
        raise ValueError('Increasing finite source stations and complete state/link arrays required')
    if any(x.ndim!=3 or x.shape[2]!=2 or not np.isfinite(x).all()or not len(x)or len(x)>256 for x in arrays)or any(c.shape!=(len(x),)or not np.isfinite(c).all()for c,x in zip(costs,arrays)):
        raise ValueError('Finite bounded configuration states and costs required')
    if not np.isfinite([maximum_curvature_per_m,maximum_primary_step_m,transition_weight,time_limit_seconds]).all()or min(maximum_curvature_per_m,maximum_primary_step_m,transition_weight,time_limit_seconds)<=0 or time_limit_seconds>20:
        raise ValueError('Positive finite curvature, step, objective and native deadline required')
    transitions=[]
    for first,second,relation in zip(arrays,arrays[1:],links):
        if not relation:raise ValueError('Missing explicit source transition')
        d=np.zeros((len(first),len(second)))
        for a,b,w in relation:
            if not(0<=a<first.shape[1]and 0<=b<second.shape[1])or w<=0:raise ValueError('Invalid source link')
            d+=w*np.sum((first[:,None,a]-second[None,:,b])**2,axis=2)
        d*=transition_weight;d[np.linalg.norm(first[:,None,0]-second[None,:,0],axis=2)>maximum_primary_step_m]=np.inf;transitions.append(d)
    running=costs[0][:,None]+costs[1][None,:]+transitions[0];parents=[]
    if measure_minimum_maximum_curvature:running=np.where(np.isfinite(running),0.,np.inf)
    for index in range(1,len(arrays)-1):
        first,mid,last=arrays[index-1:index+2]
        e1=np.zeros((len(first),len(mid),1,3));e2=np.zeros((1,len(mid),len(last),3))
        e1[:,:,:,0]=mid[None,:,None,0,0]-first[:,None,None,0,0];e1[:,:,:,2]=mid[None,:,None,0,1]-first[:,None,None,0,1];e1[:,:,:,1]=ys[index]-ys[index-1]
        e2[:,:,:,0]=last[None,None,:,0,0]-mid[None,:,None,0,0];e2[:,:,:,2]=last[None,None,:,0,1]-mid[None,:,None,0,1];e2[:,:,:,1]=ys[index+1]-ys[index]
        denominator=np.linalg.norm(e1,axis=3)*np.linalg.norm(e2,axis=3)*np.linalg.norm(e1+e2,axis=3)
        curvature=2*np.linalg.norm(np.cross(e1,e2),axis=3)/denominator
        if measure_minimum_maximum_curvature:
            proposal=np.maximum(running[:,:,None],curvature)
            proposal=np.where(np.isfinite(transitions[index])[None,:,:],proposal,np.inf)
        else:
            proposal=running[:,:,None]+transitions[index][None,:,:]+costs[index+1][None,None,:]
            proposal[curvature>maximum_curvature_per_m+1e-10]=np.inf
        parent=np.argmin(proposal,axis=0);running=np.take_along_axis(proposal,parent[None],axis=0)[0]
        if not np.isfinite(running).any():return None,dict(feasible=False,failedMiddleStation=index,curvatureBound_per_m=maximum_curvature_per_m,allStatesSearchedAtEachTriplet=True,elapsedSeconds=time.monotonic()-start)
        parents.append(parent)
        if time.monotonic()-start>time_limit_seconds:raise TimeoutError('Native second-order configuration search deadline')
    a,b=np.unravel_index(np.argmin(running),running.shape);path=[int(a),int(b)];value=float(running[a,b])
    for parent in reversed(parents):a,b=int(parent[a,b]),a;path.insert(0,a)
    actual=np.array([[arrays[i][pick,0,0],ys[i],arrays[i][pick,0,1]]for i,pick in enumerate(path)]);u=np.diff(actual,axis=0);curvature=2*np.linalg.norm(np.cross(u[:-1],u[1:]),axis=1)/(np.linalg.norm(u[:-1],axis=1)*np.linalg.norm(u[1:],axis=1)*np.linalg.norm(u[:-1]+u[1:],axis=1))
    return path,dict(feasible=True,cost=value,maximumPrimaryCurvature_per_m=float(curvature.max()),maximumPrimaryStep_m=float(np.linalg.norm(np.diff(actual[:,[0,2]],axis=0),axis=1).max()),curvatureBound_per_m=maximum_curvature_per_m,allStatesSearchedAtEachTriplet=True,wholeSurfaceSmoothnessNotEstablished=True,diagnosticOnly=bool(measure_minimum_maximum_curvature),elapsedSeconds=time.monotonic()-start)
