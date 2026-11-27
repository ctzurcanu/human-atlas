"""Bounded, capacity-checked source routes through literal section links."""
import numpy as np
from scipy import sparse
from scipy.optimize import linprog


def layered_source_routes(positions,links,*,route_count=2,time_limit_seconds=5.):
 p=[np.asarray(v,float)for v in positions]
 if len(p)<2 or len(links)!=len(p)-1 or not isinstance(route_count,int)or route_count<1 or not 0<time_limit_seconds<=5 or any(v.ndim!=2 or v.shape[1]!=2 or not np.isfinite(v).all()for v in p):raise ValueError('Finite layered points, explicit links, route count and native deadline required')
 if any(len(v)<route_count for v in p):return None,dict(feasible=False,reason='Insufficient distinct source profiles in a section')
 offsets=np.r_[0,np.cumsum([len(v)for v in p])];source=2*int(offsets[-1]);sink=source+1;nodes=sink+1;edges=[];cost=[]
 for k,v in enumerate(p):
  for i in range(len(v)):n=2*(int(offsets[k])+i);edges.append((n,n+1));cost.append(0.)
 for i in range(len(p[0])):edges.append((source,2*i));cost.append(0.)
 for i in range(len(p[-1])):edges.append((2*(int(offsets[-2])+i)+1,sink));cost.append(0.)
 transition_edges=[]
 for k,relations in enumerate(links):
  current=[]
  for relation in relations:
   a,b=int(relation[0]),int(relation[1])
   if not(0<=a<len(p[k])and 0<=b<len(p[k+1])):raise ValueError('Source route link outside its literal section')
   current.append((len(edges),a,b));edges.append((2*(int(offsets[k])+a)+1,2*(int(offsets[k+1])+b)));cost.append(float(np.sum((p[k][a]-p[k+1][b])**2)/.005**2))
  transition_edges.append(current)
 row=[];col=[];val=[]
 for j,(a,b)in enumerate(edges):row.extend([a,b]);col.extend([j,j]);val.extend([-1.,1.])
 matrix=sparse.csc_matrix((val,(row,col)),shape=(nodes,len(edges)));rhs=np.zeros(nodes);rhs[source]=-route_count;rhs[sink]=route_count
 result=linprog(cost,A_eq=matrix,b_eq=rhs,bounds=(0,1),method='highs',options=dict(time_limit=time_limit_seconds))
 proof=dict(feasible=False,status=int(result.status),message=result.message,timeLimitSeconds=time_limit_seconds,vertices=nodes,sourceEdges=len(edges),costUsesOnlySourceCourse=True,anatomicalIdentityNotInferred=True)
 if not result.success:return None,proof
 flow=np.rint(result.x);residual=float(np.max(abs(matrix@flow-rhs)));integrality=float(np.max(abs(result.x-flow)))
 if integrality>1e-8 or residual>1e-8:raise ValueError('Literal source-route capacity or flow replay failed')
 selected=[]
 for k,v in enumerate(p):
  ids=[i for i in range(len(v))if flow[int(offsets[k])+i]>.5];assert len(ids)==route_count;selected.append(ids)
 routes=[[i]for i in selected[0]]
 for current in transition_edges:
  chosen={a:b for edge,a,b in current if flow[edge]>.5};assert len(chosen)==route_count
  for route in routes:route.append(chosen[route[-1]])
 proof.update(feasible=True,flowResidual=residual,integralityError=integrality,objective=float(np.dot(cost,flow)),distinctProfilesPerSection=route_count,sourceRoutes=routes)
 return selected,proof
