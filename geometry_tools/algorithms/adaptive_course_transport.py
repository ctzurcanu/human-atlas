"""Conforming red/green refinement for a nonlinear ambient course map.

Each edge/face is measured against the actual map, with explicit error, face
and iteration budgets. Parent-face ownership is retained. This is a geometric
interpolation diagnostic, not anatomical acceptance or a lumen proof.
"""
import numpy as np

def refine_course(v,faces,effective_normals,field,*,deviation_m=1e-6,max_faces=300000,max_rounds=9,fine_parents=(),fine_deviation_m=1e-7):
    v=np.array(v,float,copy=True);f=np.array(faces,int,copy=True);n=np.array(effective_normals,float,copy=True)
    parent=np.arange(len(f));history=[]
    for iteration in range(max_rounds):
        mapped=field(v);raw=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
        # Source banks can split normals at coincident geometric vertices.
        # Refinement must propagate across geometric edges, not attribute IDs.
        _,first,canonical=np.unique(v,axis=0,return_index=True,return_inverse=True)
        geometric,inverse=np.unique(np.sort(canonical[raw],axis=1),axis=0,return_inverse=True);edges=first[geometric]
        mids=v[edges].mean(1);mid_error=np.linalg.norm(field(mids)-mapped[edges].mean(1),axis=1)
        face_error=np.linalg.norm(field(v[f].mean(1))-mapped[f].mean(1),axis=1)
        fe=inverse.reshape(3,len(f)).T;face_tolerance=np.where(np.isin(parent,fine_parents),fine_deviation_m,deviation_m)
        edge_tolerance=np.full(len(edges),deviation_m);np.minimum.at(edge_tolerance,fe.reshape(-1),np.repeat(face_tolerance,3))
        selected=mid_error>edge_tolerance
        selected[fe[face_error>face_tolerance].reshape(-1)]=True
        # Longest-edge closure prevents a green split from repeatedly moving
        # the same unresolved long chord into its neighbouring thin triangle.
        lengths=np.linalg.norm(v[edges[:,1]]-v[edges[:,0]],axis=1);longest=fe[np.arange(len(f)),lengths[fe].argmax(1)]
        for closure in range(100):
            needed=longest[np.any(selected[fe],axis=1)];before=int(selected.sum());selected[needed]=True
            if int(selected.sum())==before:break
        else:raise ValueError('Longest-edge closure budget reached')
        history.append(dict(faces=len(f),vertices=len(v),maximumEdgeDeviationMetres=float(mid_error.max()),maximumFaceDeviationMetres=float(face_error.max()),splitEdges=int(selected.sum())))
        if not np.any(selected):return v,f,n,parent,dict(passedInterpolationBound=True,deviationMetres=deviation_m,fineSourceFaces=list(map(int,fine_parents)),fineDeviationMetres=fine_deviation_m,rounds=history)
        if len(f)+int(selected[fe].sum())>max_faces:raise ValueError('Adaptive course face budget reached')
        new_id=np.full(len(edges),-1,int);new_id[selected]=np.arange(len(v),len(v)+selected.sum())
        a,b=edges[selected].T;newnormal=n[a]+n[b];length=np.linalg.norm(newnormal,axis=1);newnormal[length==0]=n[a[length==0]];newnormal/=np.linalg.norm(newnormal,axis=1)[:,None]
        v=np.r_[v,mids[selected]];n=np.r_[n,newnormal];out=[];owners=[]
        for i,(a,b,c)in enumerate(f):
            ab,bc,ca=new_id[fe[i]];mask=(ab>=0)*1+(bc>=0)*2+(ca>=0)*4
            rows={0:[[a,b,c]],1:[[a,ab,c],[ab,b,c]],2:[[a,b,bc],[a,bc,c]],4:[[a,b,ca],[ca,b,c]],
                3:[[a,ab,c],[ab,bc,c],[ab,b,bc]],5:[[a,ab,ca],[ab,b,c],[ab,c,ca]],6:[[a,b,ca],[b,bc,ca],[bc,c,ca]],
                7:[[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]}[mask]
            out.extend(rows);owners.extend([parent[i]]*len(rows))
        next_f=np.array(out,int)
        f=next_f;parent=np.array(owners,int)
    i=int(mid_error.argmax())
    raise ValueError('Adaptive course iteration budget reached: '+str(dict(history=history,worstEdge=v[edges[i]].tolist(),edgeLength=float(np.linalg.norm(v[edges[i,0]]-v[edges[i,1]])))))
