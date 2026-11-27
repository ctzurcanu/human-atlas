"""Coupled, bounded tissue-clearance proposals for rigid source section cells.

One translation variable per disjoint source cell preserves section calibre
and literal shared controls. Edge regularization couples adjoining cells; it
never removes a source branch. Actual triangle/self/interface checks remain
mandatory: oriented samples alone do not certify a valid vascular fit.
Optional barycentric edge samples expose concave bridges. Explicit face-member
guards permit existing facet contact migration while forbidding new vessel-pair
clashes, intrinsic folds and increased total crossing counts.
"""
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsmr
from .feature_normals import FeatureNormals,metric_oriented_rows


def fit_section_clearance(vertices,faces,sections,surfaces,*,protected=(),
                          gap_m=.00005,step_m=.002,maximum_steps=18,wall_seconds=18.,reference_vertices=None,face_system=None,query_movable_only=False,
                          prepare_step=None,prepared_step=None,trial_scales=(1.,.5,.25,.125,.0625),contact_feedback=None,direction_planes=None,resolve_mixed_contacts=False,coherent_mixed_sections=False,oriented_rim_proposals=(),surface_barycentrics=None,face_member=None,resolve_existing_member_contacts=True):
    started=time.monotonic();original=np.asarray(vertices,float);q=original.copy();f=np.asarray(faces,int)
    if not 0<wall_seconds<=24 or not 0<gap_m<=.001 or not 0<step_m<=.003 or not 1<=maximum_steps<=30:
        raise ValueError('Explicit finite tissue-clearance limits required')
    if not trial_scales or len(trial_scales)>5 or any(not 0<x<=1 for x in trial_scales):
        raise ValueError('Between one and five explicit trial fractions in (0,1] required')
    if (prepare_step is not None or prepared_step is not None or contact_feedback is not None or direction_planes is not None)and maximum_steps!=1:
        raise ValueError('A separately prepared direction is exactly one bounded step')
    owner=np.full(len(q),-1,int)
    for k,ids in enumerate(sections):
        ids=np.asarray(ids,int)
        if np.any(owner[ids]!=-1):raise ValueError('Rigid source cells must be disjoint')
        owner[ids]=k
    loose=np.flatnonzero(owner==-1);owner[loose]=np.arange(len(sections),len(sections)+len(loose))
    n=int(owner.max())+1;fixed=np.zeros(n,bool);fixed[owner[np.asarray(protected,int)]]=True
    section_axis=np.zeros((n,3))
    for ids in sections:
        p=original[ids];p=p-p.mean(0);axis=np.cross(p,np.roll(p,-1,axis=0)).sum(0);length_axis=np.linalg.norm(axis)
        if length_axis>1e-20:section_axis[owner[ids[0]]]=axis/length_axis
    live=np.flatnonzero(~fixed);compact=np.full(n,-1,int);compact[live]=np.arange(len(live))
    cached=prepared_step is not None and 'source_contact_ids'in prepared_step
    if cached:
        for key,value in [('vertices',q),('faces',f),('owner',owner),('fixed',fixed),('reference_vertices',reference_vertices)]:
            if not np.array_equal(prepared_step[key],value):raise ValueError('Prepared contact cache belongs to another '+key)
    feedback_proof=None
    feedback_owner=np.empty((0,2),int);feedback_normal=np.empty((0,3));feedback_rhs=np.empty(0)
    if contact_feedback is not None:
        from .future_contact_planes import future_contact_rows
        chunks=[future_contact_rows(q,f,owner,item,fixed,gap_m=gap_m,complete_pair_planes=bool(item.get('complete_pair_planes',False)))for item in (contact_feedback if isinstance(contact_feedback,list)else[contact_feedback])]
        feedback_owner=np.concatenate([c[0]for c in chunks]);feedback_normal=np.concatenate([c[1]for c in chunks]);feedback_rhs=np.concatenate([c[2]for c in chunks]);feedback_proof=[c[3]for c in chunks]
    if direction_planes is not None:
        if prepared_step is not None:raise ValueError('Explicit direction planes require a fresh source-bound solve')
        po=np.asarray(direction_planes['owner'],int);pn=np.asarray(direction_planes['normal'],float);pr=np.asarray(direction_planes['rhs'],float)
        if po.shape!=(len(pr),2)or pn.shape!=(len(pr),3)or np.any(po<0)or np.any(po>=n)or not np.isfinite(np.r_[pn.ravel(),pr]).all()or np.any(pr>1e-9)or not np.allclose(np.linalg.norm(pn,axis=1),1.,atol=1e-12,rtol=0):
            raise ValueError('Finite zero-feasible unit-normal cell direction planes required')
        feedback_owner=np.r_[feedback_owner,po];feedback_normal=np.r_[feedback_normal,pn];feedback_rhs=np.r_[feedback_rhs,pr]
    edges=np.unique(np.sort(owner[f][:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0)
    edges=edges[edges[:,0]!=edges[:,1]]
    edge_rows=[];edge_cols=[];edge_data=[]
    for axis in range(3):
        for end,sign in [(0,1.),(1,-1.)]:
            valid=~fixed[edges[:,end]]
            edge_rows.extend((np.flatnonzero(valid)*3+axis).tolist())
            edge_cols.extend((compact[edges[valid,end]]*3+axis).tolist())
            edge_data.extend(np.full(int(valid.sum()),sign*.15).tolist())
    smooth=sparse.coo_matrix((edge_data,(edge_rows,edge_cols)),shape=(len(edges)*3,len(live)*3)).tocsr()
    regular=sparse.eye(len(live)*3,format='csr')*.03
    # Coarse facets can bridge a concavity even when vertices and centroid
    # remain inside. Explicit barycentric samples retain source topology and
    # translate the same rigid profile cells; they never claim full coverage.
    bary=np.array([[1/3]*3]if surface_barycentrics is None else surface_barycentrics,float)
    if bary.ndim!=2 or bary.shape[1]!=3 or not 1<=len(bary)<=16 or not np.isfinite(bary).all()or np.any(bary<0)or not np.allclose(bary.sum(1),1,atol=1e-12,rtol=0):
        raise ValueError('One to sixteen finite barycentric surface samples required')
    if cached and not np.array_equal(prepared_step.get('surface_barycentrics',np.array([[1/3]*3])),bary):
        raise ValueError('Prepared direction belongs to another surface sample contract')
    sample_owner=np.concatenate([np.repeat(owner[:,None],3,axis=1),*([owner[f]]*len(bary))])
    sample_weights=np.concatenate([np.tile([1.,0.,0.],(len(q),1)),*[np.tile(w,(len(f),1))for w in bary]])
    moving_sample=np.any(~fixed[sample_owner],axis=1)
    query_ids=np.flatnonzero(moving_sample)if query_movable_only else np.arange(len(sample_owner))
    if any(len(entry)not in (3,4)or len(entry)==4 and not 0<entry[3]<=.001 for entry in surfaces):
        raise ValueError('Each surface needs index, side, name and optional positive gap up to 1 mm')
    feature=[FeatureNormals(entry[0].mesh)for entry in surfaces]
    face_system=np.zeros(len(f),int)if face_system is None else np.asarray(face_system,int)
    if face_system.shape!=(len(f),)or np.any((face_system<0)|(face_system>1)):
        raise ValueError('One arterial/venous class per actual source face required')
    members=None if face_member is None else np.asarray(face_member,int)
    if members is not None and (members.shape!=(len(f),)or np.any(members<0)or reference_vertices is None):
        raise ValueError('Member contact guards require explicit face membership and a reference pose')
    if not isinstance(resolve_existing_member_contacts,bool)or not resolve_existing_member_contacts and(members is None or resolve_mixed_contacts):
        raise ValueError('Retaining existing member contacts requires member guards and an explicit original-contact inventory')
    if cached and members is not None and not np.array_equal(prepared_step.get('face_member'),members):
        raise ValueError('Prepared contact cache belongs to another source-member ownership')
    if cached and bool(prepared_step.get('resolve_existing_member_contacts',True))!=resolve_existing_member_contacts:
        raise ValueError('Prepared direction belongs to another existing-member contact policy')
    source_member_pairs=set()
    source_contact_ids=set()
    source_all_pairs=np.empty((0,2),int)
    if reference_vertices is not None:
        from geometry_tools.contact import self_contact_summary
        from geometry_tools.model import Mesh
        reference=np.asarray(reference_vertices,float)
        if reference.shape!=q.shape:raise ValueError('Corresponding original source coordinates required')
        if cached:
            source_contact_ids=set(prepared_step['source_contact_ids'].tolist())
            if members is not None:
                source_all_pairs=np.asarray(prepared_step['source_all_pairs'],int).reshape(-1,2)
                if np.any(source_all_pairs<0)or np.any(source_all_pairs>=len(f)):
                    raise ValueError('Prepared source contacts contain invalid source faces')
                source_member_pairs={tuple(sorted((int(members[a]),int(members[b]))))for a,b in source_all_pairs}
        else:
            original_contacts,_=self_contact_summary(Mesh(reference,f,'diagnostic:source-cells'),max_tests=1000000,example_limit=25000)
            if not original_contacts['complete_within_pair_budget'] or original_contacts['defect_pairs']>25000:
                raise ValueError('Original source-contact inventory exceeds the finite complete budget')
            source_contact_ids={int(a)*len(f)+int(b)for a,b in original_contacts['examples']if not resolve_mixed_contacts or face_system[a]==face_system[b]}
            source_all_pairs=np.asarray(original_contacts['examples'],int).reshape(-1,2)
            if members is not None:source_member_pairs={tuple(sorted((int(members[a]),int(members[b]))))for a,b in original_contacts['examples']}
    def contact_constraints(v):
        empty_counts=dict(arterial=0,venous=0,mixed=0)
        if members is not None:empty_counts.update(intrinsic=0,newMemberPairs=0)
        if reference_vertices is None:return np.empty((0,2),int),np.empty((0,3)),np.empty(0),0,empty_counts
        if cached and np.array_equal(v,prepared_step['vertices']):
            keys=['arterial','venous','mixed']+(['intrinsic','newMemberPairs']if members is not None else[])
            values=np.asarray(prepared_step['base_counts'])
            if values.shape!=(len(keys),):raise ValueError('Prepared contact counters do not match the declared member guard')
            counts=dict(zip(keys,map(int,values)))
            return prepared_step['base_owner'],prepared_step['base_normal'],prepared_step['base_rhs'],sum(counts[k]for k in ['arterial','venous','mixed']),counts
        audit,_=self_contact_summary(Mesh(v,f,'diagnostic:source-cells'),max_tests=1000000,example_limit=25000)
        if not audit['complete_within_pair_budget']or audit['defect_pairs']>25000:
            raise ValueError('Current cell-contact inventory exceeds the finite complete budget')
        pairs=np.array([pair for pair in audit['examples']if members is not None or pair[0]*len(f)+pair[1]not in source_contact_ids],int).reshape(-1,2)
        if not len(pairs):return np.empty((0,2),int),np.empty((0,3)),np.empty(0),0,empty_counts
        types=face_system[pairs]
        counts=dict(arterial=int(np.sum(np.all(types==0,axis=1))),venous=int(np.sum(np.all(types==1,axis=1))),mixed=int(np.sum(types[:,0]!=types[:,1])))
        total_pairs=len(pairs)
        if members is not None:
            # Facet IDs naturally change when an existing crossing slides onto
            # an adjacent triangle. Count all actual contacts, forbid entirely
            # new vessel-pair clashes, and separately forbid intrinsic folds.
            current_pairs={tuple(sorted((int(members[a]),int(members[b]))))for a,b in pairs}
            counts.update(intrinsic=int(np.sum(members[pairs[:,0]]==members[pairs[:,1]])),newMemberPairs=len(current_pairs-source_member_pairs))
            pairs=pairs[types[:,0]!=types[:,1]];types=face_system[pairs]
            if not resolve_existing_member_contacts:
                # Placement can retain explicitly inventoried source defects
                # while preserving all actual count/new-member/intrinsic guards.
                # Do not turn old microscopic contacts into competing placement
                # targets. They remain counted and anatomically unaccepted.
                keep=np.array([int(a)*len(f)+int(b)not in source_contact_ids for a,b in pairs],bool)
                pairs=pairs[keep];types=face_system[pairs]
            if not len(pairs):return np.empty((0,2),int),np.empty((0,3)),np.empty(0),total_pairs,counts
        a,b=v[f[pairs[:,0]]],v[f[pairs[:,1]]]
        ea=a[:,[1,2,0]]-a;eb=b[:,[1,2,0]]-b
        axes=np.concatenate([np.cross(ea[:,0],ea[:,1])[:,None],np.cross(eb[:,0],eb[:,1])[:,None],
            np.cross(ea[:,:,None],eb[:,None]).reshape(-1,9,3)],axis=1)
        size=np.linalg.norm(axes,axis=2);valid=size>1e-20
        axes=np.divide(axes,size[:,:,None],out=np.zeros_like(axes),where=valid[:,:,None])
        origin=a[:,0];pa=np.einsum('nki,nji->nkj',a-origin[:,None],axes)
        pb=np.einsum('nki,nji->nkj',b-origin[:,None],axes)
        plus=pa.max(1)-pb.min(1)+gap_m;minus=pb.max(1)-pa.min(1)+gap_m
        depths=np.where(valid,np.minimum(plus,minus),np.inf);chosen=depths.argmin(1);row=np.arange(len(pairs))
        normal=axes[row,chosen]*np.where(plus[row,chosen]<=minus[row,chosen],1.,-1.)[:,None]
        if coherent_mixed_sections:
            # Choose one radial escape per physical section pair. Independent
            # facet separating axes can oppose each other around the same tube.
            centres=np.column_stack([np.bincount(owner,weights=v[:,k],minlength=n)/np.bincount(owner,minlength=n)for k in range(3)])
            course=[]
            for face_ids in [pairs[:,0],pairs[:,1]]:
                ids=owner[f[face_ids]];weights=np.ones(ids.shape,float)
                for k in [1,2]:weights[:,k]=np.all(ids[:,:k]!=ids[:,k,None],axis=1)
                course.append(np.einsum('nk,nkj->nj',weights,centres[ids])/weights.sum(1)[:,None])
            mixed=types[:,0]!=types[:,1];arterial_face=np.where(types[:,0]==0,pairs[:,0],pairs[:,1]);tangent=section_axis[owner[f[arterial_face]]]
            direction=course[1]-course[0]
            projection=np.mean(np.einsum('nki,nkj->nkij',tangent,tangent),axis=1)
            direction-=np.einsum('nij,nj->ni',projection,direction);size_direction=np.linalg.norm(direction,axis=1)
            good=mixed&(size_direction>1e-9);normal[good]=direction[good]/size_direction[good,None]
        # All support vertices, rather than triangle centroids, participate.
        # This linearizes actual separating-axis inequalities for each face.
        owners=[];directions=[];deficits=[]
        for i in range(3):
            for j in range(3):
                first=owner[f[pairs[:,0],i]];second=owner[f[pairs[:,1],j]]
                deficit=np.einsum('ni,ni->n',a[:,i]-b[:,j],normal)+gap_m
                take=(deficit>1e-9)&(first!=second)&(~fixed[first]|~fixed[second])
                owners.append(np.c_[first[take],second[take]]);directions.append(normal[take]);deficits.append(deficit[take])
        return np.concatenate(owners),np.concatenate(directions),np.concatenate(deficits),total_pairs,counts
    def constraints(v):
        points=np.concatenate([v,*[np.einsum('j,njk->nk',w,v[f])for w in bary]]);normals=[];deficits=[];known_rows=[];proof=[]
        for entry,fn in zip(surfaces,feature):
            index,side,name=entry[:3];surface_gap=entry[3]if len(entry)==4 else gap_m
            hit,_=index.closest_points(points[query_ids]);normal,rim=fn.at(hit)
            g,distance=metric_oriented_rows(points[query_ids],hit,normal)
            gradient=np.zeros_like(points);gradient[query_ids]=g*side
            local_rim=rim&(hit['distance_m']<=.002)if name in oriented_rim_proposals else np.zeros(len(rim),bool)
            valid=np.zeros(len(points),bool);valid[query_ids]=~rim|local_rim
            deficit=np.zeros(len(points));deficit[query_ids]=surface_gap-distance*side
            proof.append(dict(surface=name,requestedGap_m=float(surface_gap),queriedSamples=len(query_ids),protectedSamplesNotQueried=len(points)-len(query_ids),knownSamples=int((~rim).sum()),unknownRimSamples=int(rim.sum()),orientedRimProposalSamples=int(local_rim.sum()),closedContainmentInferred=False,
                violatingSamples=int(np.sum(valid&(deficit>1e-7))),
                movableViolations=int(np.sum(valid&moving_sample&(deficit>1e-7))),
                protectedViolations=int(np.sum(valid&~moving_sample&(deficit>1e-7))),
                maximumDeficit_m=float(np.maximum(deficit[valid],0).max(initial=0.))))
            normals.append(gradient);deficits.append(deficit);known_rows.append(valid)
        return normals,deficits,known_rows,proof
    history=[];total_shift=np.zeros((n,3))
    for iteration in range(maximum_steps):
        gradients,deficits,known,proof=constraints(q)
        contact_owner,contact_normal,contact_deficit,new_crossings,contact_counts=contact_constraints(q)
        energy=sum(float(np.square(np.maximum(d[k&moving_sample],0)).sum())for d,k in zip(deficits,known))
        history.append(dict(iteration=iteration,energy_m2=energy,surfaces=proof,newCrossingFacePairs=new_crossings,newCrossingClasses=contact_counts,contactFeedback=feedback_proof))
        if not any(x['movableViolations']for x in proof)and not new_crossings:break
        if prepared_step is not None:
            if not np.array_equal(prepared_step['vertices'],q)or not np.array_equal(prepared_step['owner'],owner):
                raise ValueError('Saved direction belongs to another pose or source-cell ownership')
            shift=np.asarray(prepared_step['shift'],float)
            if shift.shape!=(n,3)or not np.isfinite(shift).all()or np.any(shift[fixed]!=0)or np.linalg.norm(shift,axis=1).max(initial=0)>step_m+1e-12:
                raise ValueError('Saved direction violates fixed cells or the declared step bound')
        else:
            rows=[];cols=[];data=[];rhs=[];offset=0
            for gradient,deficit,valid in zip(gradients,deficits,known):
                selected=np.flatnonzero(valid&moving_sample&(deficit>1e-7))
                for corner in range(3):
                    group=sample_owner[selected,corner];movable=~fixed[group]
                    for axis in range(3):
                        rows.extend((np.flatnonzero(movable)+offset).tolist())
                        cols.extend((compact[group[movable]]*3+axis).tolist())
                        data.extend((100*sample_weights[selected[movable],corner]*gradient[selected[movable],axis]).tolist())
                rhs.extend((100*np.minimum(deficit[selected],step_m)).tolist());offset+=len(selected)
            for end,sign in [(0,-1.),(1,1.)]:
                group=contact_owner[:,end];movable=~fixed[group]
                for axis in range(3):
                    rows.extend((np.flatnonzero(movable)+offset).tolist())
                    cols.extend((compact[group[movable]]*3+axis).tolist())
                    data.extend((1000*sign*contact_normal[movable,axis]).tolist())
            rhs.extend((1000*np.minimum(contact_deficit,step_m)).tolist());offset+=len(contact_deficit)
            A=sparse.coo_matrix((data,(rows,cols)),shape=(offset,len(live)*3)).tocsr()
            lhs=sparse.vstack([A,smooth,regular],format='csr')
            desired=np.r_[rhs,-smooth@total_shift[live].ravel(),np.zeros(len(live)*3)]
            solution=lsmr(lhs,desired,atol=1e-7,btol=1e-7,maxiter=80)
            shift=np.zeros((n,3));shift[live]=solution[0].reshape(-1,3)
            length_shift=np.linalg.norm(shift,axis=1)
            shift*=np.minimum(1.,np.divide(step_m,length_shift,out=np.ones(n),where=length_shift>0))[:,None]
            tri=q[f];before=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);length=np.linalg.norm(before,axis=1)
            # A single tight branch cell must not reduce the step on every other
            # finger. Enforce the face guard locally before the global energy test.
            damped=set()
            if (contact_feedback is not None or direction_planes is not None)and len(feedback_rhs):
                from .clearance_direction_projection import project_cell_direction
                shift,projection=project_cell_direction(q,f,owner,fixed,shift,feedback_owner,feedback_normal,feedback_rhs,step_m=step_m,wall_seconds=2.,objective_matrix=lhs,objective_rhs=desired)
                history[-1]['hardContactProjection']=projection
            for guard in range(12):
                t=(q+shift[owner])[f];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
                bad=np.einsum('ij,ij->i',normal,before)<.05*length*length
                if not bad.any():break
                groups=np.unique(owner[f[bad]])
                if contact_feedback is not None:shift*=.5
                else:shift[groups]*=.5
                damped.update(groups.tolist())
            history[-1]['orientationDampedCells']=len(damped)
        tri=q[f];before=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);length=np.linalg.norm(before,axis=1)
        if prepare_step is not None:
            payload=dict(vertices=q.copy(),owner=owner.copy(),shift=shift.copy(),surface_barycentrics=bary.copy())
            if reference_vertices is not None:
                payload.update(faces=f.copy(),fixed=fixed.copy(),reference_vertices=np.asarray(reference_vertices).copy(),
                    source_contact_ids=np.array(sorted(source_contact_ids),np.int64),base_owner=contact_owner.copy(),
                    base_normal=contact_normal.copy(),base_rhs=contact_deficit.copy(),
                    base_counts=np.array([contact_counts[k]for k in ['arterial','venous','mixed']+(['intrinsic','newMemberPairs']if members is not None else[])],np.int64))
                if members is not None:payload.update(face_member=members.copy(),source_all_pairs=source_all_pairs.copy(),resolve_existing_member_contacts=np.array(resolve_existing_member_contacts))
            prepare_step(payload)
            history[-1]['status']='DirectionPreparedNoTrialApplied'
            break
        accepted=False
        history[-1]['trials']=[]
        for scale in trial_scales:
            trial=q+scale*shift[owner];t=trial[f]
            normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
            if np.any(np.einsum('ij,ij->i',normal,before)<.05*length*length):
                history[-1]['trials'].append(dict(scale=scale,status='RejectedOrientation'))
                continue
            _,ds,ks,_=constraints(trial)
            after=sum(float(np.square(np.maximum(d[k&moving_sample],0)).sum())for d,k in zip(ds,ks))
            if after<energy-1e-16 or (new_crossings and after<=energy*1.001+1e-16):
                _,_,_,after_crossings,after_counts=contact_constraints(trial)
                history[-1]['trials'].append(dict(scale=scale,energy_m2=after,newCrossingClasses=after_counts))
                if all(after_counts[key]<=contact_counts[key]for key in contact_counts):
                    q=trial;total_shift+=scale*shift;accepted=True
                    history[-1].update(acceptedScale=scale,afterEnergy_m2=after,afterNewCrossingFacePairs=after_crossings,afterNewCrossingClasses=after_counts)
                    break
            else:history[-1]['trials'].append(dict(scale=scale,energy_m2=after,status='RejectedTissueEnergy'))
            if time.monotonic()-started>wall_seconds:raise TimeoutError('Rigid-section clearance deadline')
        if not accepted:
            history[-1]['status']='NoEnergyDecreasingOrientationPreservingStep';break
        if time.monotonic()-started>wall_seconds:raise TimeoutError('Rigid-section clearance deadline')
    _,_,_,final=constraints(q)
    if time.monotonic()-started>wall_seconds:raise TimeoutError('Rigid-section clearance deadline')
    return q,dict(constructionComplete=True,anatomicalAcceptance=False,fullTriangleChecksRequired=True,
        rigidCells=n,protectedCells=int(fixed.sum()),sourceFacesRemoved=0,
        fixedControlsLiteral=bool(np.array_equal(q[fixed[owner]],original[fixed[owner]])),
        maximumControlTranslation_m=float(np.linalg.norm(total_shift,axis=1).max(initial=0.)),
        sampleGap_m=gap_m,surfaceBarycentrics=bary.tolist(),stepBound_m=step_m,maximumSteps=maximum_steps,
        initial=history[0]['surfaces'],final=final,steps=history,
        originalSourceContactPairsExplicitlyRetained=len(source_contact_ids),resolveOriginalArterialVenousContacts=bool(resolve_mixed_contacts),coherentMixedSectionEscape=bool(coherent_mixed_sections),
        resolveExistingMemberContacts=resolve_existing_member_contacts,
        contactCountPolicy='all-actual-pairs-with-member-and-intrinsic-guards'if members is not None else'new-source-face-pairs',
        originalContactsNotAnatomicallyAccepted=bool(source_contact_ids),elapsedSeconds=time.monotonic()-started)
