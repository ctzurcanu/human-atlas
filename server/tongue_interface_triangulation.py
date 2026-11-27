"""Finite storage-aware triangulation of declared planar tissue fragments.

Every candidate uses the exact existing boundary controls; no point, protected
shared chart, tissue thickness or gap is changed. Complete Catalan enumeration
is small here (4/7/7 boundary controls), so a failed immutable patch is an exact
finite constraint result rather than a heuristic search failure.
"""
from __future__ import annotations
from collections import Counter
from functools import lru_cache
from itertools import combinations
import time
import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature


def loop_of(faces):
    signature=boundary_signature(faces)
    if not signature or any(n!=1 for n in signature.values()):raise ValueError('One simple directed patch boundary is required')
    following={};incoming=Counter()
    for a,b in signature:
        if a in following:raise ValueError('A patch boundary branches')
        following[a]=b;incoming[b]+=1
    if any(incoming[a]!=1 for a in following):raise ValueError('Patch boundary is not a directed cycle')
    start=min(following);order=[];v=start
    while v not in order:
        order.append(v);v=following[v]
    if v!=start or len(order)!=len(following):raise ValueError('Patch needs more than one loop')
    if set(order)!=set(np.asarray(faces).ravel()):raise ValueError('Interior controls require a different explicitly bounded triangulator')
    return tuple(order)


@lru_cache(None)
def triangulations(loop):
    if len(loop)<3:return ((),)
    if len(loop)==3:return ((tuple(loop),),)
    out=[]
    for k in range(1,len(loop)-1):
        for left in triangulations(loop[:k+1]):
            for right in triangulations(loop[k:]):
                out.append(left+right+((loop[0],loop[k],loop[-1]),))
    return tuple(out)


def clear(summary):
    return summary['complete_within_pair_budget']and summary['crossing_face_pairs']==0 and summary['coplanar_overlap_face_pairs']==0


def repair(raw,stored,free_faces,parents,free_words,shared_faces,SL_faces,frame,
           *,parent_faces=(75,78,80),time_limit=45.):
    began=time.monotonic();raw=np.asarray(raw,float);stored=np.asarray(stored,float)
    original=np.asarray(free_faces,int);parents=np.asarray(parents,int)
    protected=boundary_signature(original);target=set(parent_faces)
    static=original[~np.isin(parents,list(target))]
    static_mesh=Mesh(stored,np.r_[static,shared_faces[:,[0,2,1]]],frame)
    SL=Mesh(stored,SL_faces,frame)
    histories=[];feasible={};patches={};all_counts={}
    def deadline():
        if time.monotonic()-began>time_limit:raise TimeoutError('Finite common-chart triangulation reached its explicit deadline')
    for parent in parent_faces:
        deadline();ids=np.flatnonzero(parents==parent);source=original[ids];loop=loop_of(source)
        if len(loop)>8:raise ValueError('This bounded repair supports at most8boundary controls per patch')
        vertex_words={}
        for face,wordrow in zip(source,np.asarray(free_words)[ids]):
            for vertex,word in zip(face,wordrow):
                if int(vertex)in vertex_words and not np.array_equal(vertex_words[int(vertex)],word):
                    raise ValueError('Parent-local normal field disagrees at an existing boundary control')
                vertex_words[int(vertex)]=word
        source_tri=raw[source]
        normal=np.cross(source_tri[:,1]-source_tri[:,0],source_tri[:,2]-source_tri[:,0]).sum(0)
        normal/=np.linalg.norm(normal)
        planarity=float(np.abs((raw[list(loop)]-raw[loop[0]])@normal).max())
        if planarity>2e-12:raise ValueError('A declared parent fragment is not source-planar within2pm')
        candidates=triangulations(loop);all_counts[parent]=len(candidates);feasible[parent]=[]
        for ordinal,candidate in enumerate(candidates):
            deadline();faces=np.asarray(candidate,int).reshape(-1,3)
            if len(faces)!=len(source)or boundary_signature(faces)!=boundary_signature(source):
                raise ValueError('Combinatorial triangulation lost the literal directed boundary')
            rt=raw[faces];st=stored[faces]
            rn=np.cross(rt[:,1]-rt[:,0],rt[:,2]-rt[:,0]);sn=np.cross(st[:,1]-st[:,0],st[:,2]-st[:,0])
            history=dict(parent=int(parent),ordinal=ordinal,faces=faces.tolist(),complete=True)
            if np.any(rn@normal<=0)or np.any(sn@normal<=0):
                histories.append(dict(**history,passed=False,reason='outside-concave-region-or-inverted-at-storage'));continue
            words=np.array([[vertex_words[int(v)]for v in t]for t in faces],dtype='<i2')
            nwords=words.astype(float);nwords/=np.linalg.norm(nwords,axis=2)[:,:,None]
            unit=sn/np.linalg.norm(sn,axis=1)[:,None]
            dots=np.einsum('fki,fi->fk',nwords,unit)
            if np.any(dots<=0):
                histories.append(dict(**history,passed=False,reason='nonoutward-inherited-chart',minimum_normal_dot=float(dots.min())));continue
            mesh=Mesh(stored,faces,frame)
            result,againstSL=pair_contacts(mesh,SL,0.,max_pairs=100000,max_tests=1000000)
            history['SL_contacts']=againstSL
            if not clear(againstSL):
                take=np.flatnonzero((result['relation']=='crossing')|(result['relation']=='coplanar-overlap'))
                history['SL_defect_face_pairs']=result['pairs'][take].tolist()
                histories.append(dict(**history,passed=False,reason='actual-Float32-SL-contact'));continue
            _,against_static=pair_contacts(mesh,static_mesh,0.,max_pairs=100000,max_tests=1000000)
            history['static_aponeurosis_contacts']=against_static
            if not clear(against_static):
                histories.append(dict(**history,passed=False,reason='actual-Float32-static-aponeurosis-contact'));continue
            den=np.sum((st-np.roll(st,-1,axis=1))**2,axis=(1,2))
            quality=float(np.min(2*np.sqrt(3)*np.linalg.norm(sn,axis=1)/den))
            histories.append(dict(**history,passed=True,minimum_quality=quality))
            feasible[parent].append(dict(faces=faces,words=words,ordinal=ordinal,quality=quality,mesh=mesh))
        patches[parent]=dict(parent=int(parent),face_rows=ids.tolist(),boundary_loop=list(loop),
            original_faces=source.tolist(),source_planarity_m=planarity,
            enumerated_all_boundary_triangulations=len(candidates),feasible=len(feasible[parent]))
        feasible[parent].sort(key=lambda c:-c['quality'])
    failed=[p for p in parent_faces if not feasible[p]]
    evidence=dict(method='Complete fixed-boundary combinatorial triangulation, actual Float32 contact/normal screening',
        parameters_changed=False,original_controls_moved=False,shared_interface_chart_changed=False,
        complete_enumeration=True,patches=list(patches.values()),candidate_results=histories,
        failed_parent_faces=list(failed),elapsed_seconds=time.monotonic()-began)
    if failed:
        evidence.update(passed=False,reason='No legal apo-only triangulation under the immutable shared chart',
            minimal_coupled_change_required='Retriangulate the adjacent SL-free source-parent patch jointly with the failed apo patch, or subdivide their reciprocal interface boundary in both owners; no source-coordinate motion or tissue thinning is established as necessary')
        return original.copy(),np.asarray(free_words).copy(),evidence
    compatibility={}
    def compatible(a,b):
        key=tuple(sorted(((a[0],a[1]['ordinal']),(b[0],b[1]['ordinal']))))
        if key not in compatibility:
            _,summary=pair_contacts(a[1]['mesh'],b[1]['mesh'],0.,max_pairs=100000,max_tests=1000000)
            compatibility[key]=dict(passed=clear(summary),summary=summary)
        return compatibility[key]['passed']
    chosen=[]
    def select(index):
        if index==len(parent_faces):return True
        parent=parent_faces[index]
        for candidate in feasible[parent]:
            deadline()
            if all(compatible((parent,candidate),other)for other in chosen):
                chosen.append((parent,candidate))
                if select(index+1):return True
                chosen.pop()
        return False
    possible=select(0)
    evidence['coupled_patch_pair_checks']=[dict(candidates=[list(x)for x in key],**value)for key,value in compatibility.items()]
    if not possible:
        evidence.update(passed=False,reason='Every individually clear patch combination fails an actual coupled patch contact')
        return original.copy(),np.asarray(free_words).copy(),evidence
    output=original.copy();out_words=np.asarray(free_words).copy();selected=[]
    for parent,candidate in chosen:
        ids=np.asarray(patches[parent]['face_rows'],int);output[ids]=candidate['faces'];out_words[ids]=candidate['words']
        selected.append(dict(parent=int(parent),ordinal=candidate['ordinal'],minimum_quality=candidate['quality'],
            original_faces=original[ids].tolist(),repaired_faces=output[ids].tolist()))
    if boundary_signature(output)!=protected:raise ValueError('Common free-chart boundary changed')
    evidence.update(passed=True,selected=selected,free_chart_boundary_literal=True,
        changed_face_rows=int(np.count_nonzero(np.any(output!=original,axis=1))),elapsed_seconds=time.monotonic()-began)
    return output,out_words,evidence
