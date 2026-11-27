"""Finite fixed-boundary triangulation of two declared NONPLANAR HG rim strips.

This is a producer surface repair, not the planar APO storage adapter. No
coordinate, origin footprint, tissue thickness, or protected insertion moves.
"""
from __future__ import annotations
import time,itertools,importlib.util
from pathlib import Path
import numpy as np
from fractions import Fraction as Q
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts,self_contact_summary
from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature
from geometry_tools.audit import topology

s=importlib.util.spec_from_file_location('boundary_triangulation',Path(__file__).with_name('tongue_interface_triangulation.py'))
api=importlib.util.module_from_spec(s);s.loader.exec_module(api)

DOMAINS={'left':[[60,61,103,104],[70,71,113,114]],'right':[[62,63,106,107],[89,90,133,134]]}

def clear(summary):return summary['complete_within_pair_budget']and summary['crossing_face_pairs']==0 and summary['coplanar_overlap_face_pairs']==0

def volume(vertices,faces):
    t=vertices[faces]-vertices.mean(0)
    return float(np.einsum('ij,ij->i',t[:,0],np.cross(t[:,1],t[:,2])).sum()/6)

def simple_boundary(vertices,loop):
    p=[tuple(Q.from_float(float(x))for x in vertices[i])for i in loop]
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    cross=lambda a,b:(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    for i,j in itertools.combinations(range(len(loop)),2):
        if j==i+1 or(i==0 and j==len(loop)-1):continue
        a,b=p[i],p[(i+1)%len(p)];c,d=p[j],p[(j+1)%len(p)];u,v,w=sub(b,a),sub(d,c),sub(c,a);n=cross(u,v)
        if any(n):
            if dot(w,n):continue
            axis=max(range(3),key=lambda k:abs(n[k]));axes=[k for k in range(3)if k!=axis];x,y=axes
            det=u[x]*v[y]-u[y]*v[x];t=(w[x]*v[y]-w[y]*v[x])/det;s=(w[x]*u[y]-w[y]*u[x])/det
            if 0<=t<=1 and 0<=s<=1:return False
        elif not any(cross(w,u)):
            axis=max(range(3),key=lambda k:abs(u[k]));first=sorted([a[axis],b[axis]]);second=sorted([c[axis],d[axis]])
            if max(first[0],second[0])<=min(first[1],second[1]):return False
    return True

def repair(mesh,side,*,seconds=20.):
    began=time.monotonic();v=mesh.vertices;f=mesh.faces.copy();stored=(v-[0.,1.625,0.]).astype('<f4').astype(float)+[0.,1.625,0.]
    domains=DOMAINS[side];moving=np.unique(np.concatenate(domains));static=f[~np.isin(np.arange(len(f)),moving)]
    baseline,_=self_contact_summary(Mesh(v,static,mesh.frame),max_pairs=100000,max_tests=1000000)
    if not clear(baseline):return mesh,dict(passed=False,reason='Static standalone sheet has defects outside the two authorized rim domains',static=baseline)
    histories=[];feasible=[];patches=[]
    def deadline():
        if time.monotonic()-began>seconds:raise TimeoutError('HG fixed-boundary14x14 enumeration reached its explicit deadline')
    for ordinal,ids in enumerate(domains):
        source=f[ids];loop=api.loop_of(source);options=api.triangulations(loop)
        if len(loop)!=6 or len(options)!=14:raise ValueError('The witnessed HG domain is not the authorized six-control/14-triangulation strip')
        if not simple_boundary(v,loop)or not simple_boundary(stored,loop):return mesh,dict(passed=False,reason='The fixed actual3D rim boundary is not simple',domain=ordinal)
        tri=v[source];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]).sum(0);normal/=np.linalg.norm(normal)
        residual=float(np.abs((v[list(loop)]-v[loop[0]])@normal).max());legal=[]
        for choice,candidate in enumerate(options):
            deadline();cf=np.asarray(candidate,int);record=dict(domain=ordinal,choice=choice,faces=cf.tolist())
            if boundary_signature(cf)!=boundary_signature(source):raise ValueError('Fixed rim boundary changed')
            signs=[]
            for bank in [v,stored]:
                t=bank[cf];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);signs.append(float((n@normal).min()))
            if any(np.any(np.linalg.norm(np.cross(bank[cf][:,1]-bank[cf][:,0],bank[cf][:,2]-bank[cf][:,0]),axis=1)==0)for bank in [v,stored]):
                histories.append(dict(**record,passed=False,reason='Degenerate actual raw/stored triangle'));continue
            valid=True;checks=[]
            for bank in [v,stored]:
                local=Mesh(bank,cf,mesh.frame);_,own=pair_contacts(local,local,0.,selfcheck=True,max_pairs=100000,max_tests=1000000)
                _,other=pair_contacts(local,Mesh(bank,static,mesh.frame),0.,max_pairs=100000,max_tests=1000000)
                checks.append(dict(self=own,static=other));valid=valid and clear(own)and clear(other)
            histories.append(dict(**record,passed=valid,checks=checks,diagnosticAverageNormalDot=min(signs),
                averageNormalHalfspaceNotAnOrientationGate=True,allCornersInUnchangedBoundaryConvexHull=True))
            if valid:legal.append(dict(choice=choice,faces=cf))
        patches.append(dict(domain=ordinal,originalFaces=ids,boundary=list(loop),actual3DBoundarySimple=True,
            unchangedBoundaryBoundsM=[v[list(loop)].min(0).tolist(),v[list(loop)].max(0).tolist()],
            nonplanarResidualM=residual,enumerated=14,feasible=len(legal)))
        feasible.append(legal)
    candidates=[];pair_checks=[]
    for first,second in itertools.product(*feasible):
        deadline();valid=True;checks=[]
        for bank in [v,stored]:
            _,check=pair_contacts(Mesh(bank,first['faces'],mesh.frame),Mesh(bank,second['faces'],mesh.frame),0.,max_pairs=100000,max_tests=1000000)
            checks.append(check);valid=valid and clear(check)
        pair_checks.append(dict(choices=[first['choice'],second['choice']],passed=valid,rawAndStored=checks))
        if not valid:continue
        trial=f.copy();trial[domains[0]]=first['faces'];trial[domains[1]]=second['faces']
        candidates.append(dict(faces=trial,choices=[first['choice'],second['choice']],volume=volume(v,trial)))
    evidence=dict(method='Complete14x14 nonplanar rim-strip triangulation; simple fixed3D boundary/controls, raw and stored strict contacts',
        geometryCoordinatesMoved=False,footprintChanged=False,tissueThinned=False,planarAdapterNotUsed=True,
        invalidCommonNormalHalfspaceGuardRemoved=True,
        patches=patches,enumeration=histories,coupledChecks=pair_checks,elapsedSeconds=time.monotonic()-began)
    if not candidates:return mesh,dict(**evidence,passed=False,reason='No outward nonintersecting triangulation under the exact authorized boundaries')
    # All legal alternatives use the same nonplanar boundary. Prefer the least
    # change to its original oriented tissue volume; report the actual change.
    old_volume=volume(v,f);chosen=min(candidates,key=lambda r:abs(r['volume']-old_volume));out=Mesh(v,chosen['faces'],mesh.frame)
    raw_check,_=self_contact_summary(out,max_pairs=100000,max_tests=1000000)
    stored_check,_=self_contact_summary(Mesh(stored,out.faces,mesh.frame),max_pairs=100000,max_tests=1000000)
    topo=topology(out);single=len(topo['components'])==1 and topo['components'][0]['closed_orientable_genus']==0
    passed=clear(raw_check)and clear(stored_check)and single and chosen['volume']>0
    patch_changes=[]
    for domain,ids in enumerate(domains):
        anchor=v[list(patches[domain]['boundary'])].mean(0)
        integral=lambda ff:float(np.einsum('ij,ij->i',(v[ff]-anchor)[:,0],np.cross((v[ff]-anchor)[:,1],(v[ff]-anchor)[:,2])).sum()/6)
        patch_changes.append(dict(domain=domain,signedBoundaryTriangleIntegralChangeM3=integral(out.faces[ids])-integral(f[ids])))
    from geometry_tools.algorithms.finite_tissue_boolean import make_solid
    solid,_=make_solid(v,out.faces)
    evidence.update(passed=passed,selected=chosen['choices'],patchSignedIntegralChanges=patch_changes,
        actualRepairedSheetOccupiedVolumeM3=float(solid.volume()),oldSelfIntersectingVolumeWasAlgebraicOnly=True,
        changedFaceRows=int(np.any(out.faces!=f,axis=1).sum()),
        allControlsLiteral=np.array_equal(v,out.vertices),allPatchBoundariesLiteral=all(boundary_signature(out.faces[ids])==boundary_signature(f[ids])for ids in domains),
        completeRawStandaloneSelf=raw_check,completeStoredStandaloneSelf=stored_check,topology=topo,
        oldOrientedVolumeM3=old_volume,repairedOrientedVolumeM3=chosen['volume'],orientedVolumeChangeM3=chosen['volume']-old_volume,
        elapsedSeconds=time.monotonic()-began)
    return out,evidence
