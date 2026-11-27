"""Bounded Float32 rounding repair on internal, reciprocally owned controls.

Only interior tissue interfaces may move; physiological boundary and literal
source controls are fixed. Each accepted change reduces a measured defect set,
preserves incident face orientation, and introduces no new local crossings in
any adjoining cell. Full encoded cell audits remain mandatory afterwards.
"""
import numpy as np
from .model import Mesh
from .contact import pair_contacts
from .distance import triangle_pairs

def repair(positions,surfaces,cells,origin,protected_prefix,targets,*,maximum_m=1e-7,iterations=8):
    q=np.array(positions,float,copy=True);initial=q.copy();origin=np.asarray(origin,float);edits=[]
    fs={n:np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])for n,uses in cells.items()}
    protected=set(range(protected_prefix))
    for name,uses in cells.items():
        for u in uses:
            if not u.get('neighbor')or name.split(':')[0]not in ['collecting-wall','upj','ureter-wall']or u['neighbor'].split(':')[0]not in ['collecting-wall','upj','ureter-wall']:
                protected.update(map(int,np.unique(surfaces[u['surface']])))
    def defects(name):
        arr,report=pair_contacts(Mesh(q,fs[name],'encoded'),Mesh(q,fs[name],'encoded'),0,selfcheck=True,max_pairs=1000000)
        if report['truncated']:raise ValueError('Incomplete interface rounding check')
        return arr['pairs'][np.isin(arr['relation'],['crossing','coplanar-overlap'])]
    def local_count(at,point):
        total=0
        for name,f in fs.items():
            incident=np.flatnonzero(np.any(f==point,axis=1))
            if not len(incident):continue
            arr,r=pair_contacts(Mesh(at,f[incident],'encoded'),Mesh(at,f,'encoded'),0,max_pairs=1000000)
            if r['truncated']:raise ValueError('Incomplete local interface rounding check')
            bad=np.isin(arr['relation'],['crossing','coplanar-overlap'])
            if np.any(bad):
                pairs=arr['pairs'][bad];total+=int(np.sum(incident[pairs[:,0]]!=pairs[:,1]))
        return total
    for name in targets:
        for step in range(iterations):
            bad=defects(name)
            if not len(bad):break
            candidates=sorted(set(map(int,np.unique(fs[name][bad])))-protected);best=None
            for point in candidates:
                incident=[f[np.any(f==point,axis=1)]for f in fs.values()];incident=[f for f in incident if len(f)];near=np.concatenate(incident);before=q[near];normals=np.cross(before[:,1]-before[:,0],before[:,2]-before[:,0]);previous=q[point].copy();local=None
                for axis in range(3):
                    value=np.float32(previous[axis]-origin[axis])
                    for direction in [-1,1]:
                        for ulps in [1,4,16]:
                            shifted=value
                            for _ in range(ulps):shifted=np.nextafter(shifted,np.float32(direction*np.inf))
                            q[point,axis]=float(shifted)+origin[axis]
                            if np.linalg.norm(q[point]-initial[point])>maximum_m:q[point]=previous;continue
                            p=q[near];nn=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);valid=np.all(np.einsum('ij,ij->i',nn,normals)>0)
                            if valid:
                                pairs=fs[name][bad];a=triangle_pairs(q[pairs[:,0]],q[pairs[:,1]]);count=int(np.isin(a['relation'],['crossing','coplanar-overlap']).sum())
                                if count<len(bad)and(best is None or count<best[0]):
                                    if local is None:local=local_count(np.where(np.arange(len(q))[:,None]==point,previous,q),point)
                                    changed=local_count(q,point)
                                    if changed<=local:best=(count,point,q[point].copy(),len(bad),local,changed)
                            q[point]=previous
            if best is None:raise ValueError('No bounded internal rounding correction for '+name)
            count,point,new,old_count,local,changed=best;q[point]=new;edits.append(dict(cell=name,point=point,old=initial[point].tolist(),new=new.tolist(),remaining_known_pairs=count,previous_pairs=old_count,local_before=local,local_after=changed))
        else:raise ValueError('Internal interface rounding repair exceeded its bounded iteration budget')
    return q,dict(edits=edits,maximum_change_m=float(np.linalg.norm(q-initial,axis=1).max()),external_boundary_controls_literal=True,source_controls_literal=True,full_replay_required=True)
