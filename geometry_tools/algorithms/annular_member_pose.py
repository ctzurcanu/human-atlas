"""Explicit source-strand pose overrides; never infer ownership by proximity."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
from .axial_bundle_pose import quintic


class AnnularMemberPose:
    def __init__(self, specification):
        self.spec = specification
        self.ys = np.asarray(specification['y_m'], float)
        self.support = np.asarray(specification['ySupport_m'], float)
        source = np.asarray(specification['sourceArteryCentreXZ_m'], float)
        target = np.asarray(specification['targetXZRoll'], float)
        if (self.ys.ndim != 1 or len(self.ys) < 2 or np.any(np.diff(self.ys) <= 0)
                or source.shape != (len(self.ys), 2) or target.shape != (len(self.ys), 3)
                or self.support.shape != (4,) or np.any(np.diff(self.support) <= 0)
                or self.support[1] != self.ys[0] or self.support[2] != self.ys[-1]
                or not np.isfinite(np.r_[self.ys,source.ravel(),target.ravel(),self.support]).all()):
            raise ValueError('Ordered source-bound member controls and support required')
        self.source = CubicSpline(self.ys, source, bc_type=((1,np.zeros(2)),(1,np.zeros(2))))
        self.target = CubicSpline(self.ys, target, bc_type=((1,np.zeros(3)),(1,np.zeros(3))))
        self.maximum = float(specification['maximumChange_m'])
        self.displacements=None
        if 'translationsFromReferenceXZ_m'in specification:
            self.displacements=np.asarray(specification['translationsFromReferenceXZ_m'],float)
            if self.displacements.shape!=(len(self.ys),2)or not np.isfinite(self.displacements).all()or np.any(np.linalg.norm(self.displacements,axis=1)>self.maximum):raise ValueError('Source-bound reference displacements exceed the declared member budget')
        self.anchors=np.asarray(specification.get('anchorYIntervals_m',[]),float).reshape(-1,2)
        self.anchor_transition=float(specification.get('anchorTransition_m',.012))
        if not np.isfinite(self.anchors).all()or np.any(self.anchors[:,1]<self.anchors[:,0])or not np.isfinite(self.anchor_transition)or self.anchor_transition<=0:raise ValueError('Explicit ordered source anchor intervals required')
        if not 0 < self.maximum <= .010:
            raise ValueError('Explicit member movement budget must be at most10mm')

    def pose(self, y, base_rotation, base_translation):
        base_rotation, base_translation = np.asarray(base_rotation), np.asarray(base_translation)
        a,b,c,d = self.support
        weight = float(quintic((y-a)/(b-a))*(1-quintic((y-c)/(d-c))))
        if len(self.anchors):
            distance=np.maximum(np.maximum(self.anchors[:,0]-y,y-self.anchors[:,1]),0)
            weight*=float(np.min(quintic(distance/self.anchor_transition)))
        if weight == 0:
            return base_rotation, base_translation
        station = np.clip(y, self.ys[0], self.ys[-1])
        source = self.source(station); target = self.target(station)
        center = np.array([source[0],y,source[1]])
        angle = target[2]; co,si=np.cos(angle),np.sin(angle)
        desired_rotation = np.array([[co,0,si],[0,1,0],[-si,0,co]])
        correction = Rotation.from_matrix(desired_rotation@base_rotation.T).as_rotvec()*weight
        rotation = Rotation.from_rotvec(correction).as_matrix()@base_rotation
        old_center = base_rotation@center+base_translation
        if self.displacements is None:new_center=old_center+weight*(np.array([target[0],y,target[1]])-old_center)
        else:
            interval=int(np.clip(np.searchsorted(self.ys,station,side='right')-1,0,len(self.ys)-2));u=float(quintic((station-self.ys[interval])/(self.ys[interval+1]-self.ys[interval])));delta=(1-u)*self.displacements[interval]+u*self.displacements[interval+1];new_center=old_center+weight*np.array([delta[0],0,delta[1]])
        if np.linalg.norm(new_center-old_center)>self.maximum:
            raise ValueError('Member pose exceeds declared movement budget')
        return rotation, new_center-rotation@center


def member_ownership(specification, part_id, annuli):
    rings, strips = {}, {}
    for member in specification['members']:
        if member['partId'] != part_id:
            continue
        field = AnnularMemberPose(member)
        selected = set(member['ringIds'])
        actual = {i for i,s in enumerate(annuli['strips']) if set(s['rings'])<=selected}
        if actual != set(member['stripIds']):
            raise ValueError('Source strand strip/ring ownership changed')
        if {f for i in actual for f in annuli['strips'][i]['source_face_ids']} != set(member['sourceFaceIds']):
            raise ValueError('Source strand face ownership changed')
        for rid in selected:
            if rid in rings: raise ValueError('Ambiguous source profile ownership')
            rings[rid] = field
        for sid in actual:
            if sid in strips: raise ValueError('Ambiguous source annulus ownership')
            strips[sid] = field
    return rings, strips


def transport_rigid_cells(vertices, faces, cells, targets, locked_cells):
    """Move explicitly targeted rigid source cells, preserving all other locks.

    Free junction points interpolate through literal mesh edges. This does
    not establish orientation, calibre, clearance or anatomical continuity.
    """
    v=np.asarray(vertices,float);f=np.asarray(faces,int);cell=np.asarray(cells,int)
    if cell.shape!=(len(v),) or np.any(cell<0): raise ValueError('Explicit source cell per vertex required')
    count=int(cell.max()+1);locks={int(k):np.zeros(3)for k in locked_cells};offset=np.zeros_like(v)
    for key,(rotation,translation) in targets.items():
        ids=np.flatnonzero(cell==key);R=np.asarray(rotation,float);t=np.asarray(translation,float)
        if not len(ids) or R.shape!=(3,3) or t.shape!=(3,) or not np.isfinite(np.r_[R.ravel(),t]).all() or not np.allclose(R.T@R,np.eye(3),atol=2e-12)or np.linalg.det(R)<1-2e-12:
            raise ValueError('An owned cell and finite proper pose required')
        transformed=v[ids]@R.T+t
        if key in locks and np.max(np.linalg.norm(transformed-v[ids],axis=1))>5e-8:
            raise ValueError('Member pose would move a protected source interface')
        center=v[ids].mean(0);offset[ids]=(v[ids]-center)@(R-np.eye(3)).T;locks[int(key)]=R@center+t-center
    edges=f[:,[0,1,1,2,2,0]].reshape(-1,2);end=cell[edges];length=np.linalg.norm(v[edges[:,0]]-v[edges[:,1]],axis=1);keep=(end[:,0]!=end[:,1])&(length>0);edges,end,length=edges[keep],end[keep],length[keep];weights=1/np.maximum(length,1e-5)**2
    adjacency=coo_matrix((np.r_[weights,weights],(np.r_[end[:,0],end[:,1]],np.r_[end[:,1],end[:,0]])),shape=(count,count)).tocsr()
    _,components=connected_components(adjacency,directed=False)
    for component in np.unique(components):
        ids=np.flatnonzero(components==component)
        if not any(int(i)in locks for i in ids):locks.update((int(i),np.zeros(3))for i in ids)
    force=np.zeros((count,3));difference=offset[edges[:,0]]-offset[edges[:,1]];np.add.at(force,end[:,0],weights[:,None]*difference);np.add.at(force,end[:,1],-weights[:,None]*difference)
    known=np.array(sorted(locks));free=np.setdiff1d(np.arange(count),known);delta=np.zeros((count,3));delta[known]=[locks[int(i)]for i in known];L=diags(np.asarray(adjacency.sum(1)).ravel())-adjacency;H=L[free][:,free].tocsc();rhs=-(L[free][:,known]@delta[known])-force[free]
    if len(free):delta[free]=spsolve(H,rhs)
    if not np.isfinite(delta).all():raise ValueError('Nonfinite source-cell transport')
    result=v+offset+delta[cell]
    return result,dict(targetedCells=len(targets),lockedCells=len(locked_cells),freeSourceCells=len(free),linearResidual=float(np.max(abs(H@delta[free]-rhs),initial=0)),maximumMovement_m=float(np.linalg.norm(result-v,axis=1).max()))
