"""Replay the saved, root-preserving local hepatic course field."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial.distance import cdist
from .model import AtlasReader
from .open3d_backend import library

class HepaticCourseField:
    def __init__(self,atlas,field_path,anchors_path):
        z=np.load(field_path);self.controls=z['controls'];self.coefficients=z['coefficients'];self.sigma=float(z['sigma'])
        self.anchors=json.loads(Path(anchors_path).read_text())['anchors'];reader=AtlasReader();p=next(p for p in atlas['parts']if p['id']=='TA98HEPATIC:male:region:A12.3.09.001')
        m=reader.read_part(atlas,p,frame='model:local-ta98:assembled-metres');self.o3d=library();self.cava=self.o3d.t.geometry.RaycastingScene(nthreads=1)
        # Native queries use Float32. Recenter before converting so the 1.3 m
        # body height cannot introduce submicrometre jitter into thin wall cells.
        self.query_origin=m.vertices.mean(axis=0)
        self.cava.add_triangles(self.o3d.core.Tensor((m.vertices-self.query_origin).astype(np.float32)),self.o3d.core.Tensor(m.faces.astype(np.uint32)))
    @staticmethod
    def smooth(t):
        t=np.clip(t,0,1);return t*t*(3-2*t)
    def weight(self,p):
        distance=self.cava.compute_distance(self.o3d.core.Tensor((np.asarray(p,float)-self.query_origin).astype(np.float32)),nthreads=1).numpy().astype(float)
        w=self.smooth((distance-.003)/.009)
        for a in self.anchors.values():w*=self.smooth((np.linalg.norm(p-np.asarray(a['center_m']),axis=1)-a['radius_m'])/.009)
        return w
    def __call__(self,p):
        p=np.asarray(p,float);return p+self.weight(p)[:,None]*(np.exp(-.5*cdist(p,self.controls,'sqeuclidean')/self.sigma**2)@self.coefficients)
    def jacobian(self,p,step=1e-6):
        p=np.asarray(p,float);result=np.empty((len(p),3,3))
        for axis in range(3):
            a=p.copy();b=p.copy();a[:,axis]+=step;b[:,axis]-=step;result[:,:,axis]=(self(a)-self(b))/(2*step)
        return result


class HepaticVeinShear:
    """A globally invertible, unit-determinant lift; the caval side is fixed."""
    def __init__(self,record):
        if record.get('schema')!='left-hepatic-vein-shear/v1':raise ValueError('Unexpected hepatic shear schema')
        self.start=float(record['startX']);self.end=float(record['endX']);self.shift=np.array([0.,*record['shiftYZMetres']])
        self.start_z=float(record.get('startZ',0.));self.end_z=float(record.get('endZ',1.));self.anterior_y=float(record.get('anteriorShiftYMetres',0.))
        if self.end<=self.start or not np.isfinite(np.r_[self.start,self.end,self.shift]).all():raise ValueError('Invalid finite shear controls')
    def weight(self,p):
        t=np.clip((np.asarray(p)[:,0]-self.start)/(self.end-self.start),0.,1.);return t*t*(3-2*t)
    def __call__(self,p):
        p=np.asarray(p,float);q=p+self.weight(p)[:,None]*self.shift;t=np.clip((q[:,2]-self.start_z)/(self.end_z-self.start_z),0.,1.);q[:,1]+=self.anterior_y*t*t*(3-2*t);return q
    def jacobian(self,p,step=None):
        p=np.asarray(p,float);t=np.clip((p[:,0]-self.start)/(self.end-self.start),0.,1.);slope=6*t*(1-t)/(self.end-self.start)
        j=np.tile(np.eye(3),(len(p),1,1));j[:,:,0]+=slope[:,None]*self.shift
        z=p[:,2]+self.weight(p)*self.shift[2];t=np.clip((z-self.start_z)/(self.end_z-self.start_z),0.,1.);slope_z=self.anterior_y*6*t*(1-t)/(self.end_z-self.start_z);j[:,1,:]+=slope_z[:,None]*j[:,2,:];return j
