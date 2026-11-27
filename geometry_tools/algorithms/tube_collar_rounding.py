"""Round a source tube's inlet collar while retaining its actual sections.

A quintic Hermite centreline has prescribed endpoint tangents and zero endpoint
acceleration. Each existing section is transported by a proper rigid rotation:
its polygonal area/calibre is retained, rather than shrinking the tube to fit a
tight turn. The inlet and every control at/after the last collar station stay
literal. No anatomical identity, safe clearance or acceptance is inferred.
"""
import numpy as np
from geometry_tools.algorithms.tube_station_topology import recover_station_rings


def _area(points):
    p=points-points.mean(0)
    return float(np.linalg.norm(np.cross(p,np.roll(p,-1,axis=0)).sum(0))/2)


def _rotation(a,b):
    a=a/np.linalg.norm(a);b=b/np.linalg.norm(b);axis=np.cross(a,b);cos=float(a@b)
    if cos<=-1+1e-10:raise ValueError('Antiparallel collar section needs an explicit roll reference')
    cross=np.array([[0.,-axis[2],axis[1]],[axis[2],0.,-axis[0]],[-axis[1],axis[0],0.]])
    return np.eye(3)+cross+cross@cross/(1+cos)


def round_station_collar(vertices,faces,inlet,station_count,*,origin,maximum_motion_m):
    v=np.asarray(vertices,float);f=np.asarray(faces,int);inlet=set(map(int,inlet));groups=recover_station_rings(v,f)
    if len(groups)!=1:raise ValueError('A single unbranched native tube is required')
    rings=list(groups[0]['rings'])
    if set(map(int,rings[-1]))==inlet:rings=rings[::-1]
    if set(map(int,rings[0]))!=inlet:raise ValueError('The specified inlet is not a source mouth')
    if not isinstance(station_count,int)or not 3<=station_count<len(rings)-1:raise ValueError('A bounded collar and unchanged downstream station are required')
    centers=np.array([v[q].mean(0)for q in rings]);normals=[]
    for j,ring in enumerate(rings[:station_count+1]):
        _,_,axes=np.linalg.svd(v[ring]-centers[j]);n=axes[-1]
        forward=centers[min(j+1,len(centers)-1)]-centers[max(j-1,0)]
        if n@forward<0:n=-n
        normals.append(n)
    start,stop=centers[0],centers[station_count];length=float(np.linalg.norm(np.diff(centers[:station_count+1],axis=0),axis=1).sum());d0=normals[0]*length;d1=normals[-1]*length;t=np.linspace(0.,1.,station_count+1)
    h0=1-10*t**3+15*t**4-6*t**5;h1=10*t**3-15*t**4+6*t**5;g0=t-6*t**3+8*t**4-3*t**5;g1=-4*t**3+7*t**4-3*t**5
    curve=h0[:,None]*start+h1[:,None]*stop+g0[:,None]*d0+g1[:,None]*d1
    dh0=-30*t**2+60*t**3-30*t**4;dh1=-dh0;dg0=1-18*t**2+32*t**3-15*t**4;dg1=-12*t**2+28*t**3-15*t**4
    tangent=dh0[:,None]*start+dh1[:,None]*stop+dg0[:,None]*d0+dg1[:,None]*d1;tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    output=v.copy();details=[];rotations=[];origin=np.asarray(origin,float)
    for j,ring in enumerate(rings[:station_count+1]):
        R=_rotation(normals[j],tangent[j]);rotations.append(R)
        if j in [0,station_count]:continue
        unrounded=(v[ring]-centers[j])@R.T+curve[j];stored=(unrounded-origin).astype('<f4').astype(float)+origin;output[ring]=stored;before=_area(v[ring]);after=_area(stored);delta=float(np.linalg.norm(unrounded-stored,axis=1).max());perimeter=float(np.linalg.norm(v[ring]-np.roll(v[ring],-1,axis=0),axis=1).sum());bound=2*delta*perimeter+len(ring)*delta**2+256*np.finfo(float).eps*before
        if abs(after-before)>bound:raise ValueError('Stored section area exceeds its explicit rounding bound')
        details.append(dict(station=j,areaBeforeMetersSquared=before,areaAfterMetersSquared=after,areaDifferenceMetersSquared=after-before,storageAreaBoundMetersSquared=bound,properRotationDeterminant=float(np.linalg.det(R))))
    protected=np.ones(len(v),bool)
    for ring in rings[1:station_count]:protected[ring]=False
    if not np.array_equal(output[protected],v[protected]):raise ValueError('A protected source section moved')
    motion=np.linalg.norm(output-v,axis=1)
    if motion.max()>maximum_motion_m:raise ValueError('Collar exceeds the caller\'s explicit motion domain')
    return output,dict(stations=station_count,changedControls=int(np.count_nonzero(motion)),maximumMotionMeters=float(motion.max()),maximumAllowedMotionMeters=float(maximum_motion_m),sourceInletAndDownstreamLiteral=True,sectionAreaRetainedWithinStorageRounding=True,sections=details,centerline='Quintic Hermite, prescribed source section normals at both ends, zero endpoint acceleration',anatomicallyAccepted=False)
