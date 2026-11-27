"""Transport declared complete source sections to a fixed recipient mouth.

The section centres receive a cubic falloff of the measured endpoint shift.
Every interior section is translated rigidly, preserving its literal shape
and area; only the explicitly corresponding inlet adopts the recipient curve.
No candidate section, anatomical attachment or recipient is inferred here.
"""
import numpy as np

from ..model import Mesh
from .surface_normals import transport_vertex_words


def section_area(p):
    p=np.asarray(p,float);c=p.mean(0)
    return float(np.linalg.norm(np.cross(p-c,np.roll(p,-1,axis=0)-c).sum(0))*.5)


def transport_source_collars(mesh, declarations, *, maximum_motion_m,
                             maximum_mouth_area_change,
                             maximum_facet_rotation_degrees):
    p=mesh.vertices;new=p.copy();owned=set();proof=[]
    for d in declarations:
        rings=[np.asarray(r,int)for r in d['sections']]
        target=np.asarray(d['target'],float)
        if len(rings)<3 or target.shape!=(len(rings[0]),3):
            raise ValueError('Three or more actual complete sections and a corresponding mouth required')
        if any(len(np.unique(r))!=len(r)or len(r)<3 for r in rings):
            raise ValueError('Every source section must be a non-repeating complete curve')
        points=[p[r]for r in rings];centres=np.array([q.mean(0)for q in points])
        arc=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(centres,axis=0),axis=1))]
        if np.any(np.diff(arc)<=0):raise ValueError('Source sections must have distinct centres')
        delta=target.mean(0)-centres[0]
        u=arc/arc[-1];weight=1-3*u*u+2*u*u*u
        weight[-1]=0.
        destinations=[q+delta*w for q,w in zip(points,weight)]
        destinations[0]=target
        destinations[-1]=points[-1].copy()
        change=abs(section_area(target)/section_area(points[0])-1)
        if change>maximum_mouth_area_change:
            raise ValueError('Recipient curve exceeds the physical source-mouth area bound')
        controls=[]
        for q,after in zip(points,destinations):
            for point,t in zip(q,after):
                ids=np.flatnonzero(np.all(p==point,axis=1))
                if owned.intersection(map(int,ids)):
                    raise ValueError('Declared collars share physical control ownership')
                owned.update(map(int,ids));controls.extend(map(int,ids));new[ids]=t
        for i in range(1,len(points)):
            np.testing.assert_allclose(destinations[i]-destinations[i][0],points[i]-points[i][0],atol=4e-16,rtol=0)
        proof.append(dict(sections=[r.tolist()for r in rings],sourceAreas_m2=[section_area(q)for q in points],
            candidateAreas_m2=[section_area(q)for q in destinations],centreWeights=weight.tolist(),
            maximumMouthAreaChange=change,fixedLastSection=True,interiorSectionShapesPreserved=True))
    motion=float(np.linalg.norm(new-p,axis=1).max(initial=0))
    if motion>maximum_motion_m:raise ValueError('Whole source collar exceeds its declared displacement budget')
    attrs={k:v.copy()for k,v in mesh.attributes.items()}
    attrs['normals'],normal_proof=transport_vertex_words(p,new,mesh.faces,attrs['normals'],
        maximum_facet_rotation_degrees=maximum_facet_rotation_degrees)
    return Mesh(new,mesh.faces,mesh.frame,attrs,mesh.source),dict(complete=True,collars=proof,
        maximumMotion_m=motion,normalTransport=normal_proof,originalFacesRetained=True,
        lumenDisksCreated=0,anatomicalAttachmentEstablished=False)
