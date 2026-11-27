"""Boolean occupied-volume overlap and its source-bound surface-flux derivative.

Coordinates and control lengths must share units. Returned volume has length^3
units; its derivative has length^2 units. Numerical Float64 predicates and
analysis classifiers do not establish anatomical tissue identity or acceptance.
"""
import numpy as np
import manifold3d as md
from geometry_tools.distance import barycentric
from .finite_tissue_boolean import make_solid,provenance


def overlap_volume_flux(vertices,faces,vertex_operator,fixed,*,original_id):
    v=np.asarray(vertices,float);f=np.asarray(faces,int);operator=np.asarray(vertex_operator,float)
    if operator.ndim!=2 or len(operator)!=len(v)or not np.isfinite(operator).all():
        raise ValueError('One finite source-control operator row per moving vertex required')
    moving,_=make_solid(v,f,original_id=original_id);intersection=moving^fixed
    if intersection.status()!=md.Error.NoError:raise ValueError('Boolean overlap did not finish with a valid manifold')
    volume=float(intersection.volume());gradient=np.zeros((operator.shape[1],2))
    if intersection.is_empty():return volume,gradient,dict(sourceBoundaryFragments=0,booleanStatus='NoError')
    output=intersection.to_mesh64();owner,face_ids=provenance(output);take=np.flatnonzero(owner==original_id)
    if not len(take):return volume,gradient,dict(sourceBoundaryFragments=0,booleanStatus='NoError')
    triangles=np.asarray(output.vert_properties)[:,:3][np.asarray(output.tri_verts)[take]];source_faces=face_ids[take].astype(int)
    if np.any(source_faces>=len(f)):raise ValueError('Boolean fragment has an invalid original moving-face provenance')
    weights,valid=barycentric(triangles.mean(1),v[f[source_faces]])
    if not valid.all()or np.min(weights)<-1e-7 or np.max(weights)>1+1e-7:raise ValueError('Boolean fragments do not replay on their original source faces')
    sample_operator=np.einsum('ni,nic->nc',weights,operator[f[source_faces]])
    area_normal=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])/2
    gradient[:,0]=sample_operator.T@area_normal[:,0];gradient[:,1]=sample_operator.T@area_normal[:,2]
    return volume,gradient,dict(sourceBoundaryFragments=len(take),booleanStatus='NoError',
        originalFaceProvenanceVerified=True,surfaceFluxDerivative=True,numericalFloat64=True)
