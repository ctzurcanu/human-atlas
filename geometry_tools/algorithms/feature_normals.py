"""Angle-weighted closest-feature pseudonormals on literal triangle geometry.

Independent NumPy implementation of Baerentzen & Aanaes (2005),
doi:10.1109/TVCG.2005.49. Closed, consistently oriented meshes admit the
paper's signed-distance classification. Open surfaces provide only a local
oriented side; boundary features are reported explicitly, never filled.
"""
import numpy as np


class FeatureNormals:
    def __init__(self,mesh):
        v,physical=np.unique(mesh.vertices,axis=0,return_inverse=True)
        f=physical[mesh.faces];tri=v[f]
        cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);length=np.linalg.norm(cross,axis=1)
        if np.any(length<1e-20):raise ValueError('Degenerate source triangle has no oriented feature normal')
        self.face_normal=cross/length[:,None];self.faces=f
        vertex=np.zeros_like(v)
        for k in range(3):
            a=tri[:,(k+1)%3]-tri[:,k];b=tri[:,(k+2)%3]-tri[:,k]
            angle=np.arctan2(np.linalg.norm(np.cross(a,b),axis=1),np.einsum('ij,ij->i',a,b))
            np.add.at(vertex,f[:,k],self.face_normal*angle[:,None])
        self.vertex_normal=vertex
        edges=np.sort(f[:,[1,2,2,0,0,1]].reshape(-1,2),axis=1)
        unique,ids,counts=np.unique(edges,axis=0,return_inverse=True,return_counts=True)
        self.face_edges=ids.reshape(-1,3);self.edge_counts=counts
        edge_normal=np.zeros((len(unique),3));np.add.at(edge_normal,ids,np.repeat(self.face_normal,3,axis=0))
        self.edge_normal=edge_normal;boundary_vertices=np.zeros(len(v),bool)
        boundary_vertices[unique[counts==1].ravel()]=True;self.boundary_vertices=boundary_vertices

    def _raw_at(self,hit):
        face=np.asarray(hit['face_ids'],int);bary=np.asarray(hit['barycentric'],float)
        if bary.shape!=(len(face),3):raise ValueError('Actual closest triangle IDs and barycentrics required')
        zero=bary<=1e-8;number=zero.sum(1);normal=self.face_normal[face].copy();boundary=np.zeros(len(face),bool)
        vertices=np.flatnonzero(number>=2)
        if len(vertices):
            ids=self.faces[face[vertices],bary[vertices].argmax(1)]
            normal[vertices]=self.vertex_normal[ids];boundary[vertices]=self.boundary_vertices[ids]
        edges=np.flatnonzero(number==1)
        if len(edges):
            ids=self.face_edges[face[edges],bary[edges].argmin(1)]
            normal[edges]=self.edge_normal[ids];boundary[edges]=self.edge_counts[ids]==1
        return normal,boundary

    def at(self,hit):
        normal,boundary=self._raw_at(hit)
        length=np.linalg.norm(normal,axis=1)
        if np.any(length<1e-12):raise ValueError('An ambiguous/nonorientable closest source feature needs review')
        return normal/length[:,None],boundary

    def at_masked(self,hit):
        """Return explicit invalid-feature masks without inventing a signed side.

        Strict callers retain ``at``. Ambiguous closest features have a zero
        normal, valid=False and boundary/unknown=True, so a masked distance
        consumer cannot count them as verified inside or outside samples.
        """
        normal,boundary=self._raw_at(hit)
        length=np.linalg.norm(normal,axis=1)
        valid=np.isfinite(normal).all(axis=1)&(length>=1e-12)
        unit=np.divide(normal,length[:,None],out=np.zeros_like(normal),where=valid[:,None])
        return unit,boundary|~valid,valid


def indexed_feature_normals(index,hit):
    if not hasattr(index,'_source_feature_normals'):index._source_feature_normals=FeatureNormals(index.mesh)
    return index._source_feature_normals.at(hit)


def metric_oriented_rows(points,hit,oriented_normal,*,unsigned_boundary=None):
    """Exact closest-distance magnitude and its local signed metric gradient."""
    delta=np.asarray(points,float)-hit['points'];distance=np.linalg.norm(delta,axis=1)
    polarity=np.where(np.einsum('ij,ij->i',delta,oriented_normal)>=0,1.,-1.)
    if unsigned_boundary is not None:polarity[np.asarray(unsigned_boundary,bool)]=1.
    normal=np.divide(delta,distance[:,None],out=oriented_normal.copy(),where=distance[:,None]>1e-12)
    normal*=polarity[:,None]
    return normal,distance*polarity


def classified_oriented_distance(points,hit,oriented_normal,boundary):
    """Local sheet-side distance, with unknown sign at actual open rims.

    An unsigned rim clearance is useful for avoiding a sheet edge, but its
    positive magnitude must never be counted as outside a body envelope.
    Neither a known local side nor this mask establishes a closed interior.
    """
    boundary=np.asarray(boundary,bool)
    _,distance=metric_oriented_rows(points,hit,oriented_normal)
    if boundary.shape!=distance.shape:raise ValueError('One open-rim flag per actual closest feature required')
    distance=distance.copy();distance[boundary]=np.nan
    return distance,~boundary


def inward_sheet_clearance(points,hit,oriented_normal,boundary):
    """Interior-side clearance except at rims, which have unsigned clearance."""
    normal,gap=metric_oriented_rows(points,hit,oriented_normal,unsigned_boundary=boundary)
    boundary=np.asarray(boundary,bool)
    inward=-normal;interior_gap=-gap
    inward[boundary]=normal[boundary];interior_gap[boundary]=gap[boundary]
    return inward,interior_gap,~boundary
