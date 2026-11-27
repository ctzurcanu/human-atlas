"""Conforming finite cells from explicit, piecewise-linear signed fields.

Every tetrahedron uses the same field order and edge intersection cache. Tissue
interfaces are stored once, with reciprocal ownership. This is a construction
algorithm, not an anatomical classifier: the caller supplies fields and labels.
Distances, spacing and the final Float32 replay must be independently checked.
"""
from collections import defaultdict
from itertools import permutations
import numpy as np


def has_occupied_volume(points):
    """Reject only truly lower-dimensional leaves, including exact contacts.

    A clipped contact polygon is not tissue. Very thin nonzero volumes must
    survive, so an uncertain Float64 determinant uses exact binary rationals
    rather than a thickness or volume threshold.
    """
    p=np.asarray(points,float)
    if len(p)<4:
        return False
    q=p-p[0];a=q[np.argmax(np.linalg.norm(q,axis=1))]
    normals=np.cross(a,q);b=int(np.argmax(np.linalg.norm(normals,axis=1)));n=normals[b]
    residual=np.abs(q@n)
    bound=128*np.finfo(float).eps*np.linalg.norm(a)*np.linalg.norm(q[b])*np.linalg.norm(q,axis=1)
    if np.any(residual>bound):
        return True
    from fractions import Fraction
    r=[[Fraction(float(v))for v in row]for row in p]
    r=[[x-y for x,y in zip(row,r[0])]for row in r]
    first=next((row for row in r if any(row)),None)
    if first is None:
        return False
    cross=lambda x,y:[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]
    n=next((v for row in r if any(v:=cross(first,row))),None)
    return n is not None and any(sum(x*y for x,y in zip(row,n))!=0 for row in r)


def grid_tetrahedra(lower, upper, spacing,*,domain_distance=None):
    lower=np.asarray(lower,float);upper=np.asarray(upper,float)
    shape=np.ceil((upper-lower)/spacing).astype(int)+1
    axes=[np.linspace(a,b,n)for a,b,n in zip(lower,upper,shape)]
    points=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3)
    base=np.arange(len(points)).reshape(shape)[:-1,:-1,:-1].ravel()
    strides=np.array([shape[1]*shape[2],shape[2],1]);tets=[]
    if domain_distance:
        distance=domain_distance(points);near=np.zeros(len(base),bool)
        for x in (0,1):
            for y in (0,1):
                for z in (0,1):near|=distance[base+np.dot([x,y,z],strides)]<spacing
        base=base[near]
    for order in permutations(range(3)):
        steps=np.r_[0,np.cumsum(strides[list(order)])]
        t=base[:,None]+steps
        q=points[t[:1]];det=np.linalg.det(q[:,1:]-q[:,:1])[0]
        if det<0:t=t[:,[0,2,1,3]]
        tets.append(t)
    tets=np.concatenate(tets)
    if domain_distance:
        used,inv=np.unique(tets,return_inverse=True);points=points[used];tets=inv.reshape(-1,4)
    return points,tets,dict(shape=shape.tolist(),sample_points=len(points),tetrahedra=len(tets),maximum_spacing_m=float(max(np.diff(a).max()for a in axes)))


class Arrangement:
    def __init__(self,points,values,classify,*,vertex_snap_m=0.,axis_planes=None):
        points=np.asarray(points,float);values=np.array(values,float,copy=True)
        self.axis_planes=axis_planes or {}
        if self.axis_planes and vertex_snap_m:raise ValueError('Exact axis planes cannot be combined with independent vertex snapping')
        for k,(axis,threshold,sign)in self.axis_planes.items():
            if not 0<=k<values.shape[1]or axis not in (0,1,2)or sign not in (-1,1)or not np.isfinite(threshold):
                raise ValueError('An axis plane requires a valid field, coordinate, finite threshold and orientation')
            values[:,k]=sign*(points[:,axis]-threshold)
        self.points=list(points);self.values=list(values)
        self.classify=classify;self.edge_points={};self.boundaries={};self.leaves=0
        self.vertex_snap_m=vertex_snap_m;self.point_index={};self.collapsed_faces=0;self.max_snap=0.
        self.facet_labels={};self.face_tags={};self.triangle_tags={};self.surface_source_faces={}
        self.zero_volume_leaves=0;self.repeated_owner_witnesses=[]
        if vertex_snap_m:
            for i,p in enumerate(self.points):self.point_index[tuple(np.rint(p/vertex_snap_m).astype('i8'))]=i

    def vertex(self,p,value):
        if self.vertex_snap_m:
            key=tuple(np.rint(p/self.vertex_snap_m).astype('i8'))
            if key in self.point_index:return self.point_index[key]
            snapped=np.array(key)*self.vertex_snap_m;self.max_snap=max(self.max_snap,float(np.linalg.norm(p-snapped)))
            self.point_index[key]=len(self.points);p=snapped
        if self.axis_planes:
            value=np.array(value,copy=True)
            for k,(axis,threshold,sign)in self.axis_planes.items():value[k]=sign*(p[axis]-threshold)
        result=len(self.points);self.points.append(p);self.values.append(value);return result

    def intersection(self,a,b,k):
        if self.values[a][k]==0:return a
        if self.values[b][k]==0:return b
        key=(min(a,b),max(a,b),k)
        if key not in self.edge_points:
            a,b=key[:2];va,vb=self.values[a],self.values[b]
            t=va[k]/(va[k]-vb[k]);v=va*(1-t)+vb*t;v[k]=0.
            point=self.points[a]*(1-t)+self.points[b]*t
            if k in self.axis_planes:
                axis,threshold,_=self.axis_planes[k];point[axis]=threshold
            self.edge_points[key]=self.vertex(point,v)
            # Incidences closer than the declared construction resolution use
            # the same control; enforce the newly incident affine plane.
            self.values[self.edge_points[key]][k]=0.
        return self.edge_points[key]

    def split(self,faces,k):
        halves=[[],[]];cut=set()
        for poly in faces:
            tag=self.face_tags.get(tuple(sorted(poly)))
            for positive in (False,True):
                result=[]
                for a,b in zip(poly,poly[1:]+poly[:1]):
                    da,db=self.values[a][k],self.values[b][k]
                    ia=da>=0 if positive else da<=0;ib=db>=0 if positive else db<=0
                    if ia:result.append(a)
                    if ia!=ib:
                        q=self.intersection(a,b,k);result.append(q);cut.add(q)
                result=list(dict.fromkeys(result))
                if len(result)>=3:
                    halves[positive].append(result)
                    if tag is not None:self.face_tags[tuple(sorted(result))]=tag
        if len(cut)>=3:
            ids=sorted(cut);v=np.array([self.points[i]for i in ids]);center=v.mean(0)
            # The affine field gradient is constant on this tetrahedron.
            q=v-center;normal=self.gradients[:,k];normal=normal/np.linalg.norm(normal)
            u=np.cross(normal,np.eye(3)[np.argmin(abs(normal))]);u/=np.linalg.norm(u)
            w=np.cross(normal,u);angle=np.arctan2(q@w,q@u);loop=[ids[i]for i in np.argsort(angle)]
            halves[0].append(loop);halves[1].append(loop[::-1])
        return halves

    def emit(self,faces,label):
        self.leaves+=1
        if label is None:return
        ids=sorted({i for face in faces for i in face})
        if not has_occupied_volume([self.points[i]for i in ids]):
            self.zero_volume_leaves+=1
            return
        for polygon in faces:
            tag=self.face_tags.get(tuple(sorted(polygon)))
            # Identical shared polygons triangulate identically even when their
            # outward winding is opposite. Fan apex is a global vertex ID.
            polygon=list(dict.fromkeys(polygon))
            if len(polygon)<3:continue
            if len(polygon)==3:triangles=[tuple(polygon)]
            else:
                # Preserve collinear interface controls without zero-area fan
                # triangles. Convex ears are selected by the same global ID in
                # either winding; adjacent cells therefore agree exactly.
                triangles=[];remaining=polygon[:]
                while len(remaining)>3:
                    candidates=[]
                    for j,b in enumerate(remaining):
                        a,c=remaining[j-1],remaining[(j+1)%len(remaining)]
                        area=np.linalg.norm(np.cross(self.points[b]-self.points[a],self.points[c]-self.points[a]))
                        rest=[i for i in remaining if i!=b];q=np.array([self.points[i]for i in rest]);rest_area=np.linalg.norm(np.cross(q[1:-1]-q[0],q[2:]-q[0]).sum(0))
                        if area>1e-25 and rest_area>1e-25:candidates.append((b,j,(a,b,c)))
                    if not candidates:break
                    _,j,tri=min(candidates);triangles.append(tri);remaining.pop(j)
                if len(remaining)==3:triangles.append(tuple(remaining))
            for tri in triangles:
                q=np.array([self.points[i]for i in tri]);area=np.linalg.norm(np.cross(q[1]-q[0],q[2]-q[0]))
                if not area:self.collapsed_faces+=1;continue
                self.record(tri,label)
                if tag is not None:self.triangle_tags[tuple(sorted(tri))]=tag

    def record(self,tri,label):
        key=tuple(sorted(tri));sign=1 if tri in [key,(key[1],key[2],key[0]),(key[2],key[0],key[1])]else-1
        owners=self.boundaries.setdefault(key,{})
        owners[label]=owners.get(label,0)+sign
        if abs(owners[label])>1 and len(self.repeated_owner_witnesses)<16:
            self.repeated_owner_witnesses.append(dict(triangle=list(key),label=label,
                source_tetrahedron=getattr(self,'current_tet_index',None),signed_count=owners[label]))
        if owners[label]==0:del owners[label]
        if not owners:del self.boundaries[key]

    def build(self,tets,*,outer_field=0,skip_fields=(),progress=None):
        raw=np.asarray(self.values);tets=np.asarray(tets,int)
        tets=tets[np.min(raw[tets,outer_field],axis=1)<0]
        for number,tet in enumerate(tets):
            self.current_tet_index=number
            vals=raw[tet];active=[k for k in np.flatnonzero((vals.min(0)<0)&(vals.max(0)>0))if k not in skip_fields]
            # Ordered tetrahedron has positive signed volume.
            faces=[[int(tet[i])for i in row]for row in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)]]
            for f in faces:
                key=tuple(sorted(f))
                if key in self.facet_labels:self.face_tags[key]=self.facet_labels[key]
            if not active:
                label=self.classify(vals.mean(0))
                if label is not None:
                    for f in faces:
                        self.record(tuple(f),label);tag=self.face_tags.get(tuple(sorted(f)))
                        if tag is not None:self.triangle_tags[tuple(sorted(f))]=tag
                continue
            self.gradients=np.linalg.solve(np.asarray([self.points[i]for i in tet[1:]])-self.points[tet[0]],vals[1:]-vals[0])
            for k,(axis,_,sign)in self.axis_planes.items():
                self.gradients[:,k]=0.;self.gradients[axis,k]=sign
            work=[faces]
            for k in active:
                new=[]
                for poly in work:
                    ids=list({i for p in poly for i in p});s=np.array([self.values[i][k]for i in ids])
                    if s.min()<0<s.max():new.extend(q for q in self.split(poly,int(k))if q)
                    else:new.append(poly)
                work=new
            for poly in work:
                ids=list({i for p in poly for i in p});value=np.mean([self.values[i]for i in ids],axis=0)
                self.emit(poly,self.classify(value))
            if progress and number%5000==0:progress(number,len(tets),len(self.points),self.leaves)
        surfaces=defaultdict(list);cells=defaultdict(list)
        for tri,owners in self.boundaries.items():
            owners={n:s for n,s in owners.items()if s and n is not None}
            if not owners:continue
            if any(abs(s)!=1 for s in owners.values())or len(owners)>2:
                error=ValueError('A finite interface has repeated/nonreciprocal owners: '+str(owners))
                error.details=dict(triangle=list(tri),points_m=[np.asarray(self.points[i]).tolist()for i in tri],
                    owners=owners,zero_volume_leaves_skipped=self.zero_volume_leaves,
                    witnesses=self.repeated_owner_witnesses)
                raise error
            names=list(owners);a=names[0];b=names[1]if len(names)>1 else None
            if b is not None and owners[a]!=-owners[b]:raise ValueError('Tissue interface winding is not reciprocal')
            if owners[a]<0:tri=(tri[0],tri[2],tri[1])
            key=(a,b)if b is None or a<b else(b,a)
            surfaces[key].append(tri if a==key[0]else(tri[0],tri[2],tri[1]))
        result={}
        for n,((a,b),f)in enumerate(sorted(surfaces.items(),key=lambda q:str(q[0]))):
            sid='field-interface:'+str(n);result[sid]=np.asarray(f,int)
            self.surface_source_faces[sid]=np.array([self.triangle_tags.get(tuple(sorted(q)),-1)for q in f],int)
            cells[a].append(dict(surface=sid,side=1,neighbor=b))
            if b is not None:cells[b].append(dict(surface=sid,side=-1,neighbor=a))
        return np.asarray(self.points),result,dict(cells)
