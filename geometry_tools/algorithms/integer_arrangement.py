"""Canonical rational intersections for integer-valued affine material fields.

Equivalent intersection paths share one exact source-control weight key.
Coordinates are evaluated once, from the same immutable source controls.
The geometry and anatomical interpretation are supplied by the caller.
"""
from fractions import Fraction
import numpy as np
from geometry_tools.implicit_cells import Arrangement


class IntegerArrangement(Arrangement):
    def __init__(self,points,values,classify):
        raw=np.asarray(values,float)
        if not np.isfinite(raw).all()or not np.array_equal(raw,np.rint(raw))or np.max(abs(raw),initial=0)>2**30:
            raise ValueError('Exact integer fields within 30 bits are required')
        super().__init__(points,values,classify)
        self.original_points=np.array(points,float,copy=True)
        self.original_fields=raw.astype(np.int64)
        self.affine_weights=[((i,Fraction(1)),)for i in range(len(points))]
        self.affine_index={self.key(weights):i for i,weights in enumerate(self.affine_weights)}

    @staticmethod
    def key(weights):
        return tuple((i,w.numerator,w.denominator)for i,w in weights)

    def exact_value(self,weights,k):
        return sum((w*int(self.original_fields[i,k])for i,w in weights),Fraction(0))

    def intersection(self,a,b,k):
        edge=(min(a,b),max(a,b),k)
        if edge in self.edge_points:return self.edge_points[edge]
        wa,wb=self.affine_weights[a],self.affine_weights[b]
        va,vb=self.exact_value(wa,k),self.exact_value(wb,k)
        if va==0:return a
        if vb==0:return b
        t=va/(va-vb)
        if not 0<t<1:raise ValueError('Float clipping disagreed with the exact affine sign')
        combined={}
        for i,w in wa:combined[i]=combined.get(i,Fraction(0))+(1-t)*w
        for i,w in wb:combined[i]=combined.get(i,Fraction(0))+t*w
        weights=tuple((i,w)for i,w in sorted(combined.items())if w)
        assert sum((w for i,w in weights),Fraction(0))==1
        key=self.key(weights)
        if key not in self.affine_index:
            index=len(self.points);self.affine_index[key]=index
            point=sum((self.original_points[i]*float(w)for i,w in weights),np.zeros(3))
            value=np.array([float(self.exact_value(weights,j))for j in range(self.original_fields.shape[1])])
            self.points.append(point);self.values.append(value);self.affine_weights.append(weights)
        index=self.affine_index[key];self.edge_points[edge]=index
        return index

    def build_maximum(self,tetrahedra,region_count,*,comparison_offset=1,
                      priority_field=None,priority_label='priority'):
        """Clip each winning material cell, excluding irrelevant planes.

        Fields contain every score[i]-score[j], in combinations order.
        Intersections remain canonical across adjoining tetrahedra. A
        coincident score function uses the lower region index for bulk
        ownership; a zero-volume boundary is never emitted as tissue.

        An optional positive priority halfspace overrides material ownership.
        Clip it after the material subdivision so both sides of the priority
        interface share every material intersection. Its distance scale never
        changes the underlying material comparisons. Equivalent intersections
        retain the same exact rational source-control keys.
        """
        from itertools import combinations
        if region_count<2:raise ValueError('At least two material scores are required')
        columns={pair:i+comparison_offset for i,pair in enumerate(combinations(range(region_count),2))}
        if max(columns.values())>=self.original_fields.shape[1]:raise ValueError('Missing pairwise score fields')
        if priority_field is not None:
            if not isinstance(priority_field,(int,np.integer))or not 0<priority_field<self.original_fields.shape[1]or priority_field in columns.values():
                raise ValueError('Priority needs a separate signed field, outside the domain and material comparisons')
            if not isinstance(priority_label,str)or not priority_label or priority_label in {'region:'+str(i)for i in range(region_count)}:
                raise ValueError('Priority tissue needs its own distinct nonempty label')
        raw=np.asarray(self.values)
        for number,tet in enumerate(np.asarray(tetrahedra,int)):
            self.current_tet_index=number
            values=raw[tet]
            relative=np.c_[np.zeros(4),-values[:,comparison_offset:comparison_offset+region_count-1]]
            guaranteed=relative.min(0).max()
            possible=np.flatnonzero(relative.max(0)>=guaranteed)
            # Strict nodal dominance is an affine proof on the entire
            # tetrahedron. Excluding these never-winning functions changes
            # neither the final cells nor their shared-face subdivision.
            # Retain the complete path for coincident nodal scores.
            if np.all(values[:,comparison_offset:comparison_offset+len(columns)]!=0):
                winners=np.argmax(relative,axis=1)
                if np.all(winners==winners[0]):
                    possible=np.array([winners[0]])
                else:
                    dominance=np.all(relative[:,:,None]>relative[:,None,:],axis=0)
                    possible=np.array([r for r in possible if not dominance[:,r].any()],int)
            original=[[int(tet[i])for i in row]for row in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)]]
            for face in original:
                key=tuple(sorted(face))
                if key in self.facet_labels:self.face_tags[key]=self.facet_labels[key]
            gradients=None
            for region in possible:
                if any(np.array_equal(relative[:,region],relative[:,other])for other in range(region)):
                    continue
                work=original
                for other in possible:
                    if other==region:continue
                    k=columns[tuple(sorted((int(region),int(other))))];positive=region<other
                    ids=sorted({i for face in work for i in face})
                    signed=np.array([self.values[i][k]for i in ids])*(1 if positive else-1)
                    if signed.max()<=0 and signed.min()<0:
                        work=[];break
                    if signed.min()>=0:continue
                    if gradients is None:
                        gradients=np.linalg.solve(np.asarray([self.points[i]for i in tet[1:]])-self.points[tet[0]],values[1:]-values[0])
                    self.gradients=gradients
                    work=self.split(work,k)[1 if positive else 0]
                    if not work:break
                if not work:continue
                label='region:'+str(int(region))
                if priority_field is None:
                    self.emit(work,label)
                    continue
                ids=sorted({i for face in work for i in face})
                signed=np.array([self.values[i][priority_field]for i in ids])
                if signed.min()>=0:
                    self.emit(work,priority_label)
                elif signed.max()<=0:
                    self.emit(work,label)
                else:
                    if gradients is None:
                        gradients=np.linalg.solve(np.asarray([self.points[i]for i in tet[1:]])-self.points[tet[0]],values[1:]-values[0])
                    self.gradients=gradients
                    negative,positive=self.split(work,priority_field)
                    if negative:self.emit(negative,label)
                    if positive:self.emit(positive,priority_label)
        return super().build(np.empty((0,4),int))
