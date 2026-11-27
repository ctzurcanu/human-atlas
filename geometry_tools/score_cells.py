"""Conforming score-comparison fields with explicit multi-field incidences.

Comparisons such as a-b, b-c and a-c are algebraically dependent. Independent
rounding of their intersection values creates fictitious, nearly flat cells.
This construction retains the underlying scores and inherited equality
groups. A contact constrained to a nonconstant field plane is not tissue.
This is numerical construction machinery, not an anatomical classifier.
"""
import numpy as np
from .implicit_cells import Arrangement


class ScoreArrangement(Arrangement):
    def __init__(self,points,scores,comparisons,classify,*,background_index,score_tolerance=1e-12,source_tie_tolerance=None):
        raw=np.asarray(scores,float);self.comparisons=[tuple(map(int,row))for row in comparisons]
        if raw.ndim!=2 or raw.shape[0]!=len(points)or not np.isfinite(raw).all():raise ValueError('Finite corresponding scores required')
        if not 0<=background_index<raw.shape[1]or score_tolerance<0:raise ValueError('Explicit background score and arithmetic tolerance required')
        if any(a==b or min(a,b)<0 or max(a,b)>=raw.shape[1]for a,b in self.comparisons):raise ValueError('Invalid comparison identity')
        self.score_count=raw.shape[1];self.background_index=int(background_index)
        self.background_value=float(raw[0,self.background_index]);self.score_tolerance=float(score_tolerance)
        self.source_tie_tolerance=self.score_tolerance if source_tie_tolerance is None else float(source_tie_tolerance)
        if self.source_tie_tolerance<self.score_tolerance or not np.isfinite(self.source_tie_tolerance):raise ValueError('Finite declared source tie tolerance required')
        if not np.all(raw[:,self.background_index]==self.background_value):raise ValueError('A single constant background score is required')
        self.maximum_source_projection=0.;self.maximum_derived_projection=0.;self.symbolic_contacts=0
        projected=[]
        for row in raw:
            ties=[(a,b)for a,b in self.comparisons if abs(row[a]-row[b])<=self.source_tie_tolerance]
            new=self.project(row,ties,maximum_change=4*self.source_tie_tolerance);self.maximum_source_projection=max(self.maximum_source_projection,float(abs(new-row).max()));projected.append(new)
        self.pair_offset=1+self.score_count
        values=np.array([self.fields(row)for row in projected])
        self.comparison_index={}
        for i,(a,b)in enumerate(self.comparisons):
            self.comparison_index[(a,b)]=(self.pair_offset+i,1.)
            self.comparison_index[(b,a)]=(self.pair_offset+i,-1.)
        super().__init__(points,values,classify)

    def project(self,scores,ties,maximum_change=None):
        scores=np.array(scores,float,copy=True);parent=list(range(self.score_count))
        def root(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        for a,b in ties:
            a,b=root(a),root(b)
            if a!=b:parent[max(a,b)]=min(a,b)
        groups={}
        for i in range(self.score_count):groups.setdefault(root(i),[]).append(i)
        before=scores.copy()
        for ids in groups.values():
            value=self.background_value if self.background_index in ids else scores[min(ids)]
            scores[ids]=value
        scores[self.background_index]=self.background_value
        bound=max(4*self.score_tolerance,1e-14)if maximum_change is None else max(maximum_change,1e-14)
        if abs(scores-before).max()>bound:raise ValueError('Equality incidence exceeds the declared score arithmetic bound')
        return scores

    def fields(self,scores):
        return np.r_[-1.,scores,[scores[a]-scores[b]for a,b in self.comparisons]]

    def intersection(self,a,b,k):
        if k<self.pair_offset:raise ValueError('Only declared score comparisons may cut this arrangement')
        if self.values[a][k]==0:return a
        if self.values[b][k]==0:return b
        key=(min(a,b),max(a,b),k)
        if key in self.edge_points:return self.edge_points[key]
        a,b=key[:2];va,vb=self.values[a],self.values[b];sa,sb=va[1:self.pair_offset],vb[1:self.pair_offset]
        weight=va[k]/(va[k]-vb[k]);scores=sa*(1-weight)+sb*weight
        inherited=[(i,j)for i,j in self.comparisons if sa[i]==sa[j]and sb[i]==sb[j]]
        inherited.append(self.comparisons[k-self.pair_offset]);projected=self.project(scores,inherited)
        self.maximum_derived_projection=max(self.maximum_derived_projection,float(abs(projected-scores).max()))
        point=self.points[a]*(1-weight)+self.points[b]*weight
        index=self.vertex(point,self.fields(projected));assert self.values[index][k]==0
        self.edge_points[key]=index;return index

    def emit(self,faces,label):
        ids=sorted({i for face in faces for i in face})
        if not ids:return
        values=np.array([self.values[i]for i in ids])
        for k in range(self.pair_offset,values.shape[1]):
            if np.all(values[:,k]==0.)and np.any(self.gradients[:,k]!=0.):
                self.symbolic_contacts+=1;self.zero_volume_leaves+=1;return
        super().emit(faces,label)

    def build(self,tets,**kwargs):
        if 'skip_fields'in kwargs:raise ValueError('Score fields are not geometric cut planes')
        return super().build(tets,skip_fields=range(1,self.pair_offset),**kwargs)

    def save_checkpoint(self,directory,metadata):
        """Save data only; incomplete stages are never installable cells."""
        from pathlib import Path
        from .common import atomic_npz,atomic_json,digest_file
        directory=Path(directory)
        if self.axis_planes or self.vertex_snap_m or self.facet_labels or self.face_tags or self.triangle_tags:
            raise ValueError('This score checkpoint supports untagged, unsnapped comparison fields only')
        labels=sorted({n for owners in self.boundaries.values()for n in owners});lookup={n:i for i,n in enumerate(labels)}
        faces=np.array(list(self.boundaries),int).reshape(-1,3);owners=np.full((len(faces),2),-1,int);signs=np.zeros((len(faces),2),int)
        for row,face in enumerate(faces):
            current=self.boundaries[tuple(face)]
            if len(current)>2 or any(abs(s)!=1 for s in current.values()):raise ValueError('A checkpoint has invalid finite-face ownership')
            for j,(name,sign)in enumerate(sorted(current.items())):owners[row,j]=lookup[name];signs[row,j]=sign
        edge_keys=np.array(list(self.edge_points),int).reshape(-1,3);edge_values=np.array(list(self.edge_points.values()),int)
        atomic_npz(directory/'state.npz',positions=np.asarray(self.points),values=np.asarray(self.values),
            boundary_faces=faces,boundary_owners=owners,boundary_signs=signs,edge_keys=edge_keys,edge_values=edge_values)
        atomic_json(directory/'checkpoint.json',dict(schema='score-arrangement-state/v1',metadata=metadata,labels=labels,
            state_sha256=digest_file(directory/'state.npz'),leaves=self.leaves,zero_volume_leaves=self.zero_volume_leaves,
            collapsed_faces=self.collapsed_faces,maximum_source_projection=self.maximum_source_projection,
            maximum_derived_projection=self.maximum_derived_projection,symbolic_contacts=self.symbolic_contacts,
            installed=False,not_an_installable_geometry_candidate=True))

    def restore_checkpoint(self,directory,expected_inputs):
        from pathlib import Path
        import json
        from .common import digest_file
        directory=Path(directory);record=json.loads((directory/'checkpoint.json').read_text())
        if record['schema']!='score-arrangement-state/v1'or record['metadata']['inputs']!=expected_inputs:
            raise ValueError('Checkpoint inputs or implementation changed')
        if record['state_sha256']!=digest_file(directory/'state.npz'):raise ValueError('Checkpoint data checksum failed')
        z=np.load(directory/'state.npz');points=z['positions'];values=z['values'];source_count=len(self.points)
        if not np.array_equal(points[:source_count],np.asarray(self.points))or not np.array_equal(values[:source_count],np.asarray(self.values)):
            raise ValueError('A checkpoint changed its original position or score controls')
        self.points=list(points);self.values=list(values);labels=record['labels'];self.boundaries={}
        for face,owners,signs in zip(z['boundary_faces'],z['boundary_owners'],z['boundary_signs']):
            self.boundaries[tuple(map(int,face))]={labels[i]:int(s)for i,s in zip(owners,signs)if i>=0}
        self.edge_points={tuple(map(int,key)):int(value)for key,value in zip(z['edge_keys'],z['edge_values'])}
        for name in ['leaves','zero_volume_leaves','collapsed_faces','maximum_source_projection','maximum_derived_projection','symbolic_contacts']:
            setattr(self,name,record[name])
        return record['metadata']
