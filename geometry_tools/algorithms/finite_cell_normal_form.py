"""Canonical material ownership for an already encoded finite-cell boundary.

Exact opposed faces can cancel across different charts of the same cell when
an intermediate film has zero stored thickness. Remaining physical faces must
have one exterior owner or exactly two opposite owners. This performs no
coordinate welding by distance, hole filling, face projection or remeshing.
"""
from collections import defaultdict,Counter
import numpy as np


def normalize_cell_boundaries(positions,surfaces,cells,corner_words,corner_uvs):
    v=np.asarray(positions,float)
    owners=defaultdict(list)
    for pid,uses in cells.items():
        for use in uses:owners[use['surface']].append((pid,int(use['side'])))
    physical={}
    zero=0
    for sid,faces in surfaces.items():
        for index,face in enumerate(np.asarray(faces,int)):
            key=tuple(sorted(map(int,face)))
            p=v[face]
            if np.all(np.cross(p[1]-p[0],p[2]-p[0])==0):
                zero+=1
                continue
            q=tuple(map(int,face))
            direction=1 if q in (key,key[1:]+key[:1],key[2:]+key[:2])else-1
            order=[q.index(i)for i in key]
            row=physical.setdefault(key,dict(coefficients=Counter(),attributes=[]))
            for pid,side in owners[sid]:row['coefficients'][pid]+=side*direction
            words=(np.asarray(corner_words[sid][index])[order].astype(np.int32)*direction).astype(np.int16)
            uv=np.asarray(corner_uvs[sid][index])[order]
            row['attributes'].append((sid,index,words,uv))
    grouped={}
    cancelled=Counter()
    ownerless_cancelled=0
    bad_fields=[]
    for face,row in physical.items():
        live={pid:n for pid,n in row['coefficients'].items()if n}
        for pid,n in row['coefficients'].items():
            if not n:cancelled[pid]+=1
        if not live:
            # An exactly opposed free-face pair has zero coefficient for
            # every material. Dropping it preserves the complete boundary
            # chain, just as cancelling an internal zero-thickness film does.
            ownerless_cancelled+=1
            continue
        if any(abs(n)!=1 for n in live.values()):
            raise ValueError('Repeated positive-area material coverage')
        if len(live)not in (1,2):
            raise ValueError('Physical face does not have one or two material owners')
        if len(live)==2 and sum(live.values())!=0:
            raise ValueError('Reciprocal material owners have the same orientation')
        first=sorted(live)[0]
        second=next((pid for pid in live if pid!=first),None)
        direction=live[first]
        geo=np.cross(v[face[1]]-v[face[0]],v[face[2]]-v[face[0]])
        valid=[a for a in row['attributes']if np.all(a[2].astype(float)@geo>0)]
        sid,index,words,uv=(valid or row['attributes'])[0]
        if not valid:bad_fields.append(dict(face=list(face),source_surface=sid,source_face=index))
        order=[0,1,2]if direction==1 else[0,2,1]
        f=np.asarray(face,int)[order]
        w=(words[order].astype(np.int32)*direction).astype(np.int16)
        pair=(first,second)
        group=grouped.setdefault(pair,dict(faces=[],words=[],uvs=[],sources=[]))
        group['faces'].append(f);group['words'].append(w);group['uvs'].append(uv[order])
        group['sources'].append((sid,int(index)))
    out_surfaces,out_words,out_uvs,uses={},{},{},{pid:[]for pid in cells}
    ancestry={}
    for i,((first,second),g)in enumerate(sorted(grouped.items(),key=lambda item:repr(item[0]))):
        sid='normal-form:'+str(i)
        out_surfaces[sid]=np.asarray(g['faces'],np.int64)
        out_words[sid]=np.asarray(g['words'],np.int16)
        out_uvs[sid]=np.asarray(g['uvs'],np.float32)
        ancestry[sid]=g['sources']
        uses[first].append(dict(surface=sid,side=1,**({'neighbor':second}if second else{})))
        if second:uses[second].append(dict(surface=sid,side=-1,neighbor=first))
    return dict(surfaces=out_surfaces,corner_words=out_words,corner_uvs=out_uvs,cells=uses,
        source_faces=ancestry,evidence=dict(coordinates_moved=False,new_controls=0,
            exact_zero_area_faces_removed=zero,opposed_cell_boundary_pairs_cancelled=dict(cancelled),
            ownerless_physical_faces_cancelled=ownerless_cancelled,
            each_material_oriented_triangle_coefficient_preserved=True,
            positive_volume_material_not_removed=True,physical_faces=len(physical),
            output_charts=len(out_surfaces),outward_normal_field_failures=bad_fields,
            closed_cell_or_anatomical_acceptance_inferred=False))
