"""Preserve supplied surface roles through reciprocal Boolean subdivisions.

No anatomical role is inferred here. Charts are split by existing parent
labels; positions, triangles and reciprocal owners remain unchanged.
"""
import numpy as np


def split_by_source_role(result, source_roles):
    surfaces={};charts={};roles={};children={}
    for sid,faces in result['surfaces'].items():
        chart=result['charts'][sid]
        labels=np.array([source_roles.get(str(parent),'') for parent in chart['source_surfaces']],object)
        assert len(labels)==len(faces)
        unique=sorted(set(labels));children[sid]=[]
        for k,label in enumerate(unique):
            select=np.flatnonzero(labels==label)
            name=sid if len(unique)==1 else f'{sid}/source-role-{k}'
            children[sid].append(name);surfaces[name]=faces[select];roles[name]=str(label)
            charts[name]={key:(np.asarray(value)[select].tolist() if key=='source_surfaces' else np.asarray(value)[select]) for key,value in chart.items()}
    cells={pid:[dict(use,surface=child) for use in uses for child in children[use['surface']]] for pid,uses in result['cells'].items()}
    out=dict(result, surfaces=surfaces,charts=charts,cells=cells)
    return out,roles
