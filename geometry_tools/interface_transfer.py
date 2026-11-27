"""Literal ownership transfer of a previously verified tissue interface.

Partitioning can change the two owning part IDs without moving the seam. The
historical declaration remains intact; only the proposed catalogue is checked
against the explicitly named new owners.
"""
import numpy as np
from .common import identity
from .constraints import check_constraint


def _sheet(mesh,faces):
    ids=np.asarray(faces,int)
    if ids.ndim!=1 or not len(ids) or len(np.unique(ids))!=len(ids) or ids.min()<0 or ids.max()>=len(mesh.faces):
        raise ValueError('A complete, nonduplicate interface face list is required')
    normals=mesh.attributes.get('normals')
    if normals is None:raise ValueError('Stored interface corner normal words are required')
    signs=mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices),int))
    words=normals.astype(np.int64)*signs[:,None]
    result=[]
    for face in mesh.faces[ids]:
        corners=[tuple(mesh.vertices[i])+tuple(words[i])for i in face]
        result.append(min(tuple(corners[i:]+corners[:i])for i in range(3)))
    return sorted(result)


def verify_transfer(source,first,second,new_first,new_second,parameters):
    if source.get('kind')!='shared-interface' or source.get('parameters',{}).get('verification_pending'):
        raise ValueError('Only a verified shared tissue interface can transfer ownership')
    if len({m.frame for m in (first,second,new_first,new_second)})!=1:
        raise ValueError('An ownership transfer cannot change coordinate frames')
    old_check=check_constraint(source,first,second)
    if not old_check['passed']:raise ValueError('The historical source interface does not pass')
    for role,old,new in [('first',first,new_first),('second',second,new_second)]:
        if _sheet(old,source['parameters'][role+'_faces'])!=_sheet(new,parameters[role+'_faces']):
            raise ValueError('An ownership transfer must preserve every oriented source facet and stored normal word')
    new_check=check_constraint(dict(kind='shared-interface',parameters=parameters),new_first,new_second)
    if not new_check['passed']:raise ValueError('The new interface owners do not form a reciprocal seam')
    return dict(passed=True,source_check=old_check,candidate_check=new_check,
                source_parameters_sha256=identity(source['parameters']),
                source_geometry_sha256=[first.geometry_sha,second.geometry_sha],
                candidate_geometry_sha256=[new_first.geometry_sha,new_second.geometry_sha],
                all_oriented_facets_and_corner_words_literal=True,historical_constraint_unchanged=True)
