"""Reconcile explicitly corresponding open rims without caps or face deletion.

Callers own anatomical identification and the finite movement/area budgets.
No pair is discovered or accepted merely because two vessels are nearby.
"""
import numpy as np

from ..model import Mesh
from ..source_port_dependencies import check_source_mouth
from .surface_normals import transport_vertex_words


def reconcile_declared_mouths(meshes, pairs, *, maximum_gap_m,
                              maximum_motion_m, maximum_relative_area_change,
                              maximum_facet_rotation_degrees, maximum_fixed_rim_motion_m=None):
    area = lambda p: float(np.linalg.norm(np.cross(p-p.mean(0), np.roll(p,-1,axis=0)-p.mean(0)).sum(0))*.5)
    positions = {key:m.vertices.copy() for key,m in meshes.items()}
    used = {key:set() for key in meshes}
    targets, proofs = [], []
    for pair in pairs:
        names = [pair['first'], pair['second']]
        rings = [np.asarray(pair[k],int) for k in ['firstRing','secondRing']]
        first, second = [meshes[key] for key in names]
        prior = check_source_mouth(first,second,*rings,tolerance_m=maximum_gap_m,
                                  area_tolerance=maximum_relative_area_change)
        if not prior['passed'] or rings[0].shape != rings[1].shape:
            raise ValueError('Explicit complete corresponding source mouths must satisfy the declared input bounds')
        points = [m.vertices[ids] for m,ids in zip([first,second],rings)]
        fixed = pair.get('fixedOwner')
        if fixed is not None:
            if fixed not in names or maximum_fixed_rim_motion_m is None or maximum_fixed_rim_motion_m <= 0:
                raise ValueError('Explicit fixed owner and separate positive movement budget required')
            target = points[names.index(fixed)].copy()
        else:
            target = (points[0]+points[1])*.5
            center = target.mean(0)
            target = center+(target-center)*np.sqrt((area(points[0])+area(points[1]))*.5/area(target))
        changes = [abs(area(target)/area(p)-1) for p in points]
        motion = max(float(np.linalg.norm(target-p,axis=1).max()) for p in points)
        budget = maximum_fixed_rim_motion_m if fixed is not None else maximum_motion_m
        if motion > budget or max(changes) > maximum_relative_area_change:
            raise ValueError('Common source rim exceeds a declared movement or section-area bound')
        aliases = []
        for name, ring in zip(names,points):
            owned = []
            for point,destination in zip(ring,target):
                ids = np.flatnonzero(np.all(meshes[name].vertices==point,axis=1))
                if used[name].intersection(map(int,ids)):
                    raise ValueError('Declared ports have overlapping source control ownership')
                used[name].update(map(int,ids))
                positions[name][ids] = destination
                owned.append(ids)
            aliases.append(owned)
        targets.append((names,rings,aliases))
        proofs.append(dict(first=names[0],second=names[1],maximumSourceGap_m=prior['maximumGap_m'],
                           maximumMotion_m=motion,relativeSectionAreaChanges=changes,fixedOwner=fixed,maximumAllowedMotion_m=budget))
    attributes, transport = {}, {}
    for name,m in meshes.items():
        attributes[name] = {k:v.copy() for k,v in m.attributes.items()}
        attributes[name]['normals'],transport[name] = transport_vertex_words(
            m.vertices,positions[name],m.faces,m.attributes['normals'],
            maximum_facet_rotation_degrees=maximum_facet_rotation_degrees)
    for names,rings,aliases in targets:
        fields = [attributes[name]['normals'][ids].astype(float)*
                  attributes[name]['normalSigns'][ids,None] for name,ids in zip(names,rings)]
        summed = fields[0]+fields[1]
        if np.any(np.linalg.norm(summed,axis=1)<10000):
            raise ValueError('Corresponding source outward fields cancel at the declared rim')
        direction = summed/np.linalg.norm(summed,axis=1)[:,None]
        for name,owned in zip(names,aliases):
            for i,ids in enumerate(owned):
                attributes[name]['normals'][ids] = np.rint(direction[i]*32767*
                    attributes[name]['normalSigns'][ids,None]).astype('<i2')
    result = {name:Mesh(positions[name],m.faces,m.frame,attributes[name],m.source)
              for name,m in meshes.items()}
    for pair,proof in zip(pairs,proofs):
        check = check_source_mouth(result[pair['first']],result[pair['second']],
                                  pair['firstRing'],pair['secondRing'],tolerance_m=1e-12)
        if not check['passed']:
            raise ValueError('Reconciled complete source rims do not match')
        proof['idealCommonMouth'] = check
    return result,dict(complete=True,pairs=proofs,normalTransport=transport,
                        originalFacesRetained=True,lumenDisksCreated=0,
                        anatomicalAttachmentEstablished=False,anatomicalAcceptance=False)
