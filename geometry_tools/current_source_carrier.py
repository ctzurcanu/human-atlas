"""Reconstitute a complete source from its current, active display partitions.

Suppressed parent buffers are retained provenance and may be stale after a
regional fit. Only literal source-ID/face lineage establishes this projection.
"""
import numpy as np
from geometry_tools.model import Mesh


def reconstruct_current_carrier(atlas, source, source_row, reader, *, catalogue_path):
    parts={p['id']:p for p in atlas['parts']}
    lookup={int(raw):i for i,raw in enumerate(source.source['source_vertex_ids'])}
    points=source.vertices.copy();attrs={k:v.copy() for k,v in source.attributes.items()}
    assigned=np.zeros(len(points),bool);face_seen=np.zeros(len(source.faces),bool)
    projections=[]
    for pid in source_row['activePartIds']:
        p=parts[pid]
        assert not p.get('suppressed') and pid not in atlas['defaultHidden']
        child=reader.read_part(atlas,p,frame=source.frame,catalogue_path=catalogue_path)
        ids=np.array([lookup[int(x)] for x in child.source['source_vertex_ids']])
        faces=p.get('sourceDerivedRegion',{}).get('sourceFaceIds',list(range(len(source.faces))))
        np.testing.assert_array_equal(ids[child.faces],source.faces[faces])
        same=assigned[ids]
        np.testing.assert_array_equal(points[ids[same]],child.vertices[same])
        for k in attrs:np.testing.assert_array_equal(attrs[k][ids[same]],child.attributes[k][same])
        points[ids]=child.vertices
        for k in attrs:attrs[k][ids]=child.attributes[k]
        assigned[ids]=True;face_seen[faces]=True
        projections.append(dict(partId=pid,carrierControlIds=ids.tolist(),sourceFaceIds=faces,
                                sourceGeometrySha256=child.geometry_sha,sourceAttributesSha256=child.attributes_sha))
    assert assigned.all() and face_seen.all(), 'Complete active source-face/control coverage required'
    return Mesh(points,source.faces,source.frame,attrs,source.source),projections
