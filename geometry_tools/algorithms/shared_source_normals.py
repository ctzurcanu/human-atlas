"""Normal transport across caller-declared slices of one original vertex bank.

Shared identity is the original raw bank ID inside one explicitly declared
family.  World-coordinate coincidence never creates a family.  Temporary face
deduplication affects only the normal calculation, never a part's topology.
"""
import hashlib
import json
import time

import numpy as np

from .surface_normals import transport_vertex_words


def transport_shared_source_normals(meshes, positions, families, *,
                                    maximum_facet_rotation_degrees=65,
                                    wall_seconds=5.):
    """Return ``(normal_words_by_owner, proof)`` without modifying any input.

    Families must come from caller-verified identical original bank/layout
    signatures.  Where source records are supplied, the signature is rechecked.
    Unlisted owners are independent singleton families. Original signs, UVs,
    positions and faces are retained; this function emits only normal words.
    """
    started = time.monotonic()
    bound = float(maximum_facet_rotation_degrees)
    if (not np.isfinite(wall_seconds) or not 0 < wall_seconds <= 5
            or not np.isfinite(bound) or not 0 <= bound < 90
            or not isinstance(meshes, dict) or not 1 <= len(meshes) <= 128
            or not isinstance(positions, dict) or set(positions) != set(meshes)
            or not isinstance(families, (list, tuple)) or len(families) > 128):
        raise ValueError('Bounded explicit owners, corresponding positions, families and a facet bound below 90 degrees required')

    def deadline():
        if time.monotonic() - started >= wall_seconds:
            raise TimeoutError('Shared source-normal deadline; no complete field emitted')

    def signature(mesh):
        record = mesh.source.get('record')
        if record is None:
            return None
        defaults = dict(geometryOrigin=[0., 0., 0.], geometryScale=1., sourceOffset=[0., 0., 0.])
        keys = ['chunk', 'positions', 'normals', 'normalSigns', 'uvs', 'vertexCount',
                'geometryOrigin', 'geometryScale', 'sourceOffset', 'geometryInstance']
        return json.dumps(dict(bankPath=mesh.source.get('bank_path'),
                               bankSha256=mesh.source.get('bank_sha256'),
                               layout={key: record.get(key, defaults.get(key)) for key in keys}),
                          sort_keys=True, separators=(',', ':'))

    deadline()
    if sum(len(mesh.vertices) for mesh in meshes.values()) > 150000 or sum(len(mesh.faces) for mesh in meshes.values()) > 600000:
        raise ValueError('At most 150,000 source controls and 600,000 original faces supported')
    ordered, claimed = [], set()
    for family in families:
        if (not isinstance(family, (list, tuple)) or not family
                or len(set(family)) != len(family) or set(family) - set(meshes)
                or claimed.intersection(family)):
            raise ValueError('Every explicit family must have unique, existing, singly owned members')
        ordered.append(list(family))
        claimed.update(family)
    ordered.extend([[pid] for pid in sorted(set(meshes) - claimed)])
    result, evidence = {}, []
    for family in ordered:
        deadline()
        if len({meshes[pid].frame for pid in family}) != 1:
            raise ValueError('One original source frame required within a normal family')
        signatures = [signature(meshes[pid]) for pid in family]
        if any(value is not None for value in signatures) and len(set(signatures)) != 1:
            raise ValueError('Declared family original bank/layout signatures differ')
        raw_ids, old_rows, new_rows, normal_rows, sign_rows, uv_rows = [], [], [], [], [], []
        offsets, attribute_presence = {}, set()
        count = 0
        for pid in family:
            deadline()
            mesh = meshes[pid]
            q = np.asarray(positions[pid], dtype=float)
            words = np.asarray(mesh.attributes.get('normals'))
            signs = np.asarray(mesh.attributes.get('normalSigns', np.ones(len(mesh.vertices), dtype='i1')))
            if (not len(mesh.faces) or q.shape != mesh.vertices.shape or not np.isfinite(q).all()
                    or words.shape != mesh.vertices.shape or not np.issubdtype(words.dtype, np.integer)
                    or np.any(words.astype(float) < -32768) or np.any(words.astype(float) > 32767)
                    or signs.shape != (len(mesh.vertices),) or not np.isin(signs, [-1, 1]).all()):
                raise ValueError('Actual source triangles, finite corresponding poses and literal normal words/signs required')
            master = mesh.source.get('source_vertex_ids')
            if master is None:
                if len(family) > 1:
                    raise ValueError('Every shared family member requires original source_vertex_ids')
                master = np.arange(len(mesh.vertices))
            master = np.asarray(master)
            if (master.shape != (len(mesh.vertices),) or not np.issubdtype(master.dtype, np.integer)
                    or master.min(initial=0) < 0 or master.max(initial=0) > np.iinfo(np.int64).max):
                raise ValueError('Original raw bank IDs must cover every source control')
            presence = ('normalSigns' in mesh.attributes, 'uvs' in mesh.attributes)
            attribute_presence.add(presence)
            if 'uvs' in mesh.attributes:
                uv = np.asarray(mesh.attributes['uvs'])
                if uv.ndim != 2 or len(uv) != len(master) or not np.isfinite(uv).all():
                    raise ValueError('Finite original UV rows required')
                uv_rows.append(uv)
            offsets[pid] = slice(count, count + len(master))
            count += len(master)
            raw_ids.append(master.astype(np.int64, copy=False))
            old_rows.append(mesh.vertices)
            new_rows.append(q)
            normal_rows.append(words.astype('<i2', copy=False))
            sign_rows.append(signs)
        if len(attribute_presence) != 1:
            raise ValueError('Declared source-family normal-sign/UV attribute presence differs')
        master, representative, inverse = np.unique(np.concatenate(raw_ids), return_index=True, return_inverse=True)
        arrays = [('old coordinates', np.vstack(old_rows)), ('new coordinates', np.vstack(new_rows)),
                  ('normal words', np.vstack(normal_rows)), ('normal signs', np.concatenate(sign_rows))]
        if uv_rows:
            arrays.append(('UVs', np.vstack(uv_rows)))
        for name, values in arrays:
            if not np.array_equal(values, values[representative][inverse]):
                raise ValueError('A shared original raw ID has conflicting ' + name)
        old, new, words = [values[representative] for _, values in arrays[:3]]
        family_faces = np.vstack([inverse[offsets[pid]][meshes[pid].faces] for pid in family])
        # Cyclic permutations retain orientation. Reversed faces remain distinct.
        begin = np.argmin(family_faces, axis=1)
        canonical = np.take_along_axis(family_faces, (begin[:, None] + np.arange(3)) % 3, axis=1)
        normal_faces = np.unique(canonical, axis=0)
        deadline()
        transported, proof = transport_vertex_words(old, new, normal_faces, words,
                              maximum_facet_rotation_degrees=bound)
        deadline()
        if not proof['unchanged_stars_literal']:
            raise ValueError('Unchanged complete source stars did not retain literal normal words')
        for pid in family:
            result[pid] = transported[inverse[offsets[pid]]].copy()
        evidence.append(dict(owners=family, sharedOriginalBankIds=len(master),
            originalPartFaces=len(family_faces), normalComputationFaces=len(normal_faces),
            orientedDuplicateFacesCountedOnce=len(family_faces) - len(normal_faces),
            originalBankSignatureRechecked=signatures[0] is not None,
            originalBankSignatureSha256=hashlib.sha256(signatures[0].encode()).hexdigest() if signatures[0] is not None else None,
            transport=proof))
    deadline()
    return result, dict(complete=True,
        method='one complete original-bank face union per explicitly declared normal family',
        families=evidence, sourceRawControls=sum(len(mesh.vertices) for mesh in meshes.values()),
        sourceFaces=sum(len(mesh.faces) for mesh in meshes.values()),
        maximumFacetRotationDegrees=bound, commonOriginalRawIdsAgree=True,
        crossFamilyCoordinatesNeverMerged=True, normalSignsRetained=True,
        sourcePositionsFacesUVsUntouched=True, duplicateFaceRemovalForNormalComputationOnly=True,
        unchangedCompleteSourceStarsLiteral=True, anatomicalAcceptance=False,
        wallSeconds=wall_seconds, elapsedSeconds=time.monotonic() - started)
