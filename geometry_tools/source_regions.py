"""Lossless indexed partitions of explicitly inventoried source surfaces.

This preserves an outer source surface; it does not construct finite walls,
name anatomical regions, or establish interfaces to another source mesh.
"""
import copy
import hashlib
import re
import numpy as np


class SourceRegionBank:
    """Copy source attributes once and assign every source face exactly once."""

    def __init__(self, chunk_index):
        if not isinstance(chunk_index, int) or chunk_index < 0:
            raise ValueError('An explicit destination chunk index is required')
        self.chunk_index = chunk_index
        self.data = bytearray()
        self.proofs = []

    def _append(self, raw, alignment=4):
        self.data.extend(b'\0' * (-len(self.data) % alignment))
        offset = len(self.data)
        self.data.extend(raw)
        return offset

    def partition(self, part, source_bank, regions, *, contract_sha256):
        """Return parts/leaf concepts from exhaustive, disjoint explicit faces.

        Region dictionaries require key, name, faces, qualification; term and
        concept_id are optional. Untyped remainders never acquire a whole-organ
        exact identifier. The caller supplies reviewed anatomical ownership.
        """
        if not re.fullmatch(r'[a-f0-9]{64}', contract_sha256):
            raise ValueError('A source-bound region contract hash is required')
        if part.get('geometryInstance') or part.get('colors') is not None:
            raise ValueError('Instanced or vertex-colour sources need a separately declared packing schema')
        nv, nf = part['vertexCount'], part['indexCount'] // 3
        if part['indexCount'] % 3 or nv <= 0 or nf <= 0 or len(regions) < 2:
            raise ValueError('A nonempty indexed source and at least two regions are required')
        faces = np.frombuffer(source_bank, '<u4', part['indexCount'], part['indices']).reshape(-1, 3)
        local = np.frombuffer(source_bank, '<f4', nv * 3, part['positions']).reshape(-1, 3)
        if faces.max() >= nv or not np.isfinite(local).all():
            raise ValueError('Invalid original source indices or coordinates')
        scale = part.get('geometryScale')
        scale = 1 if scale is None else scale
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError('A positive finite source scale is required')
        origin = np.asarray(part.get('geometryOrigin') if part.get('geometryOrigin') is not None else [0, 0, 0], float)
        offset = np.asarray(part.get('sourceOffset') if part.get('sourceOffset') is not None else [0, 0, 0], float)
        if origin.shape != (3,) or offset.shape != (3,) or not np.isfinite(np.r_[origin, offset]).all():
            raise ValueError('Finite source origin and offset vectors are required')
        world = local.astype(float) * scale + origin + offset
        owners = np.full(nf, -1, dtype=int)
        normalized = []
        keys = set()
        for i, region in enumerate(regions):
            ids = np.asarray(region['faces'])
            key = region['key']
            term = region.get('term')
            if not isinstance(key, str) or not key or key in keys or not region.get('name') or not region.get('qualification'):
                raise ValueError('Distinct region keys, names and qualifications are required')
            if term is not None and not re.fullmatch(r'A\d{2}\.\d\.\d{2}\.\d{3}[FM]?', term):
                raise ValueError('Invalid explicit anatomical term')
            if ids.ndim != 1 or ids.dtype.kind not in 'iu' or not len(ids) or ids.min() < 0 or ids.max() >= nf or len(np.unique(ids)) != len(ids):
                raise ValueError('Region faces must be nonempty, unique original integer indices')
            if np.any(owners[ids] != -1):
                raise ValueError('Two regions claim the same original face')
            owners[ids] = i
            normalized.append((region, ids.astype(np.uint32)))
            keys.add(key)
        if np.any(owners == -1):
            raise ValueError('Every original face must retain an owner')
        attrs, attribute_proof = {}, {}
        for field, width, alignment in [('positions', 12, 4), ('normals', 6, 2), ('normalSigns', 1, 1), ('uvs', 8, 4)]:
            if part.get(field) is None:
                continue
            raw = source_bank[part[field]:part[field] + nv * width]
            if len(raw) != nv * width:
                raise ValueError('Truncated source attribute bank')
            attrs[field] = self._append(raw, alignment)
            attribute_proof[field] = dict(sourceOffset=part[field], destinationOffset=attrs[field],
                bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
        parts, concepts, proof_regions = [], [], []
        for region, ids in normalized:
            term = region.get('term')
            pid = part['id'] + ':region:' + region['key']
            source_id = part['sourceId'] + '::SOURCE_REGION:' + region['key']
            index_bytes = faces[ids].astype('<u4').tobytes()
            index_offset = self._append(index_bytes)
            used = np.unique(faces[ids])
            p = copy.deepcopy(part)
            for key in ('suppressed', 'defaultHidden', 'replacementParts', 'tissueDrawGroups', 'ta98Term'):
                p.pop(key, None)
            p.update(id=pid, conceptId=region.get('concept_id', pid), sourceId=source_id,
                sourceOriginalId=part['sourceId'], sourceOriginalPartId=part['id'], name=region['name'],
                chunk=self.chunk_index, indices=index_offset, indexCount=len(ids)*3,
                bounds=[world[used].min(0).tolist(), world[used].max(0).tolist()], **attrs)
            if term:
                p['ta98Term'] = term
            p['sourceGroups'] = list(dict.fromkeys([*part.get('sourceGroups', []), part['conceptId']]))
            p['coverageLimitation'] = region['qualification']
            p['sourceDerivedRegion'] = dict(sourceId=part['sourceId'], sourcePartId=part['id'],
                sourceFaceIds=ids.tolist(), sourceFaceIdsSha256=hashlib.sha256(ids.astype('<u4').tobytes()).hexdigest(),
                boundaryContractSha256=contract_sha256, coordinateNormalUvMovement=False)
            p['provenance'] = dict(label='Original source, disjoint indexed surface region',
                url=part.get('provenance', {}).get('url', ''),
                detail='Original complete attributes and world pose retained; only explicitly owned original indices are drawn. ' + region['qualification'])
            c = dict(id=p['conceptId'], name=p['name'], elements=[pid], partOf=part['conceptId'],
                coverageLimitation=p['coverageLimitation'])
            if term:
                c.update(ta98Term=term, ta98Kind='exact')
            parts.append(p)
            concepts.append(c)
            proof_regions.append(dict(partId=pid, conceptId=c['id'], sourceId=source_id, term=term,
                sourceFaceIds=ids.tolist(), indexBytesSha256=hashlib.sha256(index_bytes).hexdigest()))
        self.proofs.append(dict(sourcePartId=part['id'], sourceId=part['sourceId'],
            sourceFaces=nf, originalAttributes=attribute_proof, regions=proof_regions,
            sourceFaceOwnershipExhaustive=True, sourceFacesDuplicated=0,
            newCoordinatesOrTriangles=0))
        return parts, concepts


def expand_source_ownership(atlas, ownership):
    """Expand existing source selectors and guarded manifests through a partition."""
    parts = {p['id']: p for p in atlas['parts']}
    def expand(ids):
        return [q for pid in ids for q in ownership.get(pid, [pid])]
    for c in atlas['concepts']:
        c['elements'] = expand(c['elements'])
    for b in atlas['ta98TermBindings']:
        b['partIds'] = expand(b['partIds'])
    for manifest in atlas.get('sourceCollectionIdentities', []):
        updated = []
        for member in manifest['members']:
            if member['partId'] not in ownership:
                updated.append(member)
                continue
            for pid in ownership[member['partId']]:
                p = parts[pid]
                updated.append(dict(partId=pid, sourceId=p['sourceId'], conceptId=p['conceptId'],
                    **({'ta98Term': p['ta98Term']} if p.get('ta98Term') else {}),
                    sourceOriginalId=p['sourceOriginalId'], sourceOriginalPartId=p['sourceOriginalPartId']))
        manifest['members'] = updated
