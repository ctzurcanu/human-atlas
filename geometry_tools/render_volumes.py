"""Audit closed cut-face classifiers without promoting them to tissue meshes."""
import numpy as np

from .audit import audit
from .common import digest_file
from .model import AtlasReader, Mesh


def validate_liver_sections(atlas, *, reader=None, frame=None):
    marker = atlas.get('liverSectionIntegration')
    if marker is None:
        return []
    liver = atlas.get('liverIntegration') or {}
    if (atlas.get('sex') not in ['male','female'] or marker.get('schema') != 'liver-section-volumes/v1'
            or marker.get('renderOnly') is not True
            or marker.get('bankUrl') != liver.get('bankUrl')
            or marker.get('bankSha256') != liver.get('bankSha256')):
        raise ValueError('Liver section classifiers require their explicit local assembly bank')
    volumes, exclusions = marker.get('volumes', []), marker.get('exclusions', [])
    expanded=liver.get('schema','').endswith('/v2')
    expected={'right','left','quadrate','caudate'}|({'fibrous-appendix','papillary-process','caudate-process'}if expanded else set())
    if (len(volumes) != len(expected) or {v.get('id') for v in volumes} != expected
            or len(exclusions) != 5 or {v.get('id') for v in exclusions} != {'artery','portal','biliary','hepatic_vein','retained-native-parents'}):
        raise ValueError('Incomplete liver bulk or conduit section classifiers')
    reader = reader or AtlasReader()
    path = reader.asset(marker['bankUrl'])
    sha = digest_file(path)
    if sha != marker['bankSha256']:
        raise ValueError('Liver section classifier bank changed')
    bank = np.memmap(path, dtype='u1', mode='r')
    parts = {p['id']: p for p in atlas['parts']}
    checks = []
    for v in volumes + exclusions:
        if v in volumes:
            pid = 'TA98PART:ZA:Liver:' + v['id']
            part = parts.get(pid)
            finite_subdivision=expanded and v['id']not in {'right','quadrate'}
            if (not part or v.get('partId') != pid or part.get('sectionSolid') is not finite_subdivision
                    or atlas['chunks'][part['chunk']]['url'] != marker['bankUrl']):
                raise ValueError('A section classifier cannot replace or close the visible source hilum')
        origin = np.asarray(v.get('origin'), float)
        if origin.shape != (3,) or not np.isfinite(origin).all():
            raise ValueError('Invalid section classifier storage origin')
        for key in ['positions','indices','vertexCount','indexCount']:
            if type(v.get(key)) is not int or v[key] < 0:
                raise ValueError('Invalid section classifier array descriptor')
        if v['positions'] % 4 or v['indices'] % 4 or v['indexCount'] % 3:
            raise ValueError('Misaligned section classifier arrays')
        scale=float(v.get('geometryScale',1.))
        if not np.isfinite(scale)or scale<=0:raise ValueError('Invalid classifier scale')
        vertices = reader.view(bank, v['positions'], v['vertexCount'], '<f4', 3).astype(float)*scale + origin
        faces = reader.view(bank, v['indices'], v['indexCount'], '<u4').reshape(-1, 3)
        mesh = Mesh(vertices, faces, frame or 'model:local-ta98:assembled-metres')
        if mesh.geometry_sha != v.get('geometrySha256'):
            raise ValueError('Section classifier geometry identity changed')
        result = audit(mesh, require_closed=True, max_pairs=2000000, max_tests=16000000)
        if not result['passed']:
            raise ValueError('Invalid emitted closed section classifier: ' + v['id'])
        checks.append(dict(id=v['id'], representation='render-only-section-volume',
                           bank_path=str(path), bank_sha256=sha, audit=result, passed=True,
                           new_tissue_solid=False, added_ta98_terms=0, anatomical_acceptance=False))
    return checks
