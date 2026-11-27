"""Source-bound tongue surfaces, margin curves and apex; no duplicated tissue."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import igl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader
from tongue_whole_input import whole_input


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='surface-controls-v1')
    parser.add_argument('--whole-tongue',type=Path)
    parser.add_argument('--fit-result',type=Path)
    parser.add_argument('--carrier-part-id',default='ZA:Tongue')
    args = parser.parse_args()
    out = ROOT / '.local-models/geometry-workbench/tongue-completion' / args.output
    if out.exists():
        raise ValueError('Use a new output stage; reviewed evidence is immutable')
    inventory_path = ROOT / 'reports/ta98-tongue-inventory.json'
    inventory = json.loads(inventory_path.read_text())
    source = Path(inventory['source_catalogues']['male']['path'])
    assert digest_file(source) == inventory['source_catalogues']['male']['sha256']
    mesh,mesh_path,binding=whole_input(args.whole_tongue,args.fit_result)
    corrected=binding['mode']!='native-source-default'
    if not corrected and args.carrier_part_id!='ZA:Tongue':
        raise ValueError('A nonnative carrier needs an explicitly reviewed corrected whole')
    v, f, tri = mesh.vertices, mesh.faces, mesh.triangles
    centers = tri.mean(1)
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norms = np.linalg.norm(normals, axis=1)
    assert np.all(norms > 0)
    normals /= norms[:, None]
    bank = mesh.source['bank_sha256']
    host = args.carrier_part_id
    groups = []

    def anchor(face, bary=None):
        bary = np.asarray([1/3] * 3 if bary is None else bary, float)
        point = bary @ tri[face]
        assert abs(bary.sum() - 1) < 1e-14 and np.all(bary >= 0)
        return dict(partId=host, face=int(face), barycentric=bary.tolist(), point=point.tolist())

    def group(key, term, name, representation):
        return dict(key=key, term=term, name=name, conceptId=f'TA98TONGUE:male:feature:{key}',
                    kind='exact', parent='tongue', partIds=[host], representation=representation,
                    boundaryProvisional=True, anatomicalAcceptance=False)

    # Dorsum follows actual superior and upper posterior source facets. The
    # posterior root is not equivalent to the complete rear-facing surface.
    dorsal = (normals[:, 1] > .25) | ((normals[:, 2] < -.5) & (centers[:, 1] > np.median(v[:, 1])))
    # The terminal V is a reference-estimated extent, not a newly sculpted
    # groove, and is deliberately not counted as the terminal-sulcus term.
    z0, z1 = v[:, 2].min(), v[:, 2].max()
    terminal = z1 - (z1-z0) * 2/3 + .4 * np.abs(centers[:, 0])
    anterior = centers[:, 2] >= terminal
    masks = [
        ('dorsum', 'A05.1.04.004', 'Dorsum of tongue', dorsal),
        ('anterior-dorsum', 'A05.1.04.005', 'Anterior part of dorsum', dorsal & anterior),
        ('posterior-dorsum', 'A05.1.04.006', 'Posterior part of dorsum', dorsal & ~anterior),
        ('inferior-surface', 'A05.1.04.007', 'Inferior surface of tongue', (normals[:, 1] < -.25) & anterior),
    ]
    for key, term, name, mask in masks:
        faces = np.flatnonzero(mask)
        assert len(faces)
        g = group(key, term, name, 'surface-region')
        g['surface'] = dict(partId=host, faces=faces.tolist(), sourceBankSha256=bank)
        centre = centers[faces].mean(0)
        g['anchor'] = anchor(faces[np.argmin(np.linalg.norm(centers[faces] - centre, axis=1))])
        g['coverageLimitation'] = 'Actual stored tongue facets; regional boundaries are reference estimates, not donor-specific sulcal segmentation. No additional tissue is created.'
        groups.append(g)
    assert set(groups[0]['surface']['faces']) == set(groups[1]['surface']['faces']) | set(groups[2]['surface']['faces'])
    assert not set(groups[1]['surface']['faces']) & set(groups[2]['surface']['faces'])

    # The apex is a real source vertex, with an exact source-face attachment.
    tip = int(np.argmax(v[:, 2]))
    face, corner = np.argwhere(f == tip)[0]
    g = group('apex', 'A05.1.04.010', 'Apex of tongue', 'point-landmark')
    g['anchor'] = anchor(face, np.eye(3)[corner])
    g['boundaryProvisional'] = False
    g['coverageLimitation'] = 'Anterior-most vertex of the literal native tongue; a point reference, without an invented solid.'
    groups.append(g)

    # At each AP station take the actual lateral source-section extremum.
    # Keep each sample on its carrier; these curves have no solid occupancy.
    edges = np.unique(np.sort(np.r_[f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]], axis=1), axis=0)
    curves = []
    for side in (-1, 1):
        points = []
        for z in np.linspace(z1-(z1-z0)*.64, z1-1e-7, 45):
            e = v[edges]
            d = e[:, :, 2] - z
            crossing = d[:, 0] * d[:, 1] < 0
            e, d = e[crossing], d[crossing]
            if not len(e):
                continue
            t = d[:, 0] / (d[:, 0] - d[:, 1])
            values = e[:, 0] + t[:, None] * (e[:, 1] - e[:, 0])
            points.append(values[np.argmax(side * values[:, 0])])
        points.append(v[tip])
        points = np.asarray(points)
        distances, _, projected = igl.point_mesh_squared_distance(points, v, f)
        assert distances.max() < 1e-22
        curves.append(dict(partId=host, points=projected.tolist(), closed=False, sourceBankSha256=bank))
    g = group('margin', 'A05.1.04.009', 'Margin of tongue', 'boundary-curve')
    g['curves'] = curves
    g['anchor'] = groups[-1]['anchor']
    g['coverageLimitation'] = 'Reference-estimated lateral margins sampled on actual source triangles. No tube or extra anatomical solid is created.'
    groups.append(g)
    out.mkdir(parents=True)
    if corrected:
        # A construction NPZ is not the original runtime bank. The coupled
        # packer must verify the literal mucosa outer triangles, then bind
        # these references to that actual bank without inventing a carrier.
        for g in groups:
            if g.get('surface'):
                g['surface'].pop('sourceBankSha256')
                g['surface']['sourceMeshSha256']=digest_file(mesh_path)
            for curve in g.get('curves',[]):
                curve.pop('sourceBankSha256')
                curve['sourceMeshSha256']=digest_file(mesh_path)
    atomic_json(out / 'surface-controls.json', dict(schema='tongue-construction-features/v1' if corrected else 'tongue-source-features/v1',
        sourceCatalogueSha256=digest_file(source), inventorySha256=digest_file(inventory_path),
        hostPartId=host, carrierBankHashes={} if corrected else {host:bank}, termGroups=groups,
        wholeInput=binding,carrierMesh=str(mesh_path),carrierMeshSha256=digest_file(mesh_path),
        carrierBankBindingPending=corrected,
        sourceGeometryChanged=False, additionalSolidMeshes=0, additionalImageAssets=0,
        sourceFacets=len(f), terms=len(groups), anatomicalAcceptance=False,
        references=['https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/'],
        elapsedSeconds=time.monotonic()-started))
    print(json.dumps(dict(output=str(out), terms=len(groups), additionalMeshes=0,
                         elapsedSeconds=time.monotonic()-started)))


if __name__ == '__main__':
    main()
