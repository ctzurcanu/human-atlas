"""Construct a bounded, source-bound finite tongue mucosa proposal.

Thickness is an explicit histology-informed modeling field, not a donor
measurement. The literal source outer surface is reused; no source is edited.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import time
import numpy as np
import igl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import Mesh
from geometry_tools.tissue import layer_arrangement, normal_words, pack

WORK = ROOT / '.local-models/geometry-workbench/tongue-completion'
REFERENCE = 'https://onlinelibrary.wiley.com/doi/10.1111/joa.13637'


def main(output: Path):
    began = time.monotonic()
    output = output.resolve()
    if output.exists() or not output.is_relative_to(WORK) or output == WORK:
        raise ValueError('Use a new stage in the local tongue workbench')
    source = WORK / 'source-v1/meshes/native-00.npz'
    native = Mesh.load(source)
    catalogue = Path(native.source['catalogue_path'])
    if digest_file(catalogue) != native.source['catalogue_sha256']:
        raise ValueError('Native source catalogue changed')
    output.mkdir(parents=True)
    # May 2022, table 1: central reticular tissue 300–350 µm and
    # epithelium 80–100 µm; inferior reticular tissue 140–190 µm and
    # epithelium 120–200 µm. Papillary projections are separate relief.
    # This broad field stays within those base-tissue ranges and does not
    # claim a measured local thickness or model the separate aponeurosis.
    normals = native.attributes['normals'].astype(float)
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    dorsal = np.clip(normals[:, 1], 0, 1)
    ventral = np.clip(-normals[:, 1], 0, 1)
    thickness = .00035 + .00007 * dorsal - .00005 * ventral
    inner = native.vertices - thickness[:, None] * normals
    plan = dict(schema='tongue-finite-mucosa-proposal/v1', term='A05.1.04.011',
        source=str(source), source_sha256=digest_file(source),
        source_geometry_sha256=native.geometry_sha,
        source_catalogue=str(catalogue), source_catalogue_sha256=digest_file(catalogue),
        source_bank_sha256=native.source['bank_sha256'],
        representation='finite-wall-layer', reference=REFERENCE,
        reference_scope='Base epithelial and reticular connective tissue ranges in human histology; papillary projections and aponeurosis are distinct.',
        thickness_scope='Histology-informed modeling field, not donor-specific measurement',
        minimum_thickness_m=float(thickness.min()), maximum_thickness_m=float(thickness.max()),
        excludes=['Lingual aponeurosis', 'Papillary projections', 'Muscle tissue'],
        adjoining_interfaces=[dict(first='Mucosa', second='Lingual aponeurosis or underlying connective tissue', status='Inner surface constructed; adjoining tissue pending')],
        source_geometry_modified=False, anatomical_acceptance=False, installed=False,
        wall_seconds=40)
    atomic_json(output / 'construction-input.json', plan)
    positions, surfaces, cells, prism = layer_arrangement(
        np.array([native.vertices, inner]), native.faces, ['mucosa'])
    # Retain the original artist normal words on the literal outer chart.
    words = normal_words(positions, surfaces)
    words['interface:0'] = native.attributes['normals'][native.faces].copy()
    result = pack(positions, surfaces, cells, {'mucosa':'TA98TONGUE:male:mucosa'},
        native.frame, output / 'packed', origin=(0., 1.5, 0.),
        corner_words=words, construction=prism, max_pairs=250000, max_tests=2000000)
    if time.monotonic()-began >= 40:
        raise TimeoutError('Mucosa construction exceeded its 40 second worker budget')
    # Verify the actual emitted inner chart is inside the literal source.
    encoded_inner = (inner - [0, 1.5, 0]).astype(np.float32).astype(float) + [0, 1.5, 0]
    queries = np.r_[encoded_inner, encoded_inner[native.faces].mean(axis=1)]
    winding = np.abs(igl.winding_number(native.vertices, native.faces, queries))
    encoded_error = float(np.linalg.norm(encoded_inner-inner, axis=1).max())
    passed = bool(result['checks']['passed'] and np.all(winding > .99))
    report = dict(**plan, packed_manifest=str(output / 'packed/tissue.json'),
        passed=passed, finite_prisms=prism,
        emitted_checks=result['checks'], stored_inner_contained=bool(np.all(winding>.99)),
        minimum_sampled_inner_winding=float(winding.min()),
        maximum_position_encoding_error_m=encoded_error,
        total_bytes=result['total_bytes'], new_image_assets=0,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output / 'construction.json', report)
    print({key:report[key] for key in ('passed','total_bytes','stored_inner_contained','elapsed_seconds')}, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output)
