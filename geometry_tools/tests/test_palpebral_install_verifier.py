"""Independent rejection tests for the bounded ten-owner palpebral installer."""
import copy
import importlib.util
from pathlib import Path
import unittest

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.tests.test_coupled_source_displacement import pair, tube


spec = importlib.util.spec_from_file_location(
    'palpebral_install_verifier',
    Path(__file__).resolve().parents[2] / 'scripts/verify-palpebral-source-fit.py')
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)


def header():
    owners = verifier.EXPECTED_OWNERS
    primary = [pid for pid in owners if '_palpebral_' in pid]
    contract = dict(courseCase='explicit-palpebral-source-lid-fit', complete=True,
        sourceAtlas=str(Path('/test') / verifier.SOURCE_REVISION / 'circulation-atlas.json'),
        sourceAtlasSha256=verifier.SOURCE_SHA256, wholeAssembly=owners.copy(),
        primaryOwners=primary, recipientOwners=[p for p in owners if p not in primary],
        interfaces=[{} for _ in range(16)], externalInterfaces=[{} for _ in range(8)],
        pendingOffsetPairs=[{} for _ in range(6)], sectionScale=1.,
        solverEnergy='source-triangle-fem', **verifier.LIMITS)
    report = dict(complete=True, passedGeometryGates=True, rejected=[],
        sourceAtlasSha256=verifier.SOURCE_SHA256, rows=[dict(partId=p) for p in owners])
    return contract, report


def source_tube(*, alias=False):
    mesh = tube(alias=alias)
    return Mesh(mesh.vertices, mesh.faces, mesh.frame, mesh.attributes,
                dict(source_vertex_ids=list(range(len(mesh.vertices)))))


def moved(mesh, points, *, faces=None, attributes=None):
    return Mesh(points, mesh.faces if faces is None else faces, mesh.frame,
                mesh.attributes if attributes is None else attributes, mesh.source)


def source_row():
    return dict(fixedStarRawIds=[], sourceGauges=[list(range(i, i+6)) for i in range(0, 24, 6)])


class PalpebralInstallVerifierTests(unittest.TestCase):
    def test_exact_current_revision_complete_assembly_and_original_gates_are_required(self):
        contract, report = header()
        verifier.verify_header(contract, report)
        for field, replacement in [('sourceAtlasSha256', 'not-current'),
                                   ('wholeAssembly', contract['wholeAssembly'][:-1]),
                                   ('maximumSingularScale', 1.6),
                                   ('sectionScale', 1.01),
                                   ('solverEnergy', 'inverse-edge-length')]:
            candidate = copy.deepcopy(contract)
            candidate[field] = replacement
            with self.subTest(field=field), self.assertRaises(AssertionError):
                verifier.verify_header(candidate, report)

    def test_incomplete_missing_owner_or_rejected_candidate_cannot_pass_header(self):
        contract, report = header()
        for field, value in [('complete', False), ('passedGeometryGates', False),
                             ('rejected', ['normal gate']), ('rows', report['rows'][:-1]),
                             ('newIntrinsicPairs', 1), ('anatomicalAcceptance', True)]:
            candidate = copy.deepcopy(report)
            candidate[field] = value
            with self.subTest(field=field), self.assertRaises(AssertionError):
                verifier.verify_header(contract, candidate)

    def test_source_sections_accept_rigid_motion_without_appearance_or_face_changes(self):
        old = source_tube(alias=True)
        angle = .1
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0.],
                             [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
        new = moved(old, old.vertices @ rotation.T + [.0002, -.0001, .0003])
        proof = verifier.verify_mesh_geometry(old, new, source_row(), verifier.LIMITS)
        self.assertLess(proof['maximumSourceSectionPairDistanceError_m'], 1e-17)
        self.assertAlmostEqual(proof['metrics']['minimumSingularScale'], 1.)
        self.assertAlmostEqual(proof['metrics']['maximumGaugeAreaRatio'], 1.)
        topology = verifier.verify_intrinsic_topology(old, new)
        self.assertEqual(topology['addedSelfPairs'], [])

    def test_deleted_face_changed_attribute_or_split_alias_is_rejected(self):
        old = source_tube(alias=True)
        deleted = moved(old, old.vertices, faces=old.faces[:-1])
        with self.assertRaises(AssertionError):
            verifier.verify_mesh_geometry(old, deleted, source_row(), verifier.LIMITS)
        attributes = {key: value.copy() for key, value in old.attributes.items()}
        attributes['uv'][0, 0] += 1
        altered = moved(old, old.vertices, attributes=attributes)
        with self.assertRaises(AssertionError):
            verifier.verify_mesh_geometry(old, altered, source_row(), verifier.LIMITS)
        points = old.vertices.copy();points[-1, 0] += 1e-9
        with self.assertRaisesRegex(AssertionError, 'exact-alias partition'):
            verifier.verify_intrinsic_topology(old, moved(old, points))

    def test_fixed_stars_are_literal_not_within_a_large_packing_tolerance(self):
        old = source_tube()
        row = source_row();row['fixedStarRawIds'] = [0]
        new = moved(old, old.vertices + [1e-12, 0., 0.])
        with self.assertRaises(AssertionError):
            verifier.verify_mesh_geometry(old, new, row, verifier.LIMITS,
                                          agreement_m=1e-6, quantization_m=1e-6)

    def test_section_shape_and_motion_bounds_are_independently_measured(self):
        old = source_tube()
        for factor in [1.09, .94]:
            points = old.vertices.copy();points[:, :2] *= factor
            with self.subTest(factor=factor), self.assertRaises(AssertionError):
                verifier.verify_mesh_geometry(old, moved(old, points), source_row(), verifier.LIMITS)
        with self.assertRaises(AssertionError):
            verifier.verify_mesh_geometry(old, moved(old, old.vertices+[.013, 0., 0.]),
                                          source_row(), verifier.LIMITS)

    def test_complete_source_mouths_are_replayed_on_actual_new_points(self):
        originals, interfaces = pair(gap=2e-7)
        new = {pid: moved(mesh, mesh.vertices+[.0001, 0., 0.])
               for pid, mesh in originals.items()}
        proof = verifier.verify_known_mouths(originals, new, interfaces)
        self.assertTrue(proof[0]['before']['passed'] and proof[0]['after']['passed'])
        new['b'] = moved(originals['b'], originals['b'].vertices+[.0002, 0., 0.])
        with self.assertRaisesRegex(AssertionError, 'complete source mouth'):
            verifier.verify_known_mouths(originals, new, interfaces)
        incomplete = dict(interfaces[0], firstRing=list(range(18, 23)), secondRing=list(range(5)))
        with self.assertRaisesRegex(ValueError, 'complete external source mouth'):
            verifier.verify_known_mouths(originals, originals, [incomplete])

    def test_added_self_crossing_is_rejected_despite_unchanged_component_topology(self):
        points = np.array([[-1., -1., 0.], [1., -1., 0.], [0., 1., 0.],
                           [3., -.5, -1.], [3., .5, 1.], [3., 1., -1.]])
        old = Mesh(points, np.array([[0, 1, 2], [3, 4, 5]]), 'test:metres')
        new_points = points.copy();new_points[3:, 0] = 0.
        new = moved(old, new_points)
        with self.assertRaisesRegex(AssertionError, 'Added actual intrinsic crossing pairs'):
            verifier.verify_intrinsic_topology(old, new)

    def test_interface_identity_retains_correspondence_while_ignoring_common_loop_start(self):
        pair = dict(first='b', second='a', firstRing=[6, 7, 8], secondRing=[1, 2, 3])
        rotated = dict(first='a', second='b', firstRing=[2, 3, 1], secondRing=[7, 8, 6])
        self.assertEqual(verifier.interface_identity(pair), verifier.interface_identity(rotated))
        wrong = dict(rotated, secondRing=[8, 7, 6])
        self.assertNotEqual(verifier.interface_identity(pair), verifier.interface_identity(wrong))


if __name__ == '__main__':
    unittest.main()
