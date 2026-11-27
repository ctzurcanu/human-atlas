import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.tissue import validate_cells


def tetra(offset=(0,0,0),size=1.):
    vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)*size+offset
    return Mesh(vertices,np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]]),'test:metres')


class TissueCheckpoints(unittest.TestCase):
    def test_only_the_same_actual_geometry_reuses_its_closed_cell_audit(self):
        with tempfile.TemporaryDirectory() as directory,patch('geometry_tools.tissue.audit',wraps=audit) as spy:
            first=validate_cells({'cell':tetra()},checkpoint_directory=directory)
            second=validate_cells({'cell':tetra()},checkpoint_directory=directory)
            self.assertTrue(first['passed']);self.assertTrue(second['passed']);self.assertEqual(spy.call_count,1)
            validate_cells({'cell':tetra(size=1.1)},checkpoint_directory=directory)
            self.assertEqual(spy.call_count,2)

    def test_an_edited_checkpoint_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            validate_cells({'cell':tetra()},checkpoint_directory=directory)
            file=next(Path(directory).glob('*.json'));record=json.loads(file.read_text());record['result']['passed']=False;file.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError,'Mismatched'):
                validate_cells({'cell':tetra()},checkpoint_directory=directory)

    def test_overlap_rejection_cannot_be_reused_for_a_different_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            overlapping={'a':tetra(),'b':tetra((.1,.1,.1))}
            first=validate_cells(overlapping,checkpoint_directory=directory)
            self.assertFalse(first['passed']);self.assertTrue(any(x.startswith('volume_overlap:') for x in first['failures']))
            self.assertFalse(validate_cells(overlapping,checkpoint_directory=directory)['passed'])
            self.assertTrue(validate_cells({'a':tetra(),'b':tetra((2.,0,0))},checkpoint_directory=directory)['passed'])

    def test_triangle_work_budget_is_explicit_and_part_of_the_cache_signature(self):
        with tempfile.TemporaryDirectory() as directory,patch('geometry_tools.tissue.audit',wraps=audit) as spy:
            short=validate_cells({'cell':tetra()},max_tests=1,checkpoint_directory=directory)
            self.assertFalse(short['passed'])
            self.assertFalse(short['cells']['cell']['self_intersections']['complete_within_pair_budget'])
            complete=validate_cells({'cell':tetra()},max_tests=100,checkpoint_directory=directory)
            self.assertTrue(complete['passed'])
            self.assertEqual(complete['cells']['cell']['self_intersections']['triangle_test_budget'],100)
            self.assertEqual(spy.call_count,2)
            self.assertTrue(validate_cells({'cell':tetra()},max_tests=100,checkpoint_directory=directory)['passed'])
            self.assertEqual(spy.call_count,2)


if __name__=='__main__':unittest.main()
