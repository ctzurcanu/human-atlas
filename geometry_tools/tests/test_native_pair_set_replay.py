"""Native enumeration order must not weaken complete geometric pair replay."""
import copy,importlib.util,unittest
from pathlib import Path

spec=importlib.util.spec_from_file_location('radicular_pair_replay',Path(__file__).resolve().parents[2]/'scripts/verify-radicular-source-cleanup.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def diagnostic():
    return dict(geometry_sha256='literal',source_face_pairs_complete=True,self_intersection_pairs=2,
        source_face_pairs=[[1,3],[4,7]],source_face_pair_examples=[[1,3]],native_pair_replay_examples=[])

class NativePairSetReplayTests(unittest.TestCase):
    def test_different_enumeration_and_example_order_is_accepted(self):
        a=diagnostic();b=copy.deepcopy(a);b['source_face_pairs']=[[7,4],[3,1]];b['source_face_pair_examples']=[[4,7]]
        module.same_diagnostic(a,b)

    def test_changed_pair_is_rejected_even_with_the_same_count(self):
        a=diagnostic();b=copy.deepcopy(a);b['source_face_pairs']=[[1,3],[4,8]]
        with self.assertRaises(AssertionError):module.same_diagnostic(a,b)

    def test_incomplete_pair_evidence_is_rejected(self):
        a=diagnostic();a['source_face_pairs_complete']=False
        with self.assertRaises(AssertionError):module.same_diagnostic(a,copy.deepcopy(a))

if __name__=='__main__':unittest.main()
