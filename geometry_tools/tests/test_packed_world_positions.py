import unittest
import numpy as np
from geometry_tools.algorithms.packed_world_positions import encode_world_positions,decode_source_positions


class PackedWorldPositionTests(unittest.TestCase):
    def test_nonzero_hand_origin_is_applied_once_and_source_words_replay(self):
        part=dict(geometryOrigin=[.3,.9,.04])
        raw=np.array([[.01,-.02,.03],[-.04,-.1,.02]],dtype='<f4')
        world=decode_source_positions(raw,part)
        np.testing.assert_array_equal(encode_world_positions(world,part),raw)
        moved=world+[[0.,.001,0.],[.001,0.,0.]]
        replay=decode_source_positions(encode_world_positions(moved,part),part)
        self.assertLess(float(np.linalg.norm(replay-moved,axis=1).max()),1e-8)
        self.assertGreater(float(np.linalg.norm(decode_source_positions(moved.astype('<f4'),part)-moved,axis=1).min()),.9)

    def test_scale_offset_and_instance_are_inverted_in_the_declared_order(self):
        part=dict(geometryOrigin=[.25,.75,.04],geometryScale=1.25,sourceOffset=[.003,-.002,.001],geometryInstance=dict(rotation=[0.,0.,2**-.5,2**-.5],scale=.8))
        raw=np.array([[.01,-.03,.02],[-.02,.04,.05]],dtype='<f4')
        world=decode_source_positions(raw,part)
        np.testing.assert_array_equal(encode_world_positions(world,part),raw)
        for invalid in [dict(geometryScale=0),dict(geometryOrigin=[float('nan'),0,0])]:
            with self.assertRaises(ValueError):encode_world_positions(world,invalid)


if __name__=='__main__':unittest.main()
