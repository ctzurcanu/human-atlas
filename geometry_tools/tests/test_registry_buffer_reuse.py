import unittest
from geometry_tools.registry_delta import reusable_geometry


class LiteralRegistryReuseTests(unittest.TestCase):
    def test_reuses_only_the_same_frame_complete_layout_path_and_bank_digest(self):
        record=dict(vertexCount=24,indexCount=60,positions=0,indices=288,geometryScale=.001,geometryOrigin=[0,0,0])
        prior=dict(error=None,frame_id='model:local-ta98:assembled-metres',bank_path='/local/geometry.bin',bank_sha256='checked-digest',record=record)
        self.assertTrue(reusable_geometry(prior,record,prior['frame_id'],prior['bank_path'],prior['bank_sha256']))
        self.assertFalse(reusable_geometry(prior,record,'model:registry-alias:assembled-metres',prior['bank_path'],prior['bank_sha256']))
        self.assertFalse(reusable_geometry(prior,record,prior['frame_id'],prior['bank_path'],'changed-digest'))
        self.assertFalse(reusable_geometry(prior,dict(record,geometryOrigin=[0,.001,0]),prior['frame_id'],prior['bank_path'],prior['bank_sha256']))
        self.assertFalse(reusable_geometry(prior,record,prior['frame_id'],'/other/geometry.bin',prior['bank_sha256']))


if __name__=='__main__':unittest.main()
