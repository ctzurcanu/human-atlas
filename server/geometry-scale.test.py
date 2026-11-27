import tempfile,unittest
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import AtlasReader

class GeometryScaleTest(unittest.TestCase):
 def test_literal_coordinates_scale_once_before_origin_and_reload(self):
  with tempfile.TemporaryDirectory()as d:
   root=Path(d);(root/'public').mkdir();v=np.array([[.001,.002,.003],[.01,.02,.03],[-.02,.01,.007]],dtype='<f4');n=np.tile(np.array([0,0,32767],'<i2'),(3,1));f=np.array([0,1,2],'<u4');raw=v.tobytes()+n.tobytes()+b'\0\0'+f.tobytes();(root/'public/mesh.bin').write_bytes(raw)
   p=dict(id='test',chunk=0,positions=0,normals=36,indices=56,vertexCount=3,indexCount=3,geometryOrigin=[0,1.2,0],geometryScale=.6258474041913706)
   a=dict(chunks=[dict(url='/mesh.bin')]);r=AtlasReader(root);m=r.read_part(a,p,frame='female');np.testing.assert_array_equal(m.vertices,v.astype(float)*p['geometryScale']+p['geometryOrigin']);np.testing.assert_array_equal(m.attributes['normals'],n);self.assertEqual((root/'public/mesh.bin').read_bytes(),raw)
   np.testing.assert_array_equal(r.reload(m.source,m.frame).vertices,m.vertices)
   for value in [0,-1,float('nan')]:
    with self.assertRaises(ValueError):r.read_part(a,{**p,'geometryScale':value},frame='female')

if __name__=='__main__':unittest.main()
