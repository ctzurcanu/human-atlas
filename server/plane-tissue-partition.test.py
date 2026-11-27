"""Physical regression: a two-coat tube retains its lumen through three cuts."""
import unittest
import numpy as np
from geometry_tools.tissue import layer_arrangement
from geometry_tools.algorithms.plane_tissue_partition import partition
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

class PlanePartition(unittest.TestCase):
    def test_shared_caps_preserve_two_wall_layers_and_open_lumen(self):
        n=20;theta=np.arange(n)*2*np.pi/n;faces=[]
        for j in range(n):
            k=(j+1)%n;faces.extend([[j,k,n+k],[j,n+k,n+j]])
        fields=[]
        for radius in [.3,.25,.2]:
            ring=np.c_[radius*np.cos(theta),np.zeros(n),radius*np.sin(theta)]
            fields.append(np.r_[ring,ring+[0,1,0]])
        fields=np.asarray(fields);faces=np.asarray(faces)
        q=fields[0][faces];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
        if np.sum(normal[:,[0,2]]*q.mean(1)[:,[0,2]])<0:faces=faces[:,[0,2,1]]
        v,s,c,_=layer_arrangement(fields,faces,['outer','inner'])
        bands=[('a',-.1,.237123),('b',.237123,.511234),('c',.511234,.842345),('d',.842345,1.1)]
        result=partition(v,s,c,bands,origin=[0,.5,0])
        vv=result['positions'];volume=0.
        for name,uses in result['cells'].items():
            ff=np.concatenate([result['surfaces'][u['surface']]if u['side']==1 else result['surfaces'][u['surface']][:,[0,2,1]]for u in uses])
            t=topology(Mesh(vv,ff,'fixture:metres'))
            self.assertTrue(t['closed_oriented_manifold'],name);volume+=t['oriented_triangle_integral_m3']
        expected=n*np.sin(2*np.pi/n)/2*(.3**2-.2**2)
        self.assertAlmostEqual(volume,expected,places=8)
        for sid,ff in result['surfaces'].items():
            if not sid.startswith('cut:'):continue
            # No triangle may reach across the central lumen.
            centers=vv[ff].mean(1);self.assertGreater(np.linalg.norm(centers[:,[0,2]],axis=1).min(),.19)
            owners=[u for uses in result['cells'].values()for u in uses if u['surface']==sid]
            self.assertEqual(sorted(u['side']for u in owners),[-1,1])

if __name__=='__main__':unittest.main()
