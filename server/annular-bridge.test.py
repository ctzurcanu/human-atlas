from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.algorithms.annular_bridge import attach_bridge
from geometry_tools.algorithms.tube_continuity import loft_faces
from geometry_tools.algorithms.finite_vessel_wall import finite_wall
from geometry_tools.model import Mesh
from geometry_tools.audit import audit

def cylinder(n,z):
 a=np.arange(n)*2*np.pi/n;outer=np.array([np.c_[np.cos(a)*.002,np.sin(a)*.002,np.full(n,x)]for x in z]).reshape(-1,3);inner=outer.copy();inner[:,:2]*=.8
 wall,_,caps=finite_wall(outer,inner,loft_faces(2,n));return wall['PairedPositions'].astype(float),wall['FiniteFaces'].astype(int),caps

class AnnularBridgeTest(unittest.TestCase):
 def test_unequal_real_rims_preserve_lumen_and_controls(self):
  pv,pf,pc=cylinder(5,[0,.01]);cv,cf,cc=cylinder(12,[.02,.03])
  def ids(f,cap):
   keys={tuple(sorted(map(int,t)))for t in cap};return np.array([i for i,t in enumerate(f)if tuple(sorted(map(int,t)))in keys])
  parent_cap=max(pc,key=lambda c:pv[np.unique(c),2].mean());child_cap=min(cc,key=lambda c:cv[np.unique(c),2].mean())
  q=attach_bridge(pv,pf,ids(pf,parent_cap),cv,cf,ids(cf,child_cap))
  self.assertTrue(q['parent_source_controls_literal']);self.assertTrue(q['child_source_controls_literal']);self.assertEqual(q['lumen_disks_added'],0)
  m=Mesh(q['vertices'],q['faces'],'fixture:metres');r=audit(m,require_closed=True,max_pairs=100000,max_tests=1000000)
  self.assertTrue(r['passed'],r['failures']);self.assertEqual(len(r['topology']['components']),1)
  # The cut annulus retains its hole: no emitted face intersects the axis.
  tri=m.triangles[q['interface_faces']];self.assertGreater(np.linalg.norm(tri[:,:,:2],axis=2).min(),.0015)
if __name__=='__main__':unittest.main()
