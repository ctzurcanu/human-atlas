import numpy as np
import shapely
import unittest
from geometry_tools.algorithms.section_pair_packing import minimum_cost_clear_pair


def square(x,z,size=.0002):
 return shapely.polygons([[x-size,z-size],[x+size,z-size],[x+size,z+size],[x-size,z+size]])


def test_nearest_shortlist_must_not_establish_infeasibility():
 first=np.array([square(0,0)]*80+[square(.002,0)],object)
 second=np.array([square(0,0)]*80+[square(-.003,0)],object)
 cost=np.r_[np.full(80,.001),.002]
 result=minimum_cost_clear_pair(first,second,cost,np.r_[np.full(80,.001),.003],wall_seconds=2)
 assert result['complete']and result['feasible']
 assert result['indices']==[80,0]
 assert result['minimumMaximumCost']==.002


def test_finite_optimum_matches_full_pair_matrix():
 rng=np.random.default_rng(17);a=rng.uniform(-.003,.003,(20,2));b=rng.uniform(-.003,.003,(25,2))
 first=shapely.polygons(a[:,None]+np.array([[-.0002,-.0002],[.0002,-.0002],[.0002,.0002],[-.0002,.0002]]))
 second=shapely.polygons(b[:,None]+np.array([[-.0002,-.0002],[.0002,-.0002],[.0002,.0002],[-.0002,.0002]]))
 ac=np.linalg.norm(a,axis=1);bc=np.linalg.norm(b,axis=1)
 valid=shapely.disjoint(first[:,None],second[None])&(shapely.distance(first[:,None],second[None])>=.0001)
 expected=np.where(valid,np.maximum(ac[:,None],bc[None]),np.inf).min()
 result=minimum_cost_clear_pair(first,second,ac,bc,wall_seconds=2)
 assert result['complete']and np.isclose(result['minimumMaximumCost'],expected,rtol=0,atol=1e-15)


def test_infeasible_finite_set_and_expired_deadline_are_distinct():
 a=np.array([square(0,0)],object)
 result=minimum_cost_clear_pair(a,a,[1.],[1.],wall_seconds=2)
 assert result['complete']and not result['feasible']
 result=minimum_cost_clear_pair(a,a,[1.],[1.],wall_seconds=1e-15)
 assert not result['complete']


def load_tests(loader,tests,pattern):
 return unittest.TestSuite(unittest.FunctionTestCase(fn)for name,fn in globals().items()if name.startswith('test_'))
