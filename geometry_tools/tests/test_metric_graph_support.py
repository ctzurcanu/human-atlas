"""Keep metric source support independent of tessellation and spatial proximity."""
import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.metric_graph_support import metric_support


class MetricSupport(unittest.TestCase):
    def test_irregular_source_edge_lengths_replace_arbitrary_hop_count(self):
        c=np.c_[[0,.001,.002,.03],np.zeros(4),np.zeros(4)]
        graph=sparse.coo_matrix((np.ones(6),([0,1,1,2,2,3],[1,0,2,1,3,2])),shape=(4,4)).tocsr()
        ids,_=metric_support(graph,c,[0],.005);np.testing.assert_array_equal(ids,[0,1,2])

    def test_spatially_close_disconnected_source_is_not_added(self):
        c=np.array([[0.,0,0],[.001,0,0],[.0001,0,0]])
        graph=sparse.coo_matrix((np.ones(2),([0,1],[1,0])),shape=(3,3)).tocsr()
        ids,_=metric_support(graph,c,[0],.005);np.testing.assert_array_equal(ids,[0,1])

    def test_coincident_dependency_groups_are_not_welded_or_dropped(self):
        c=np.zeros((2,3));graph=sparse.coo_matrix((np.ones(2),([0,1],[1,0])),shape=(2,2)).tocsr()
        ids,_=metric_support(graph,c,[0],1e-5);np.testing.assert_array_equal(ids,[0,1])

if __name__=='__main__':unittest.main()
