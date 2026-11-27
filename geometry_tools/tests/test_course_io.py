import json,tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.course_io import sha,load_course_spec,write_json


class CourseInputProtectionTest(unittest.TestCase):
    def fixture(self,folder):
        p=Path(folder);c=np.array([[0.,0.,0.],[0.,.01,0.],[0.,.02,0.]])
        np.savez(p/'graph.npz',centres=c,cells=np.arange(3),moving_cells=np.arange(3),fixed_cells=np.array([0,2]),edges=np.array([[0,1],[1,2]]),profile_cells=np.arange(3))
        np.savez(p/'reference.npz',centres=c);write_json(p/'catalogue.json',{})
        write_json(p/'source.json',dict(sourceProfileSha256=sha(p/'catalogue.json'),records=[dict(partId='tube',offset=0,vertices=2),dict(partId='protected-neighbour',offset=2,vertices=1)]))
        write_json(p/'course.json',dict(rows=[dict(index=i)for i in range(3)]))
        plan=dict(complete=True,completedRequestedSections=True,unsupportedSupportPlanes=[],coveredSourceIndices=[0,1,2],graphArraysSha256=sha(p/'graph.npz'),sourceGraphSha256=sha(p/'source.json'),courseReportSha256=sha(p/'course.json'),maximumMotion_m=.025,calibreCurvaturePlan=dict(limits_per_m={'1':80.}))
        write_json(p/'constraints.json',plan);spec=dict(schema='constrained-tube-course/v1',moving_parts=['tube'],maximum_motion_m=.025,fixed_axis=1)
        for key,name in dict(graph_arrays='graph.npz',source_graph='source.json',constraints='constraints.json',course='course.json',catalogue='catalogue.json',reference='reference.npz',source_curvature_arrays='graph.npz').items():spec[key]=str(p/name);spec[key+'_sha256']=sha(p/name)
        write_json(p/'spec.json',spec);return p,spec,plan

    def test_declared_fixed_and_unselected_neighbour_survive_explicit_moving_selection(self):
        with tempfile.TemporaryDirectory()as folder:
            p,spec,plan=self.fixture(folder);*_,active,fixed,edges,reference=load_course_spec(p/'spec.json')
            self.assertEqual(active,{0,1});self.assertEqual(fixed,{0,2})

    def test_stale_geometry_is_rejected_before_fitting(self):
        with tempfile.TemporaryDirectory()as folder:
            p,spec,plan=self.fixture(folder);(p/'graph.npz').write_bytes(b'changed geometry')
            with self.assertRaisesRegex(ValueError,'Stale source'):load_course_spec(p/'spec.json')

    def test_partial_station_coverage_is_not_a_complete_course(self):
        with tempfile.TemporaryDirectory()as folder:
            p,spec,plan=self.fixture(folder);plan['coveredSourceIndices']=[0,1];write_json(p/'constraints.json',plan);spec['constraints_sha256']=sha(p/'constraints.json');write_json(p/'spec.json',spec)
            with self.assertRaisesRegex(ValueError,'whole proposed course'):load_course_spec(p/'spec.json')


if __name__=='__main__':unittest.main()
