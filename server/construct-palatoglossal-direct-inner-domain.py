"""One direct classifier from actual mucosal inner sheet and reflection rim.

The classifier is construction data only. Every inner/rim corner and raw word
is the literal opposite of the current finite M cell. No anatomy is installed.
"""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
W=ROOT/'.local-models/geometry-workbench/tongue-completion';O=np.array([0.,1.625,0.]);FRAME='model:local-ta98:assembled-metres'
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def words(m):return (m.attributes['normals'][m.faces].astype(np.int32)*m.attributes.get('normalSigns',np.ones(len(m.vertices),dtype='i1'))[m.faces][:,:,None]).astype('<i2')
def mesh_from_triangles(t,n):return Mesh(t.reshape(-1,3),np.arange(t.size//3).reshape(-1,3),FRAME,dict(normals=n.reshape(-1,3)))
def main(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use one new direct inner classifier stage')
    feasibility_path=W/'palatoglossal-inner-domain-feasibility-v1/evidence.json';feasibility=read(feasibility_path);d=feasibility['direct_inner_domain'];require(d['feasible_literal_closure'],'Direct source closure has not passed its exact edge census')
    for path,sha in feasibility['source_bindings'].items():require(digest_file(path)==sha,'Feasibility input changed: '+path)
    plan=read(W/'palatoglossal-original-parent-v2/plan.json');mu=Mesh.load(plan['source_files']['mucosa']);mt=mu.triangles;mw=words(mu)
    shared_path=W/'shared-coverings-exposed-apo-sl-v3/construction.json';construction=read(shared_path);ep=Path(construction['exposedSource']);fields=np.load(ep.with_name('fields.npz'));excluded=Mesh.load(ep.with_name('excluded-basal-carrier.npz'))
    encoded=(fields['fields'][0]-O).astype('<f4').astype(float)+O;inner=np.asarray(d['actual_inner_sheet_faces']);rim=np.asarray(d['actual_reflection_rim_faces']);literal=np.r_[inner,rim]
    anatomical=mt[literal][:,[0,2,1]];anatomical_words=-mw[literal][:,[0,2,1]];basal=encoded[excluded.faces]
    n=np.cross(basal[:,1]-basal[:,0],basal[:,2]-basal[:,0]);n/=np.linalg.norm(n,axis=1)[:,None];bn=np.repeat(np.rint(n*32767).astype('<i2')[:,None,:],3,axis=1)
    negative=mesh_from_triangles(np.r_[anatomical,basal],np.r_[anatomical_words,bn]);cap=mesh_from_triangles(basal,bn)
    output.mkdir(parents=True);negative.save(output/'inner-classifier.npz');cap.save(output/'basal-classifier-only.npz')
    actual=Mesh.load(output/'inner-classifier.npz');aa=audit(actual,require_closed=True,max_pairs=250000,max_tests=2000000);solid,_=make_solid(actual.vertices,actual.faces)
    require(np.array_equal(actual.triangles[:len(literal)],anatomical),'Classifier changed literal M inner/rim coordinates')
    require(np.array_equal(words(actual)[:len(literal)],anatomical_words),'Classifier changed literal opposed M inner/rim words')
    pair_rows=[]
    for name,host in [('arch_left',Mesh.load(plan['source_files']['arch_left'])),('arch_right',Mesh.load(plan['source_files']['arch_right'])),('mucosa',mu)]:
        require(time.monotonic()-began<45,'Direct classifier audit approached its finite deadline')
        witness,summary=pair_contacts(cap,host,0.,max_pairs=250000,max_tests=2000000);require(summary['complete_within_pair_budget'],'Classifier pair audit incomplete')
        path=output/('basal--'+name+'.npz');atomic_npz(path,**witness)
        literal_touch=[]
        for pair,relation in zip(witness['pairs'],witness['relation']):
            a,b=map(int,pair);common=set(map(tuple,cap.triangles[a]))&set(map(tuple,host.triangles[b]))
            literal_touch.append(dict(first=a,second=b,relation=str(relation),literal_shared_controls=len(common)))
        pair_rows.append(dict(host=name,summary=summary,witnesses=str(path),witnesses_sha256=digest_file(path),touch_ancestry=literal_touch,
            all_contacts_literal_source_boundary_only=all(r['relation']=='touching'and r['literal_shared_controls']>0 for r in literal_touch)))
    top=aa['topology'];positive=float(solid.volume())>0;strict_outside=all(not r['summary']['returned_triangle_pairs']for r in pair_rows)
    report=dict(schema='PG-direct-actual-inner-classifier/v1',source_feasibility=str(feasibility_path),source_feasibility_sha256=digest_file(feasibility_path),
        negative_mesh=str(output/'inner-classifier.npz'),negative_sha256=digest_file(output/'inner-classifier.npz'),negative_geometry_sha256=actual.geometry_sha,
        classifier_only_basal_faces=list(range(len(literal),len(actual.faces))),original_inner_source_faces=inner.tolist(),original_rim_source_faces=rim.tolist(),
        literal_opposed_inner_rim_coordinates=True,literal_opposed_inner_rim_raw_words=True,apo_chart_preserved_literal=True,
        exact_boundary_ancestry=d['feasible_literal_closure'],audit=aa,positive_volume_m3=float(solid.volume()),closed_components=len(top['components']),
        basal_pair_checks=pair_rows,basal_faces_strictly_outside_all_hosts=strict_outside,
        classifier_geometry_passed=aa['passed']and positive,strict_classifier_cut_gate_passed=aa['passed']and positive and strict_outside,
        source_geometry_modified=False,anatomical_meshes_added=0,qualifying_candidate=False,installed=False,anatomical_acceptance=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'evidence.json',report);print(json.dumps(dict(output=str(output),closed_components=len(top['components']),geometry_passed=report['classifier_geometry_passed'],volume_m3=float(solid.volume()),strict_cut_gate=report['strict_classifier_cut_gate_passed'],basal_pairs=[dict(host=r['host'],summary=r['summary'],only_literal_boundary=r['all_contacts_literal_source_boundary_only'])for r in pair_rows],elapsed_seconds=report['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
