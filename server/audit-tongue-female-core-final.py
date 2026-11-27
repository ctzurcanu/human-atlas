"""Independent bank read and finite geometry audit of the six-part female core."""
import argparse,itertools,json,sys,time
from pathlib import Path
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.constraints import check_constraint
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def read(p):return json.loads(Path(p).read_text())
def require(x,msg):
    if not x:raise ValueError(msg)
def main(ep,out):
    began=time.monotonic();e=read(ep);require(out.is_relative_to(W)and not out.exists(),'Use a new final audit stage');out.mkdir(parents=True);cp=Path(e['candidate']['path']);csha=digest_file(cp);require(csha==e['candidate']['sha256'],'Final emitted catalogue changed');atlas=read(cp);records={r['id']:r for r in atlas['parts']};ids=e.get('new_part_ids',e.get('active_part_ids'));require(len(ids)==6,'Audit requires exactly six active core/IL parts');require(all(any(q in pid for q in [':mucosa',':aponeurosis',':intrinsic:'])for pid in ids),'An extrinsic slipped into the active finite core');rows={r['part_id']:r for r in e['actual_world_meshes']};reader=AtlasReader(ROOT);meshes={};selfrows=[];trees={};solids={};origin=np.array([0.,1.5,0.]);bindings={str(ep):digest_file(ep),str(cp):csha,str(Path(__file__).resolve()):digest_file(__file__)}
    for i,pid in enumerate(ids):
        require(time.monotonic()-began<45,'Final core audit deadline approached');row=rows[pid];p=Path(row['mesh']['path']);require(digest_file(p)==row['mesh']['sha256'],'Emitted world packet changed');packet=Mesh.load(p);m=reader.read_part(atlas,records[pid],frame=packet.frame,catalogue_path=cp,catalogue_sha=csha);require(m.geometry_sha==row['geometry_sha256']==packet.geometry_sha and m.attributes_sha==packet.attributes_sha,'Independent actual bank read disagrees with decoded packet');require(np.array_equal(m.triangles,packet.triangles),'Actual face order differs from packet');meshes[pid]=m;trees[pid]=BVH(m);bindings[str(p)]=digest_file(p);bindings[m.source['bank_path']]=m.source['bank_sha256'];check=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000)
        signs=m.attributes.get('normalSigns',np.ones(len(m.vertices),dtype='i1'));words=m.attributes['normals'][m.faces];effective=words.astype(float)*signs[m.faces][:,:,None];g=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0]);dots=(effective*g[:,None,:]).sum(2);backward=int((dots<=0).sum());solid,_=make_solid(m.vertices-origin,m.faces);require(solid.status()==md.Error.NoError,'Actual closed core failed solid conversion');solids[pid]=solid
        item=dict(part_id=pid,mesh=str(p),mesh_sha256=digest_file(p),actual_geometry_sha256=m.geometry_sha,actual_attributes_sha256=m.attributes_sha,actual_bank_sha256=m.source['bank_sha256'],independent_bank_read_matches_decoded_packet=True,actual_self_topology=check,normalSigns_applied=True,nonpositive_effective_corner_dots=backward,minimum_effective_corner_dot=float(dots.min()),passed=check['passed']and backward==0);selfrows.append(item);atomic_json(out/('self-'+str(i)+'.json'),item)
    interfaces=[];shared={}
    def chart_faces(a,b):return [f for group in records[a].get('tissueDrawGroups',[])if group.get('neighbor')==b for f in range(group['start']//3,(group['start']+group['count'])//3)]
    for a,b in itertools.combinations(sorted(ids),2):
        af,bf=chart_faces(a,b),chart_faces(b,a)
        if not af and not bf:continue
        require(af and bf,'One-sided declared core tissue interface');spec=dict(kind='shared-interface',parameters=dict(first_faces=af,second_faces=bf));proof=check_constraint(spec,meshes[a],meshes[b]);shared[(a,b)]=(set(af),set(bf));interfaces.append(dict(first=a,second=b,first_faces=af,second_faces=bf,first_geometry_sha256=meshes[a].geometry_sha,second_geometry_sha256=meshes[b].geometry_sha,first_attributes_sha256=meshes[a].attributes_sha,second_attributes_sha256=meshes[b].attributes_sha,actual_check=proof))
    require(len(interfaces)==2 and sorted(len(p['first_faces'])for p in interfaces)==[62,122],'Final core lost its two literal reciprocal interface families')
    pairs=[]
    for i,(a,b)in enumerate(itertools.combinations(sorted(ids),2)):
        require(time.monotonic()-began<45,'Final pair audit deadline approached');arrays,summary=pair_contacts(meshes[a],meshes[b],0.,first_bvh=trees[a],second_bvh=trees[b],max_pairs=250000,max_tests=2000000);p=out/('pair-'+str(i).zfill(2)+'.npz');atomic_npz(p,**arrays);coplanar=arrays['pairs'][arrays['relation']=='coplanar-overlap'];declared=shared.get((a,b));undeclared=[q.tolist()for q in coplanar if declared is None or int(q[0])not in declared[0]or int(q[1])not in declared[1]];intersection=solids[a]^solids[b];require(intersection.status()==md.Error.NoError,'Actual core intersection computation failed');volume=float(intersection.volume());passed=summary['complete_within_pair_budget']and summary['crossing_face_pairs']==0 and not undeclared and volume==0.
        item=dict(first=a,second=b,first_geometry_sha256=meshes[a].geometry_sha,second_geometry_sha256=meshes[b].geometry_sha,first_attributes_sha256=meshes[a].attributes_sha,second_attributes_sha256=meshes[b].attributes_sha,summary=summary,witnesses=str(p),witnesses_sha256=digest_file(p),literal_shared_interface=(a,b)in shared,undeclared_coplanar_pairs=undeclared,actual_occupied_intersection_m3=volume,occupied_volume_exactly_zero=volume==0.,passed=passed);pairs.append(item);atomic_json(out/('progress-pair-'+str(i).zfill(2)+'.json'),item)
    passed=all(p['passed']for p in selfrows+pairs)and all(p['actual_check']['passed']and p['actual_check']['literal_opposed_normal_words']for p in interfaces);require(all(digest_file(p)==sha for p,sha in bindings.items()),'Bound actual packet/bank changed during audit');report=dict(schema='female-final-six-core-independent-audit/v1',packet=str(ep),packet_sha256=digest_file(ep),candidate=str(cp),candidate_sha256=csha,actual_float32_banks_independently_read=True,source_bindings=bindings,self_checks=selfrows,reciprocal_interfaces=interfaces,all_fifteen_pair_constraints=pairs,all_actual_jobs_complete=True,passed=passed,nonpositive_corner_normals=sum(p['nonpositive_effective_corner_dots']for p in selfrows),self_defect_pairs=sum(p['actual_self_topology']['self_intersections']['defect_pairs']for p in selfrows),proper_pair_crossings=sum(p['summary']['crossing_face_pairs']for p in pairs),anatomical_attachment_acceptance_inferred=False,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began);atomic_json(out/'evidence.json',report);print(json.dumps({k:v for k,v in report.items()if k not in ['source_bindings','self_checks','reciprocal_interfaces','all_fifteen_pair_constraints']}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--packet',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.packet.resolve(),a.output.resolve())
