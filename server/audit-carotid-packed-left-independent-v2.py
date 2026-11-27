"""Independent actual finite-bank audit. No catalogue, geometry or DB changes."""
from pathlib import Path
from collections import Counter
import hashlib,json,time,sys,itertools,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
from geometry_tools.algorithms.finite_cell_bank import PackedFiniteCells
from geometry_tools.tissue import decode_packed
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
import manifold3d as md
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1';O=B/'packed-left-v2';FRAME='model:local-ta98:assembled-metres'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,default=lambda v:v.tolist()if hasattr(v,'tolist')else v.item())+'\n')
def canon(f):return min(tuple(f),tuple(np.roll(f,1)),tuple(np.roll(f,2)))
def main():
 start=time.monotonic();rh=sha(O/'readiness.json');ready=json.loads((O/'readiness.json').read_text());data=(O/'cells.bin').read_bytes();bh=hashlib.sha256(data).hexdigest();assert bh==ready['bank_sha256'];planpath=B/'common-bank-prepacking-v1/plan.json';assert sha(planpath)==ready['plan_sha256'];plan=json.loads(planpath.read_text());source=ROOT/plan['source_catalogue'];assert sha(source)==ready['source_catalogue_sha256'];cat=json.loads(source.read_text());origin=np.array(ready['origin']);rows=ready['parts'];first=next(iter(rows.values()));n=first['vertexCount'];pos=np.ndarray((n,3),'<f4',buffer=data,offset=first['positions']);words=np.ndarray((n,3),'<i2',buffer=data,offset=first['normals']);uv=np.ndarray((n,2),'<f4',buffer=data,offset=first['uvs']);packed=PackedFiniteCells(data,pos,words,first['positions'],first['normals'],rows,{});meshes=decode_packed(packed,FRAME,origin);assert len(meshes)==10;partresult={};solids={};rawfaces={};rawsign={};allids=set(ready['part_ids'].values());fail=[]
 for role,m0 in meshes.items():
  row=rows[role];assert np.array_equal(row['geometryOrigin'],origin);m=Mesh(m0.vertices,m0.faces,FRAME,m0.attributes,dict(kind='actual-carotid-finite-bank',bank_path=str(O/'cells.bin'),bank_sha256=bh,readiness_sha256=rh,source_catalogue_sha256=ready['source_catalogue_sha256'],part_id=row['id']));ff=np.ndarray((row['indexCount']//3,3),'<u4',buffer=data,offset=row['indices']);sg=np.ndarray((n,),'i1',buffer=data,offset=row['normalSigns']);rawfaces[role]=ff;rawsign[role]=sg;norm=m.attributes['normals'].astype(float)*m.attributes['normalSigns'][:,None];tr=m.triangles;cross=np.cross(tr[:,1]-tr[:,0],tr[:,2]-tr[:,0]);bad=np.flatnonzero(np.any(np.einsum('ij,ikj->ik',cross,norm[m.faces])<=0,axis=1));a=audit(m,require_closed=True,max_pairs=600000,max_tests=2000000);volume=None
  if a['passed']:
   solid,_=make_solid((m.vertices-origin)*1000,m.faces);assert solid.status()==md.Error.NoError;volume=float(solid.volume());solids[role]=solid
  passed=a['passed']and not len(bad)and volume is not None and volume>0;path=O/('audited-'+role.replace(':','_')+'.npz');m.save(path);partresult[row['id']]=dict(passed=passed,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,role=role,decoded_path=str(path),actual_audit=a,nonforward_corner_faces=bad.tolist(),positive_tissue_volume_mm3=volume)
  if not passed:fail.append('cell:'+role)
  write(O/'independent-cell-stage.json',dict(complete=len(partresult)==10,readiness_sha256=rh,bank_sha256=bh,parts=partresult,failures=fail,elapsed_seconds=time.monotonic()-start));print('actual cell',role,'PASS'if passed else'FAIL',flush=True)
 assert not fail,fail
 interfaces=[];seen=set();byid={x['id']:k for k,x in rows.items()};kept=[];signs=np.zeros(n,'i1')
 for role,row in rows.items():
  for g in row['tissueDrawGroups']:
   faces=rawfaces[role][g['start']//3:(g['start']+g['count'])//3];neighbor=g.get('neighbor')
   if neighbor in allids:
    peer=byid[neighbor];key=tuple(sorted([role,peer]))
    if key in seen:continue
    seen.add(key);pg=[q for q in rows[peer]['tissueDrawGroups']if q.get('neighbor')==row['id']];assert len(pg)==1;q=pg[0];other=rawfaces[peer][q['start']//3:(q['start']+q['count'])//3];assert Counter(canon(f)for f in faces)==Counter(canon(f[[0,2,1]])for f in other);used=np.unique(faces);assert np.array_equal(rawsign[role][used],-rawsign[peer][used]);interfaces.append(dict(first=role,second=peer,triangles=len(faces),same_coordinate_normal_UV_indices=True,opposed_winding_and_signs=True))
   else:
    kept.append(faces);used=np.unique(faces);assert np.all((signs[used]==0)|(signs[used]==rawsign[role][used]));signs[used]=rawsign[role][used]
 assert len(interfaces)==9;signs[signs==0]=1;assembled=Mesh(pos.astype(float)+origin,np.concatenate(kept),FRAME,dict(normals=words,normalSigns=signs,uvs=uv),dict(kind='actual-carotid-and-lingual-assembled-wall',bank_sha256=bh));ass=audit(assembled,require_closed=True,max_pairs=1000000,max_tests=3000000);assert ass['passed'],ass['failures'];assert len(ass['topology']['components'])==1;genus=ass['topology']['components'][0]['closed_orientable_genus'];assembled.save(O/'audited-assembled-wall.npz')
 # Existing retained lingual courses retain ordered faces and shader fields; only authorised storage error is allowed.
 oldparts={x['id']:x for x in cat['parts']};reader=AtlasReader();retained=[]
 for role,m in meshes.items():
  if not role.startswith('retained-lingual:'):continue
  id=rows[role]['id'];old=reader.read_part(cat,oldparts[id],frame=FRAME,catalogue_path=source,catalogue_sha=ready['source_catalogue_sha256']);assert m.faces.shape==old.faces.shape;err=float(np.linalg.norm(m.vertices[m.faces]-old.vertices[old.faces],axis=2).max());assert err<=ready['maximum_retained_lingual_position_error_m']+1e-15;assert np.array_equal(m.attributes['normals'][m.faces],old.attributes['normals'][old.faces]);assert np.array_equal(m.attributes['uvs'][m.faces],old.attributes['uvs'][old.faces]);assert np.array_equal(m.attributes['normalSigns'][m.faces],old.attributes['normalSigns'][old.faces]);retained.append(dict(id=id,maximum_authorized_storage_error_m=err,normal_words_UV_signs_literal=True))
 assert sha(O/'readiness.json')==rh and sha(O/'cells.bin')==bh and sha(source)==ready['source_catalogue_sha256'];report=dict(schema='carotid-packed-catalogue-audit/v1',complete=True,passed=True,source_proof_passed=True,readiness_sha256=rh,bank_sha256=bh,source_catalogue_sha256=ready['source_catalogue_sha256'],parts=partresult,actual_bank_bindings=dict(path=str(O/'cells.bin'),bytes=len(data),origin_m=origin.tolist(),frame=FRAME,readiness_parts_exactly_decoded=True),shared_interfaces=interfaces,assembled_audit=ass,assembled_genus=genus,expected_blood_ports_if_genus7=8,assembled_genus_matches_expected=genus==7,retained_source_preservation=retained,source_plan_sha256=ready['plan_sha256'],clinical_limits='Existing IX/IJV overlaps and source open SCM volume limits remain unresolved; these are finite tissue/packing checks, not full anatomical acceptance.',anatomical_acceptance=False,occupied_pair_stage_pending=True,elapsed_seconds=time.monotonic()-start);write(O/'independent-catalogue-audit.json',report);print(json.dumps(dict(complete=True,passed=True,cells=len(partresult),interfaces=len(interfaces),assembled_genus=genus,elapsed=report['elapsed_seconds'])))
if __name__=='__main__':main()
