"""Guarded genuine branch exterior constructor. No work runs on import.

Usage: --course-spec FILE --approval FILE --output NEW_DIRECTORY
The exact reviewed course spec and its source arrays must be approved first.
"""
import argparse,copy,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit,winding_contains,topology
from geometry_tools.algorithms.native_tree_exterior import construct
from geometry_tools.algorithms.tube_continuity import boundary_loops

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME='model:local-ta98:assembled-metres';PARENT='ZA:Thoraco-acromial artery.l'
CHILD='TA98PART:ZA:Thoraco-acromial artery.l:A12.2.09.008';TERM='A12.2.09.008'
H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(path,obj):path.write_text(json.dumps(obj,indent=2,default=lambda x:x.tolist()if hasattr(x,'tolist')else str(x))+'\n')
def generated_mesh(z,prefix):
 v=np.asarray(z[prefix+'_vertices'],float);f=np.asarray(z[prefix+'_faces'],int);fields=np.asarray(z[prefix+'_fields'],float);centres=np.asarray(z[prefix+'_centres'],float)
 assert np.array_equal(v,fields.reshape(-1,3)),'Explicit source fields must index literal course vertices'
 radial=fields-centres[:,None];radial/=np.linalg.norm(radial,axis=2)[...,None];n=np.rint(radial.reshape(-1,3)*32767).astype(np.int16)
 # The independently generated ring chart may use either traversal. Orient
 # its triangles against the explicit outward radial field before closing
 # diagnostic cells. This never changes native parent faces or any position.
 face_normal=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])
 radial_face=radial.reshape(-1,3)[f].mean(1)
 alignment=np.einsum('ij,ij->i',face_normal,radial_face)
 reversed_winding=bool(alignment.sum()<0)
 if reversed_winding:f=f[:,::-1].copy();alignment=-alignment
 assert np.all(alignment>0),'Generated tube has a locally inverted radial chart'
 distance=np.r_[0,np.cumsum(np.linalg.norm(np.diff(centres,axis=0),axis=1))];u=np.arange(fields.shape[1])/fields.shape[1];uv=np.stack(np.meshgrid(u,distance,indexing='xy'),axis=-1).reshape(-1,2).astype(np.float32)
 return Mesh(v,f,FRAME,dict(normals=n,uvs=uv),dict(kind='independently-authored-illustrative-course',vendor_mesh_or_UV_copied=False,outward_winding_corrected=reversed_winding)),np.arange(fields.shape[1])
def run(spec_path,approval_path,out,side='left'):
 global PARENT,CHILD
 if side=='right':PARENT='ZA:Thoraco-acromial artery.r';CHILD='TA98PART:ZA:Thoraco-acromial artery.r:A12.2.09.008'
 start=time.monotonic();out=out.resolve();spec=json.loads(spec_path.read_text());approval=json.loads(approval_path.read_text())
 assert approval['course_spec_sha256']==H(spec_path) and approval['course_passed']is True and approval['root_surface_course_review_passed']is True,'Root course approval required before any union'
 if side=='right':
  assert approval.get('right_context_review_passed')is True,'Right native recipients, veins/nerves and fixed context require their own review'
  context=Path(approval['right_context_evidence']);assert H(context)==approval['right_context_evidence_sha256']
 assert H(CAT)==spec['source_catalogue_sha256'];archive=Path(spec['course_archive']);assert H(archive)==spec['course_archive_sha256']
 out.mkdir(exist_ok=False);inv=Path(spec['authority_inventory']);assert H(inv)==spec['authority_inventory_sha256']
 inventory=json.loads(inv.read_text());assert len(inventory['terms'])==7
 (out/'authority-inventory.json').write_bytes(inv.read_bytes());write(out/'assembly-plan-before-union.json',dict(schema='native-parent-independent-clavicular-ostium/v1',source_catalogue_sha256=H(CAT),course_spec=spec,course_spec_sha256=H(spec_path),approval_sha256=H(approval_path),parts=[PARENT,CHILD],face_ownership='Parent role0; new main role1 and muscular twig role2 are one child part. Every output face owned exactly once.',interfaces=['Retain native parent proximal/distal profiles exactly.','Open both main and muscular tips as original continuations.','Internal child roots closed only with explicit completely buried conical anchors; reject any emitted continuation cone chart.','Parent-child shared interface is an open arterial ostium, shared geometric boundary and locally averaged packed normal words; no lumen disk.'],sourceOuterSurfaceOnly=True,finite_wall=False,perfusion_or_histology_complete=False,all_outside_anatomy_fixed=True,illustrative_variant=True,installed=False))
 atlas=json.loads(CAT.read_text());parent_record=next(p for p in atlas['parts']if p['id']==PARENT);parent=AtlasReader().read_part(atlas,parent_record,frame=FRAME,catalogue_path=CAT,catalogue_sha=H(CAT));z=np.load(archive)
 if side=='right':
  # Reflect only the finished independent course. The native right parent stays
  # literal; an exact geometric face correspondence maps the local ostium patch.
  z={k:z[k].copy()for k in z.files};mirror=np.array([-1.,1.,1.])
  for prefix in ['main','muscular']:
   for key in ['vertices','fields','centres']:z[prefix+'_'+key]*=mirror
   z[prefix+'_faces']=z[prefix+'_faces'][:,[0,2,1]]
  frozen_path=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracoacromial-clavicular-branch-v1/part_32.npz';frozen=Mesh.load(frozen_path,frame=FRAME)
  def signature(t):return frozenset(map(tuple,t))
  face_lookup={signature(t):i for i,t in enumerate(parent.triangles)};mapped=[]
  for index in spec['allowed_parent_face_ids']:
   target=signature(frozen.triangles[index]*mirror);assert target in face_lookup,'Right parent has no exact reflected source chart';mapped.append(face_lookup[target])
  spec=copy.deepcopy(spec);spec['allowed_parent_face_ids']=mapped;spec['conical_anchors_m']={k:(np.asarray(v)*mirror).tolist()for k,v in spec['conical_anchors_m'].items()};spec['blood_probe_points_m']=(np.asarray(spec['blood_probe_points_m'])*mirror).tolist();spec['construction_side']='right';spec['reflection_plane']='X=0; source course reflection only'
  plan=json.loads((out/'assembly-plan-before-union.json').read_text());plan['course_spec']=spec;plan['reflection_transfer']='Literal X=0 course reflection onto independently checked native right context';write(out/'assembly-plan-before-union.json',plan)
 main,mroot=generated_mesh(z,'main');twig,troot=generated_mesh(z,'muscular');sources=[parent,main,twig];nv=np.cumsum([0]+[len(m.vertices)for m in sources]);nf=np.cumsum([0]+[len(m.faces)for m in sources]);composite=Mesh(np.concatenate([m.vertices for m in sources]),np.concatenate([m.faces+nv[i]for i,m in enumerate(sources)]),FRAME,{k:np.concatenate([m.attributes[k]for m in sources])for k in ['normals','uvs']},dict(source_catalogue_sha256=H(CAT),independent_course_archive_sha256=H(archive)))
 mouths=[];anchors={}
 for role,m in enumerate(sources):
  loops=boundary_loops(m.faces);assert len(loops)==2
  for loop in loops:
   internal=role>0 and set(map(int,loop))==set(map(int,mroot if role==1 else troot));ids=loop+nv[role]
   mouths.append(dict(positions_m=composite.vertices[ids],source_control_ids=ids,kind='internal root'if internal else 'original'if role==0 else 'original continuation'))
   if internal:anchors[len(mouths)-1]=np.asarray(spec['conical_anchors_m']['main'if role==1 else 'muscular'],float)
 assert len(anchors)==2
 composite.save(out/'composite-input.npz');write(out/'explicit-mouths-and-anchors-before-union.json',dict(mouths=mouths,anchors=anchors,component_face_ranges=nf))
 arrays,proof=construct(composite,[np.arange(nf[i],nf[i+1])for i in range(3)],mouths,anchors,wall_seconds=35,max_pairs=600000,max_tests=2000000)
 np.savez_compressed(out/'union-source-provenance.npz',**arrays);write(out/'native-tree-proof.json',proof);assert proof['passed'],proof['failures'];assert proof['derived_fragments']==0,'Every conical root closure must be fully buried'
 changed=np.asarray(arrays['changed_source_parent_face_ids']);changed=changed[changed<nf[1]];allowed=np.asarray(spec['allowed_parent_face_ids'],int);assert set(changed)<=set(allowed),'Native parent changed outside declared ostium patch'
 v=arrays['vertices'];f=arrays['faces'];roles=arrays['source_component_ids'];owner=(roles!=0).astype(int);assert set(owner)=={0,1}
 normals=arrays['corner_normal_words'].copy();uvs=arrays['corner_uvs'];seam=[]
 for vi in np.unique(f):
  where=np.argwhere(f==vi)
  if len(set(owner[where[:,0]]))<2:continue
  value=normals[where[:,0],where[:,1]].astype(float).sum(0);assert np.linalg.norm(value)>0;word=np.rint(value/np.linalg.norm(value)*32767).astype(np.int16)
  for fi,ci in where:
   if roles[fi]==0:
    pid=int(arrays['source_parent_face_ids'][fi]);assert pid in set(allowed),'Shared normal change reaches undeclared native parent chart'
   normals[fi,ci]=word
  seam.append(int(vi))
 # One dyadic storage origin preserves literal native Float32 controls exactly.
 # Center X on a nearby dyadic value so new ostium intersections retain
 # sub-nanometre precision. Native Float32 X controls remain literal under
 # exact subtraction; an absolute-X encoding lost precision at this fork.
 origin=np.array([.1015625 if side=='left'else-.1015625,1.5,0.]);stored=(v-origin).astype('<f4');actual=stored.astype(float)+origin;compiled=Mesh(actual,f,FRAME)
 check=audit(compiled,require_closed=False,expected_components=1,max_pairs=600000,max_tests=2000000);write(out/'compiled-merged-audit.json',check);assert check['passed'];loops=boundary_loops(f);assert len(loops)==4 and topology(compiled)['components'][0]['euler']==-2
 for p in arrays['source_vertices'][:len(parent.vertices)]:
  matching=np.flatnonzero(np.all(v==p,axis=1))
  assert all(np.array_equal(actual[i],p)for i in matching),'Retained native source control moved on storage'
 capv=actual.tolist();capf=f.tolist()
 for loop in loops:
  center=len(capv);capv.append(actual[loop].mean(0).tolist());capf.extend([[center,int(loop[(j+1)%len(loop)]),int(loop[j])]for j in range(len(loop))])
 classifier=Mesh(np.asarray(capv),np.asarray(capf),FRAME);c=audit(classifier,require_closed=True,expected_components=1,expected_genus=[0],max_pairs=600000,max_tests=2000000);write(out/'diagnostic-connected-blood-classifier-audit.json',c);assert c['passed']
 # Original open-mouth centers belong to the diagnostic cap plane, not the
 # blood interior. Bind these endpoint witnesses explicitly instead of making
 # a false strict-interior assertion or accepting arbitrary boundary probes.
 supplied=np.asarray(spec['blood_probe_points_m']);mouth_centers=np.asarray([actual[loop].mean(0)for loop in loops]);gaps=np.linalg.norm(supplied[:,None]-mouth_centers[None],axis=2);endpoint=np.min(gaps,axis=1)<1e-7
 endpoint_rows=[dict(probe_id=int(i),original_mouth_index=int(np.argmin(gaps[i])),literal_mouth_center_gap_m=float(np.min(gaps[i])))for i in np.flatnonzero(endpoint)]
 probes=np.r_[arrays['communication_probe_points_m'],supplied[~endpoint]];inside=winding_contains(classifier,probes);assert all(inside['inside'])and not any(inside['boundary_ambiguous']);write(out/'explicit-blood-probes-and-endpoint-witnesses.json',dict(interior_probe_count=len(probes),all_strictly_interior=True,original_open_mouth_endpoint_witnesses=endpoint_rows,endpoint_match_tolerance_m=1e-7,temporary_caps_diagnostic_only=True));classifier.save(out/'diagnostic-only-blood-classifier.npz')
 data=bytearray();after=copy.deepcopy(atlas);changedrecords=[]
 def append(a):
  a=np.ascontiguousarray(a)
  while len(data)%a.dtype.itemsize:data.append(0)
  offset=len(data);data.extend(a.tobytes());return offset
 for group,pid in [(0,PARENT),(1,CHILD)]:
  faceids=np.flatnonzero(owner==group);cp=stored[f[faceids]].reshape(-1,3);cn=normals[faceids].reshape(-1,3);cu=uvs[faceids].reshape(-1,2);dtype=np.dtype([('p','<f4',(3,)),('n','<i2',(3,)),('u','<f4',(2,))]);keys=np.empty(len(cp),dtype);keys['p'],keys['n'],keys['u']=cp,cn,cu;_,first,indices=np.unique(keys,return_index=True,return_inverse=True);q,n,uv=cp[first],cn[first],cu[first];ff=indices.reshape(-1,3).astype('<u4');record=next(p for p in after['parts']if p['id']==pid)if group==0 else dict(id=pid,conceptId=pid,name=f'Clavicular branch of thoraco-acromial artery ({side})',groups=parent_record['groups'],system='arterial',material=parent_record['material'],sourceId=PARENT.removeprefix('ZA:')+':A12.2.09.008',license='CC BY-SA 4.0',provenance=dict(label='Independently authored illustrative course with source-derived native parent ostium',url='https://github.com/Z-Anatomy/Models-of-human-anatomy',detail='Native Z-Anatomy parent ostium is source-derived; remainder independently constructed. CA2024 coarse course evidence used only for orientation. No CA2024 mesh, normals, UVs or textures copied.'))
  for key in ['sourceOffset','geometryScale','normalSigns','geometryInstance','pairedTissue','sectionSolid','tissueWallAssembly','tissueWallLayer','tissueDrawGroups']:record.pop(key,None)
  record.update(positions=append(q),normals=append(n),uvs=append(uv),indices=append(ff),chunk=len(atlas['chunks']),geometryOrigin=origin.tolist(),vertexCount=len(q),indexCount=int(ff.size),bounds=[actual[f[faceids]].min((0,1)).tolist(),actual[f[faceids]].max((0,1)).tolist()],sourceOuterSurfaceOnly=True,nativeArterialBranchAssembly=dict(schema='independent-clavicular-parent-ostium/v1',sourceCatalogueSha256=spec['source_catalogue_sha256'],role='parent ostium'if group==0 else 'clavicular main and muscular twig',courseIllustrativeVariant=True,finiteWall=False,perfusionOrHistologyComplete=False),coverageLimitation='Illustrative macroscopic arterial exterior course with a communicating parent ostium and open distal continuations. Variant origin/recipient territories remain provisional; no finite wall, histology or complete perfusion claim.')
  if group:record.update(ta98Term=TERM,uvMethod='Procedural longitudinal/circumferential coordinates on independently constructed tubes');after['parts'].append(record)
  changedrecords.append(record);np.savez_compressed(out/('parent-owned.npz'if group==0 else 'clavicular-owned.npz'),vertices=actual,faces=f[faceids],merged_face_ids=faceids,corner_normal_words=normals[faceids],corner_uvs=uvs[faceids])
 bank=out/'clavicular-exterior.bin';bank.write_bytes(data);after['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=len(data)));after['concepts'].append(dict(id=CHILD,name=f'Clavicular branch of thoraco-acromial artery ({side})',elements=[CHILD],partOf=PARENT,parentId=PARENT,anatomicalParentId=PARENT,ta98Term=TERM,ta98Kind='exact',coverageLimitation='Independently authored illustrative macroscopic variant course. Recipient supply territories and exact origin are provisional; distal continuations remain open. No finite wall, histology or complete perfusion claim.'));after['ta98TermBindings'].append(dict(conceptId=CHILD,partIds=[CHILD],term=TERM,kind='exact'));after['partsCount']=len(after['parts']);after['triangles']=int(atlas['triangles']-len(parent.faces)+len(f));write(out/'proposed-atlas.json',after)
 for rec in changedrecords:
  m=AtlasReader().read_part(after,rec,frame=FRAME);m.save(out/('parent-packed.npz'if rec['id']==PARENT else 'clavicular-packed.npz'))
 np.savez_compressed(out/'compiled-shared-geometry.npz',vertices=actual,faces=f,face_owner=owner,corner_normal_words=normals,corner_uvs=uvs,shared_interface_control_ids=np.asarray(seam))
 already_bound=any(x.get('term')==TERM for x in atlas['ta98TermBindings'])or any(x.get('ta98Term')==TERM for x in atlas['concepts'])
 write(out/'packing-proof.json',dict(complete=True,passed=True,source_catalogue_sha256=spec['source_catalogue_sha256'],proposed_catalogue_sha256=H(out/'proposed-atlas.json'),bank_sha256=H(bank),bank_bytes=len(data),changed_parts=[PARENT,CHILD],new_TA98_terms=[]if already_bound else[TERM],represented_term=TERM,construction_side=side,parent_changed_faces=changed,allowed_parent_faces=allowed,shared_interface_controls=len(seam),all_faces_owned_once=True,all_four_outer_mouths_open=True,temporary_visible_classifiers=0,derived_cones_visible=0,closed_blood_classifier_components=1,sourceOuterSurfaceOnly=True,all_outside_parts_literal=True,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--course-spec',type=Path,required=True);p.add_argument('--approval',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--side',choices=['left','right'],default='left');a=p.parse_args();run(a.course_spec,a.approval,a.output,a.side)
