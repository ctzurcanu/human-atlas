"""Bounded, read-only native-source inventory and physical-port diagnostics."""
from pathlib import Path
import json,sqlite3,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.ta98_authority import corrected_terms
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/posterior-intercostal-native-v1'
C=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME='model:local-ta98:assembled-metres'
SOURCES=[f'ZA:{name}.{side}'for side in ['l','r']for name in ['First posterior intercostal artery','Second posterior intercostal artery','Posterior intercostal arteries']]
REFS=[dict(url='https://ifaa.unifr.ch/Public/EntryPage/PDF/TA98%20Chapter%2012.pdf',scope='Authoritative .007-.019 inventory and .011 footnote: distinct foraminal routes and root/cord targets'),dict(url='https://medicine.uams.edu/neuroscience/education/medical-school-courses/human-structure-module/anatomy-tables/artery-tables/arteries-of-the-thorax/',scope='First two courses arise from supreme intercostal; third through eleventh from thoracic aorta; recipient classes'),dict(url='https://pubmed.ncbi.nlm.nih.gov/21819734/',scope='Primary anatomical variation study; no universal identical dorsal/spinal branch pattern assumed')]

def loaded():
 a=json.loads(C.read_bytes());return a,{p['id']:p for p in a['parts']},AtlasReader()

def inventory():
 started=time.monotonic();B.mkdir(exist_ok=True,parents=True);assert not(B/'inventory.json').exists();a,parts,reader=loaded()
 db=sqlite3.connect(ROOT/'.local-models/ta98-authority.sqlite');db.row_factory=sqlite3.Row
 terms=corrected_terms({r['id']:dict(r)for r in db.execute('select * from ta98')});ids={'A12.2.11.007'}
 while True:
  children={q for q,t in terms.items()if t.get('parent_id')in ids}
  if children<=ids:break
  ids|=children
 assert len(ids)==13
 census_path=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-blood-controls-v2/census-installed.json';census=json.loads(census_path.read_bytes());resolved={r['id']:r for r in census['terms']}
 routes={8:'Dorsal course from posterior intercostal stem towards deep back; physical takeoff not yet established',9:'Medial skin/perforating recipient; dorsal branch must exist first',10:'Lateral back skin recipient; not the separate .018 lateral thoracic branch',11:'Intervertebral-foramen channel; avoid fixed pedicles and accompany existing nerve route',12:'Canal-side route behind vertebral body, distinct from prelaminar surface',13:'Canal-side route anterior to lamina, distinct from postcentral surface',14:'Posterior spinal root/rootlet recipient; no mandatory cord continuation',15:'Anterior spinal root/rootlet recipient; no mandatory cord continuation',16:'Selected level-dependent contribution to actual anterior spinal artery; do not repeat at every level without evidence',17:'Collateral course near costal angle along upper margin of next lower rib',18:'Lateral thoracic skin takeoff from intercostal course near lateral chest',19:'Lateral mammary recipient branches; source and organ ports require sex-specific verification'}
 sources=[]
 for pid in SOURCES:
  p=parts[pid];m=reader.read_part(a,p,frame=FRAME)
  sources.append(dict(id=pid,name=p['name'],geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,bank_sha256=m.source['bank_sha256'],bank_path=m.source['bank_path'],source_record=m.source['record'],expected_scope='Unbranched course carriers; named child branch identities not inferred',license=p.get('license'),provenance=p.get('provenance')))
 rows=[]
 for code in sorted(ids):
  t=terms[code];r=resolved[code];num=int(code.split('.')[-1]);rows.append(dict(authority=t,representation_class='whole-organ',representation_kind='vessel-or-vascular-network',current_direct_parts=r['resolved_direct_part_ids'],existing_source_parts=SOURCES if num==7 else[],named_source_candidates=r['named_but_unbound_source_candidates'],missing_distinct_representation=num!=7,source_geometry_anatomically_accepted=False,planned_interface=routes.get(num,'Retained bilateral first-to-eleventh course assembly'),new_generic_branches_allowed=False))
 landmarks=[]
 for pid in ['ZA:Body of sternum:Bone-7','ZA:Vertebra T6:Bone-5','ZA:Vertebra T9:Bone-5','ZA:Sacrum:Bone-5']:
  m=reader.read_part(a,parts[pid],frame=FRAME);v=m.vertices;landmarks.append(dict(id=pid,geometry_sha256=m.geometry_sha,bounds_m=np.array([v.min(0),v.max(0)]).tolist(),vertex_mean_m=v.mean(0).tolist()))
 sternum=np.array(landmarks[0]['vertex_mean_m']);back=np.mean([landmarks[1]['vertex_mean_m'],landmarks[2]['vertex_mean_m']],axis=0);direction=sternum-back;direction[0:2]=0;direction/=np.linalg.norm(direction)
 axis=dict(frame=FRAME,landmarks=landmarks,anterior_axis=direction.tolist(),posterior_axis=(-direction).tolist(),superior_axis=[0,1,0],left_axis=[1,0,0],basis='Native fixed sternal body lies wholly anterior of native T6/T9 vertebral bones along Z; anatomical sign is not inferred from a camera preset',source_catalogue_sha256=digest_file(C),geometry_changed=False)
 atomic_json(B/'axis-proof.json',axis)
 inv=dict(schema='posterior-intercostal-complete-authority-inventory/v1',source_catalogue=str(C),source_catalogue_sha256=digest_file(C),authority_sha256=digest_file(ROOT/'.local-models/ta98-authority.sqlite'),census_sha256=digest_file(census_path),root='A12.2.11.007',inclusive_term_count=13,current_direct_terms=1,missing_named_branch_terms=12,terms=rows,native_sources=sources,frame=FRAME,skeleton_fixed=True,geometry_changed=False,master_changed=False,wall_lumen_and_anatomical_acceptance=False,references=REFS,elapsed_seconds=time.monotonic()-started)
 atomic_json(B/'inventory.json',inv)
 report=ROOT/'reports/mesh-audit/REVIEW.md';marker='### Posterior intercostal artery source inventory — 2026-10-26'
 assert marker not in report.read_text()
 with report.open('a')as f:
  f.write('\n\n'+marker+'\n\nRead-only native-source audit; skeleton and every parent/recipient stay fixed. Full authoritative A12.2.11.007 closure contains 13 vessel/network terms. Only the whole posterior intercostal family currently resolves directly; all 12 named descendant branches remain unrepresented. Six ZA carriers must be inspected for actual serial courses before any construction. Coarse whole-course aliases do not establish children. The first two courses have a supreme-intercostal parent, while courses 3–11 have a thoracic-aortic parent. Root contacts alone do not prove an ostium or lumen. Local evidence: `'+str(B.relative_to(ROOT))+'/inventory.json`.\n\n| TA98 | Parent | Term | Existing distinct geometry | Required work |\n|---|---|---|---|---|\n')
  for row in rows:
   t=row['authority'];f.write('| '+t['id']+' | '+str(t.get('parent_id'))+' | '+t['name_en']+' | '+('whole-course source carriers'if t['id'].endswith('.007')else'absent')+' | '+row['planned_interface']+' |\n')
 print(json.dumps(dict(terms=13,current_direct=1,native_carriers=len(sources),missing_named_branches=12,anterior_axis=direction.tolist(),elapsed_seconds=inv['elapsed_seconds'])))

def courses(side):
 from geometry_tools.algorithms.tube_station_topology import recover_station_rings
 started=time.monotonic();a,parts,reader=loaded();inventory=json.loads((B/'inventory.json').read_bytes());bindings={p['id']:p for p in inventory['native_sources']};arrays={};rows=[]
 for pid in [p for p in SOURCES if p.endswith('.'+side)]:
  m=reader.read_part(a,parts[pid],frame=FRAME);assert m.geometry_sha==bindings[pid]['geometry_sha256']and m.attributes_sha==bindings[pid]['attributes_sha256']
  recovered=recover_station_rings(m.vertices,m.faces)
  for c in recovered:
   index=len(rows);rings=np.array(c['rings']);centres=m.vertices[rings].mean(1);prefix=f'course_{index:02}'
   arrays[prefix+'_rings']=rings;arrays[prefix+'_centres']=centres;arrays[prefix+'_source_face_ids']=c['face_ids']
   ends=[]
   for endpoint in [0,-1]:
    ring=rings[endpoint];v=m.vertices[ring];center=v.mean(0);vector=np.cross(v-center,np.roll(v-center,-1,axis=0)).sum(0)
    ends.append(dict(source_vertex_ids=ring.tolist(),world_points_m=v.tolist(),center_m=center.tolist(),section_area_m2=float(np.linalg.norm(vector)/2),diameter_m=float(np.linalg.norm(v[:,None]-v[None,:],axis=2).max()),normal=vector.tolist()))
   rows.append(dict(index=index,source_id=pid,source_component=c['source_component'],source_geometry_sha256=m.geometry_sha,source_attributes_sha256=m.attributes_sha,station_count=len(rings),section_controls=len(rings[0]),source_faces=len(c['face_ids']),root_endpoint_not_assigned=True,source_mouths=ends,course_length_m=float(np.linalg.norm(np.diff(centres,axis=0),axis=1).sum()),source_course_has_fork=False,named_descendant_branch_identity_inferred=False,centreline_scope='Means of actual topologically recovered source rings; no new path, fitting or physiological lumen is inferred'))
 arrays['frame']=np.array(FRAME);atomic_npz(B/f'{side}-courses.npz',**arrays)
 result=dict(schema='posterior-intercostal-literal-native-courses/v1',side=side,source_catalogue_sha256=digest_file(C),source_inventory_sha256=digest_file(B/'inventory.json'),courses=rows,course_count=len(rows),all_recovered_from_actual_source_connectivity=True,all_named_branches_absent_from_these_unbranched_carriers=True,atlas_changed=False,elapsed_seconds=time.monotonic()-started)
 assert len(rows)==11
 atomic_json(B/f'{side}-courses.json',result);print(json.dumps(dict(side=side,courses=len(rows),stations=[r['station_count']for r in rows],ring_controls=sorted({r['section_controls']for r in rows}),elapsed_seconds=result['elapsed_seconds'])))

def ports(side):
 from geometry_tools.spatial import BVH
 from geometry_tools.algorithms.tube_station_topology import recover_station_rings
 started=time.monotonic();a,parts,reader=loaded();inventory=json.loads((B/'inventory.json').read_bytes());data=json.loads((B/f'{side}-courses.json').read_bytes());bank=np.load(B/f'{side}-courses.npz',allow_pickle=False)
 bindings={p['id']:p for p in inventory['native_sources']};meshes={};bvhs={}
 def mesh(pid):
  if pid not in meshes:meshes[pid]=reader.read_part(a,parts[pid],frame=FRAME)
  return meshes[pid]
 def bvh(pid):
  if pid not in bvhs:bvhs[pid]=BVH(mesh(pid))
  return bvhs[pid]
 ordinal=['First','Second','Third','Fourth','Fifth','Sixth','Seventh','Eighth','Ninth','Tenth','Eleventh','Twelfth']
 # The retained intercostal nerve roots are potential foraminal guide ports,
 # not surrogate meshes for separate anterior/posterior radicular arteries.
 nerve_id=f'ZA:Intercostal nerves.{side}';nerve_ports=[];nerve_error=None
 if nerve_id in parts:
  try:
   nm=mesh(nerve_id)
   for c in recover_station_rings(nm.vertices,nm.faces):
    # Both posterior root and anterior sternal end can be medial. The fixed
    # sternum/vertebral proof, rather than absolute X, identifies posterior.
    first,last=c['rings'][0],c['rings'][-1];r=first if nm.vertices[first].mean(0)[2]<nm.vertices[last].mean(0)[2]else last
    nerve_ports.append(dict(source_vertex_ids=r.tolist(),points_m=nm.vertices[r].tolist(),center_m=nm.vertices[r].mean(0).tolist(),geometry_sha256=nm.geometry_sha,source_component=c['source_component']))
   nerve_ports.sort(key=lambda q:-q['center_m'][1])
  except ValueError as e:nerve_error=str(e)
 root_rows=[]
 for c in data['courses']:
  pid=c['source_id'];m=mesh(pid);assert m.geometry_sha==bindings[pid]['geometry_sha256']and m.attributes_sha==bindings[pid]['attributes_sha256']
  parent=f'ZA:Supreme intercostal artery.{side}'if'First posterior'in pid or'Second posterior'in pid else'TA98WALL:male:thoracic-aorta:externa'
  if parent not in parts:raise ValueError('Actual installed parent missing: '+parent)
  ends=[]
  for endpoint,q in enumerate(c['source_mouths']):
   points=np.asarray(q['world_points_m']);nearest,evidence=bvh(parent).closest_points(points,backend='bvh');d=nearest['distance_m']
   ends.append(dict(endpoint=endpoint,minimum_sample_distance_m=float(d.min()),mean_sample_distance_m=float(d.mean()),maximum_sample_distance_m=float(d.max()),closest_parent_points_m=nearest['points'].tolist(),closest_parent_face_ids=nearest['face_ids'].tolist(),query=evidence))
  root_end=int(np.argmin([q['mean_sample_distance_m']for q in ends]));q=c['source_mouths'][root_end];root=np.asarray(q['center_m']);shared=sum(any(np.array_equal(x,y)for y in mesh(parent).vertices)for x in np.asarray(q['world_points_m']))
  centres=bank[f"course_{c['index']:02}_centres"].copy()
  if root_end:centres=centres[::-1]
  level=1 if 'First posterior' in pid else 2 if 'Second posterior' in pid else None
  root_rows.append(dict(course_index=c['index'],source_id=pid,source_component=c['source_component'],source_geometry_sha256=m.geometry_sha,source_attributes_sha256=m.attributes_sha,expected_parent=parent,parent_geometry_sha256=mesh(parent).geometry_sha,parent_attributes_sha256=mesh(parent).attributes_sha,root_endpoint=root_end,root_mouth=q,distal_mouth=c['source_mouths'][1-root_end],parent_endpoint_queries=ends,shared_exact_parent_root_controls=int(shared),communicating_ostium_verified=False,root_is_actual_source_boundary_ring=True,level=level,level_assignment_basis='Explicit source name'if level else'Pending root-height order and fixed rib review',centres=centres))
 group=sorted([q for q in root_rows if q['level']is None],key=lambda q:-q['root_mouth']['center_m'][1]);assert len(group)==9
 for level,q in enumerate(group,3):q['level']=level;q['level_assignment_basis']='Superior-to-inferior order of nine aortic source root rims; fixed matching rib context listed for review, not accepted anatomy'
 for q in root_rows:
  level=q['level'];centres=q.pop('centres');lateral=int(np.argmax(abs(centres[:,0])));probe=centres[lateral]
  rib_candidates=[p for p in parts if p.startswith(f'ZA:{ordinal[level-1]} rib.{side}:Bone')];assert len(rib_candidates)==1
  rib=rib_candidates[0];near,evidence=bvh(rib).closest_points(probe[None],backend='bvh')
  q['lateral_course_landmark']=dict(station_from_root=lateral,actual_source_center_m=probe.tolist(),fixed_rib_id=rib,fixed_rib_geometry_sha256=mesh(rib).geometry_sha,closest_rib_face=int(near['face_ids'][0]),closest_rib_point_m=near['points'][0].tolist(),center_to_rib_surface_distance_m=float(near['distance_m'][0]),anatomical_branch_port_exists=False,scope='Actual lateral-most course station and fixed rib surface; not a fabricated .018 branch origin or certified rib-groove fit')
  vertebrae=[p for p in parts if p.startswith(f'ZA:Vertebra T{level}:Bone')or p.startswith(f'ZA:Vertebra T{level+1}:Bone')]
  q['fixed_foramen_boundaries']=[dict(id=p,geometry_sha256=mesh(p).geometry_sha,bounds_m=np.array([mesh(p).vertices.min(0),mesh(p).vertices.max(0)]).tolist())for p in vertebrae]
  if len(nerve_ports)==11:
   port=nerve_ports[level-1];centre=np.asarray(port['center_m']);station=int(np.argmin(np.linalg.norm(centres-centre,axis=1)))
   q['foraminal_nerve_guide']=dict(id=nerve_id,root_port=port,root_assignment_basis='Actual posterior-most source boundary endpoint under fixed +Z-anterior landmark proof; level numbering by descending root height',closest_actual_arterial_station_from_root=station,arterial_station_center_m=centres[station].tolist(),centre_distance_m=float(np.linalg.norm(centres[station]-centre)),physical_spinal_branch_port_exists=False,foramen_canal_clearance_verified=False,scope='Existing posterior source nerve boundary; numbering remains inferred by height. Not separate anterior/posterior rootlets or a free-space foramen certificate.')
  else:q['foraminal_nerve_guide']=dict(source_id=nerve_id,recovered_source_components=len(nerve_ports),error=nerve_error,physical_spinal_branch_port_exists=False)
  q['new_named_branch_ports_verified']=False
 spinal=[dict(id=p['id'],bounds=p.get('bounds'),source_name=p.get('name'),suppressed=bool(p.get('suppressed')),default_hidden=p['id']in set(a['defaultHidden']))for p in a['parts']if p['system']=='arterial'and any(w in p.get('name','').lower()for w in ['anterior spinal','posterior spinal','radicular','segmental medullary'])]
 proof=dict(schema='posterior-intercostal-actual-source-port-audit/v1',side=side,source_catalogue_sha256=digest_file(C),inventory_sha256=digest_file(B/'inventory.json'),course_bank_sha256=digest_file(B/f'{side}-courses.npz'),frame=FRAME,axis_proof_sha256=digest_file(B/'axis-proof.json'),courses=sorted(root_rows,key=lambda q:q['level']),fixed_nerve_source=nerve_id,fixed_nerve_source_components=len(nerve_ports),spinal_recipient_sources=spinal,all_named_descendant_branch_geometry_missing=True,all_13_terms_represented=False,source_parent_rims_are_physical=True,source_parent_ostia_verified=False,full_course_packing_not_audited=True,skeleton_fixed=True,geometry_modified=False,new_geometry_created=False,elapsed_seconds=time.monotonic()-started)
 proof['supersedes_unverified_medial_nerve_endpoint_inference']=str(B/f'{side}-physical-ports.json')
 atomic_json(B/f'{side}-physical-ports-v2.json',proof)
 print(json.dumps(dict(side=side,source_roots=len(root_rows),fixed_nerve_components=len(nerve_ports),nerve_recovery_error=nerve_error,root_gaps_mm=[round(q['parent_endpoint_queries'][q['root_endpoint']]['mean_sample_distance_m']*1000,5)for q in proof['courses']],exact_shared_root_controls=[q['shared_exact_parent_root_controls']for q in proof['courses']],elapsed_seconds=proof['elapsed_seconds'])))

if __name__=='__main__':
 if sys.argv[1]=='inventory':inventory()
 elif sys.argv[1]=='courses':courses(sys.argv[2])
 elif sys.argv[1]=='ports':ports(sys.argv[2])
