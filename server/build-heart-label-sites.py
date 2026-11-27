"""Compact, explicitly qualified cardiac reference points and curves.

These controls are NOT additional tissue meshes or histological layers. The
existing 42 geometry-backed term controls are counted separately. Every new
label has an individual source/derivation, carrier, pose binding and limitation.
"""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
import igl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.geometry_instance import instance_pose
from geometry_tools.common import atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--reference-stage',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();started=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
reference_path=h/args.reference_stage/'report.json';reference=json.loads(reference_path.read_text());ref={r['cid']:r for r in reference['landmarks']};reader=AtlasReader();paths={s:ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{s}/atlas.json'for s in ['male','female']};atlases={s:json.loads(p.read_text())for s,p in paths.items()};male=atlases['male'];marker=male['finiteCardiacIntegration'];parts={p['id']:p for p in male['parts']};groups={g['key']:g for g in marker['termGroups']};terms=json.loads((h/'male-cardiopulmonary-final-install-v525/inventory.json').read_text())['terms'];authority={r['authority']['id']:r['authority']for r in terms if r['authority']['id'].startswith('A12.1.')};assert len(authority)==113
existing={g['term']for g in marker['termGroups']if g['kind']=='exact'};assert len(existing)==42
unfinished={'A12.1.00.008':'A distinct apical cleft remains unconstructed.','A12.1.00.024':'Variable false chordae remain unconstructed.','A12.1.00.028':'A distinct infundibular tendon remains unconstructed.','A12.1.00.029':'A distinct tendon of the caval valve remains unconstructed.','A12.1.00.030':'Koch triangle: bounded endocardial surface remains provisional; the TA98 sinu-atrial wording is corrected to atrioventricular in later terminology.','A12.1.01.004':'Individual microscopic venous channels remain unconstructed.','A12.1.01.007':'Closed adult site; no patent fetal opening is added.','A12.1.01.014':'Superior interatrial-fold reference; a discrete tubercle is not established.','A12.1.04.006':'Commissural leaflet regions remain unconstructed.','A12.1.08.004':'A separate bronchopericardial sheet remains unconstructed.','A12.1.08.008':'A separate serosal stratum remains unconstructed.','A12.1.08.009':'A separate subserosal stratum remains unconstructed.','A12.1.08.010':'The variable left caval fold remains unconstructed.'}
new_terms=sorted(set(authority)-existing);assert len(new_terms)==71
meshes={pid:reader.read_part(male,parts[pid],frame='assembled-metres',catalogue_path=paths['male'])for pid in set(sum([marker.get(k,[])for k in ['wallPartIds','valvePartIds','papillaryPartIds','rootPartIds']],[]))}
def ids(key):return groups[key]['partIds']
def center(pids):
 points=[];weights=[]
 for pid in pids:
  t=meshes[pid].triangles;w=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1);points.append(w@t.mean(1)/w.sum());weights.append(w.sum())
 return np.average(points,axis=0,weights=weights)
def c(key):return center(ids(key))
def landmark(cid):return np.array(ref[cid]['centroid_m'])
ra='TA98HEART:male:myocardium:right-atrium';la='TA98HEART:male:myocardium:left-atrium';rv='TA98HEART:male:myocardium:right-ventricle';lv='TA98HEART:male:myocardium:left-ventricle';endo='TA98HEART:male:endocardium';epi='TA98HEART:male:epicardium';ias='TA98HEART:male:myocardium:interatrial-septum';ivs='TA98HEART:male:myocardium:interventricular-muscular-septum'
def root_point(pid,atrium):
 m=meshes[pid];host=meshes[atrium];_,_,q=igl.point_mesh_squared_distance(m.vertices,host.vertices,host.faces);distance=np.linalg.norm(q-m.vertices,axis=1);return m.vertices[np.argmin(distance)]
svc=root_point('ZA:Superior vena cava',ra);ivc=root_point('TA98HEPATIC:male:region:A12.3.09.001',ra);tv=c('tricuspid-valve');mv=c('mitral-valve');pv=c('pulmonary-valve');av=c('aortic-valve');cs=.55*ivc+.45*tv;cs[2]-=.004
def point_on(pid,q):
 m=meshes[pid];d,face,point=igl.point_mesh_squared_distance(np.asarray(q,float).reshape(-1,3),m.vertices,m.faces);return point[0],dict(partId=pid,face=int(face[0]),projectionDistanceMetres=float(np.sqrt(d[0])))
def nearby(pid,q):return point_on(pid,q)[0]
spec={}
def add(term,carrier,q,kind,derivation,cids=(),points=None):
 assert term in new_terms and term not in spec
 spec[term]=dict(carrier=carrier,point=np.array(q,float),representation=kind,derivation=derivation,referenceCids=list(cids),points=None if points is None else np.asarray(points,float))
def sourced(term,cid,carrier,kind='anatomical-point'):
 q,proof=point_on(carrier,landmark(cid));add(term,carrier,q,kind,'Named local reference landmark projected onto the current carrier; source-to-carrier registration is provisional.',[cid]);spec[term]['projection']=proof
epicardial_features={g['term']:np.array(g['anchor']['point'])for g in male['cardiacFeatureIntegration']['termGroups']};apex=epicardial_features['A12.1.00.007']
# Regions and boundary curves use actual chamber ownership as their carrier.
right_edge=meshes[ra].vertices;edge=right_edge[np.argmin(right_edge[:,0])];add('A12.1.00.006',ra,edge,'boundary-curve','Rightmost native right-atrial border, with provisional curve extent.',points=[edge+np.array([0,.006,0]),edge,edge-np.array([0,.006,0])])
for term,anterior in [('A12.1.00.009',True),('A12.1.00.010',False)]:
 septal=c('interventricular-septum');q=.55*apex+.45*septal;q[2]+=.018 if anterior else -.018;q=nearby(epi,q);add(term,epi,q,'boundary-curve','Interventricular ownership boundary projected to the anterior/posterior coat; groove relief is not added.',points=[nearby(epi,.25*apex+.75*septal+np.array([0,0,.018 if anterior else -.018])),q,apex])
q=nearby(epi,(tv+mv)/2+np.array([0,0,.014]));add('A12.1.00.011',epi,q,'boundary-curve','Atrioventricular surface junction guided by the two measured annuli.',points=[nearby(epi,tv),q,nearby(epi,mv)])
q=(landmark(20454)+landmark(23809))/2;add('A12.1.00.018',ra,nearby(ra,landmark(20454)),'collective','Collective reference for the right and left auricles.',[20454,23809],points=[nearby(ra,landmark(20454)),nearby(la,landmark(23809))])
add('A12.1.00.020',endo,nearby(endo,(landmark(22800)+landmark(18687))/2),'collective','Apical trabecular reference sites; fine trabecular relief remains unconstructed.',[22800,18687],points=[nearby(endo,landmark(22800)),nearby(endo,landmark(18687))])
add('A12.1.00.021',lv,nearby(lv,apex),'anatomical-point','Apical myocardial vortex location; individual helical muscle fibres are not resolved.')
add('A12.1.00.022',marker['papillaryPartIds'][0],center(marker['papillaryPartIds']),'collective','Collective control for the five actual stored native papillary tissues.',[18980,20592,18331,24354,23271],points=[center([pid])for pid in marker['papillaryPartIds']])
for term,cid,carrier in [('A12.1.00.025',20644,'TA98HEART:male:fibrous:right-annulus'),('A12.1.00.026',23777,'TA98HEART:male:fibrous:left-annulus'),('A12.1.01.002',20454,ra),('A12.1.01.005',22621,ias),('A12.1.01.006',21663,ias),('A12.1.01.008',18054,ra),('A12.1.02.008',16459,rv),('A12.1.02.020',19988,rv),('A12.1.02.021',18687,rv),('A12.1.03.002',23809,la),('A12.1.04.007',23271,marker['papillaryPartIds'][-1]),('A12.1.04.008',24354,'ZA:Inferior papillary muscle of left ventricle'),('A12.1.04.011',22800,lv)]:sourced(term,cid,carrier,'regional-reference')
# Caval/valve mouths are annotated spaces and add no cap or plug.
for term,carrier,q,derivation in [('A12.1.01.009',ra,cs,'Posteroinferior right-atrial coronary-sinus ostium reference; exact recipient join is not established.'),('A12.1.01.010',ra,ivc,'Measured closest native IVC/right-atrial interface.'),('A12.1.01.011',ra,svc,'Measured native SVC/right-atrial interface.'),('A12.1.02.002','TA98HEART:male:fibrous:right-annulus',tv,'Right AV opening bounded by the measured fibrous support and tricuspid leaflets.'),('A12.1.02.009',rv,pv,'Pulmonary outflow opening at the actual semilunar apparatus.'),('A12.1.04.002','TA98HEART:male:fibrous:left-annulus',mv,'Left AV opening bounded by the measured fibrous support and mitral leaflets.'),('A12.1.04.010',lv,av,'Aortic outflow opening at the actual semilunar apparatus.')]:add(term,carrier,q,'opening-rim',derivation)
fossa=spec['A12.1.01.005']['point'];pectinate=spec['A12.1.01.008']['point'];crest=nearby(ra,.55*(svc+ivc)/2+.45*pectinate)
add('A12.1.01.003',ra,crest,'boundary-curve','Terminal-crest reference between the caval sinus and pectinate region; distinct crest relief is not constructed.',points=[nearby(ra,svc),crest,nearby(ra,ivc)])
add('A12.1.01.012',ra,nearby(ra,(svc+ivc)/2),'regional-reference','Posterior smooth right-atrial region between the two caval mouths.')
add('A12.1.01.013',ra,nearby(ra,crest+np.array([0,0,-.004])),'boundary-curve','External terminal-sulcus reference corresponding to the terminal crest.')
add('A12.1.01.015',ra,nearby(ra,.75*ivc+.25*fossa),'anatomical-point','Inferior caval valve location; a distinct valve leaflet has not been constructed.')
add('A12.1.01.016',ra,nearby(ra,.8*cs+.2*tv),'anatomical-point','Coronary-sinus valve location; a distinct valve leaflet has not been constructed.')
add('A12.1.02.007',rv,nearby(rv,.6*pv+.4*tv),'regional-reference','Supraventricular-crest region between inflow and pulmonary outflow; distinct crest relief remains provisional.')
add('A12.1.03.003',la,nearby(la,landmark(23809)+np.array([0,-.003,0])),'regional-reference','Left-auricular pectinate region; fine pectinate relief is not separately constructed.',[23809])
veins=[pid for pid in marker['rootPartIds']if 'pulmonary vein'in parts[pid]['name'].lower()];mouths=[root_point(pid,la)for pid in veins];add('A12.1.03.004',la,mouths[0],'opening-rim','Four actual pulmonary-vein/left-atrial interface locations; no opening disks.',points=mouths)
add('A12.1.03.005',ias,nearby(ias,fossa+np.array([0,0,-.001])),'regional-reference','Left septal flap/closure reference opposite the fossa; no patent adult foramen is implied.')
add('A12.1.04.009',lv,nearby(lv,.7*av+.3*mv),'regional-reference','Left-ventricular outflow vestibule adjoining the actual aortic apparatus.')
# Semilunar fine anatomy is a point/curve reference on the real leaflet
# carriers. It is not silently counted as a constructed microscopic nodule.
for key,base in [('pulmonary-valve','A12.1.02'),('aortic-valve','A12.1.04')]:
 valve=ids(key);mid=center(valve);nodes=[];lunules=[];commissures=[]
 for pid in valve:
  m=meshes[pid];d=np.linalg.norm(m.vertices-mid,axis=1);node=m.vertices[d.argmin()];nodes.append(node);target=node+.002*(center([pid])-node)/max(np.linalg.norm(center([pid])-node),1e-12);lunules.append(nearby(pid,target));commissures.append(m.vertices[d.argmax()])
 numbers=(14,15,16)if base.endswith('02')else(16,17,18)
 for number,points,kind,name in zip(numbers,[nodes,lunules,commissures],['anatomical-point','boundary-curve','boundary-curve'],['Free-margin central nodule sites; nodule volumes are not constructed.','Coaptation-margin lunule sites; microscopic bands remain unconstructed.','Peripheral leaflet attachment/commissural reference sites.']):add(f'{base}.{number:03d}',valve[0],points[0],kind,name,points=points)
# Registered source locations for the conducting system remain schematic;
# the tissue network is not added to or carved out of the myocardium.
conduction=[('A12.1.06.003',18229,ra),('A12.1.06.004',17776,ias),('A12.1.06.005',22137,ivs),('A12.1.06.006',20044,ivs),('A12.1.06.007',18445,ivs),('A12.1.06.008',23412,endo)]
for term,cid,carrier in conduction:sourced(term,cid,carrier,'anatomical-point'if '003'in term or '004'in term else'boundary-curve')
points=[spec[t]['point']for t,cid,carrier in conduction];add('A12.1.06.002',ivs,points[1],'collective','Registered schematic nodal/bundle locations; no conduction-tissue volumes are claimed.',[cid for t,cid,carrier in conduction],points=points)
add('A12.1.07.001',epi,(pv+av+c('left-atrium'))/3,'anatomical-point','Transverse pericardial recess reference behind the arterial roots and anterior to the atrial/caval structures.')
add('A12.1.07.002',epi,c('left-atrium')+np.array([0,0,-.009]),'anatomical-point','Oblique pericardial recess reference posterior to the left atrium between pulmonary-venous reflections.')
for term,cid in [('A12.1.08.002',23932),('A12.1.08.006',23251)]:
 q=landmark(cid);q[2]=ref[cid]['bounds_m'][1][2];q=nearby(epi,q);add(term,epi,q,'regional-reference','Named pericardial covering reference location only; a separate finite tissue layer is not constructed.',[cid])
add('A12.1.08.001',epi,spec['A12.1.08.002']['point'],'collective','Pericardial complex reference; only its visceral/epicardial component currently has a finite tissue coat.',[23932,23251,20879])
add('A12.1.08.005',epi,spec['A12.1.08.006']['point'],'collective','Serous-pericardial complex reference; the separate parietal wall remains unconstructed.',[23251,20879])
front=epicardial_features['A12.1.00.003'];add('A12.1.08.003',epi,front+np.array([0,0,.002]),'boundary-curve','Sternopericardial attachment direction reference, not a constructed ligament.',points=[front,front+np.array([0,0,.006])])
# Complete the catalogue now with explicitly provisional sites. These entries
# do not install the failed serosal prototype, create holes, or imply that a
# point is the final representation of a surface, band or tissue volume.
add('A12.1.00.008',epi,nearby(epi,apex),'regional-reference','Apical-notch neighborhood at the current ventricular apex.')
add('A12.1.00.024',lv,nearby(lv,.5*center(marker['papillaryPartIds'])+.5*c('interventricular-septum')),'regional-reference','Left-ventricular false-tendon neighborhood; no particular variable band is asserted.')
add('A12.1.00.028',rv,nearby(rv,.65*pv+.35*tv),'regional-reference','Infundibular fibrous-attachment neighborhood between the right AV and pulmonary apparatus.')
add('A12.1.00.029',ra,nearby(ra,.5*fossa+.5*ivc),'boundary-curve','Todaro-tendon reference between the inferior caval valve and septal AV region.',points=[nearby(ra,ivc),nearby(ra,.5*fossa+.5*tv)])
koch=[nearby(ra,cs),nearby(ra,tv),nearby(ra,.5*fossa+.5*ivc)];add('A12.1.00.030',ra,np.mean(koch,axis=0),'regional-reference','Koch/AV-node triangle neighborhood, bounded provisionally by coronary-sinus, septal tricuspid and Todaro reference sites.',points=koch)
add('A12.1.01.004',ra,nearby(ra,.6*c('right-atrium')+.4*ivc),'regional-reference','Right-atrial endocardial neighborhood for smallest cardiac-vein openings; no individual opening is localized.')
add('A12.1.01.007',ias,fossa,'regional-reference','Foramen-ovale developmental site at the interatrial septum; the default adult septum remains closed.')
add('A12.1.01.014',ra,nearby(ra,.6*fossa+.4*svc),'regional-reference','Superior interatrial-fold neighborhood above the fossa ovalis.')
add('A12.1.04.006',ids('mitral-valve')[0],nearby(ids('mitral-valve')[0],mv),'regional-reference','Mitral commissural-leaflet neighborhood; separate accessory cusp anatomy is provisional.')
add('A12.1.08.004',epi,nearby(epi,c('left-atrium')+np.array([0,.012,-.008])),'regional-reference','Dorsal pericardial neighborhood below the bronchial bifurcation; membrane extent remains unbuilt.')
add('A12.1.08.008',epi,nearby(epi,front),'regional-reference','Serosal covering reference on the existing epicardial carrier; no new tissue stratum is installed.')
add('A12.1.08.009',epi,nearby(epi,front+np.array([.004,0,0])),'regional-reference','Subserosal covering reference on the existing epicardial carrier; no new tissue stratum is installed.')
add('A12.1.08.010',epi,nearby(epi,.5*landmark(23809)+.5*max(mouths,key=lambda q:q[1])),'regional-reference','Left caval/Marshall-fold neighborhood between the left auricle and superior pulmonary-venous reflection.')
assert set(spec)==set(new_terms),(sorted(set(new_terms)-set(spec)),sorted(set(spec)-set(new_terms)))
ep=parts[epi];R,scale=instance_pose(ep);origin=np.array(ep['geometryOrigin'])
local=lambda q:(np.asarray(q)-origin)@R/scale
records=[]
for term,s in sorted(spec.items()):
 record={k:v for k,v in s.items()if k not in ['point','points']};record.update(term=term,name=authority[term]['name_en'],authorityParent=authority[term]['parent_id'],point_local=local(s['point']).tolist(),points_local=None if s['points']is None else local(s['points']).tolist(),tissueGeometryEstablished=False,provisional=True);records.append(record)
atomic_json(out/'reference-plan.json',dict(schema='cardiac-label-reference-plan/v1',authority_terms=len(authority),geometry_backed_term_controls=len(existing),qualified_reference_controls=len(records),labeled_terms=len(existing)+len(records),reference_statistics=str(reference_path),reference_statistics_sha256=digest_file(reference_path),registration_rms_metres=reference['fit']['landmark_rms_m'],unfinished=unfinished,sites=records,additional_meshes=0,additional_textures=0))
for sex,a in atlases.items():
 a=copy.deepcopy(a);pmap={p['id']:p for p in a['parts']};m=a['finiteCardiacIntegration'];owned=set(sum([m.get(k,[])for k in ['wallPartIds','valvePartIds','papillaryPartIds','chordalPartIds','rootPartIds']],[]));ep=pmap[f'TA98HEART:{sex}:epicardium'];rot,ss=instance_pose(ep);oo=np.array(ep['geometryOrigin']);forward=lambda q:(np.asarray(q)@rot.T)*ss+oo;site_rows=[];concepts={c['id']:c for c in a['concepts']};existing_terms={g['term']for g in m['termGroups']if g['kind']=='exact'};assert existing_terms==existing
 key_by_term={g['term']:g['key']for g in m['termGroups']if g['kind']=='exact'}
 for row in records:
  term=row['term'];cid=f'TA98HEART:{sex}:label-site:{term}';carrier=row['carrier'].replace(':male:',f':{sex}:')
  if carrier not in owned:
   # Native ZA papillary IDs are preserved as donor-qualified identifiers in
   # the female transfer. Resolve the actual term-group source, never invent
   # an absent mirrored ID.
   matches=[g for g in marker['termGroups']if row['carrier']in g['partIds']];assert matches;original_group=min(matches,key=lambda g:len(g['partIds']));carrier=next(g for g in m['termGroups']if g['key']==original_group['key'])['partIds'][0]
  assert carrier in owned and carrier in pmap and not pmap[carrier].get('suppressed'),(sex,term,carrier)
  collective={'A12.1.00.018':[f'TA98HEART:{sex}:myocardium:right-atrium',f'TA98HEART:{sex}:myocardium:left-atrium'],'A12.1.00.020':[f'TA98HEART:{sex}:endocardium',f'TA98HEART:{sex}:myocardium:right-ventricle',f'TA98HEART:{sex}:myocardium:left-ventricle'],'A12.1.00.022':m['papillaryPartIds'],'A12.1.06.002':[f'TA98HEART:{sex}:myocardium:interventricular-muscular-septum',f'TA98HEART:{sex}:myocardium:interatrial-septum',f'TA98HEART:{sex}:myocardium:right-atrium',f'TA98HEART:{sex}:endocardium']};carriers=collective.get(term,[carrier]);assert carrier in carriers
  title=row['name'][0].upper()+row['name'][1:]
  if term in ['A12.1.02.014','A12.1.02.015','A12.1.02.016']:title+=' of pulmonary valve'
  if term in ['A12.1.04.016','A12.1.04.017','A12.1.04.018']:title+=' of aortic valve'
  if term in ['A12.1.02.021','A12.1.04.011']:title+=' of '+('right'if '.02.'in term else'left')+' ventricle'
  if term in ['A12.1.01.008','A12.1.03.003']:title='Pectinate muscles of '+('right'if '.01.'in term else'left')+' atrium'
  if term in ['A12.1.04.007','A12.1.04.008']:title+=' of left ventricle'
  if term=='A12.1.00.030':title='Triangle of atrioventricular node (Koch)'
  parent=key_by_term.get(row['authorityParent'],'heart');key='label-site:'+term;point=forward(row['point_local']).tolist();points=None if row['points_local']is None else forward(row['points_local']).tolist();coverage='Provisional '+row['representation'].replace('-',' ')+' representation. '+row['derivation']+' '+unfinished.get(term,'');concept=dict(id=cid,name=title,elements=carriers,ta98Term=term,ta98Kind='exact',partOf=m['wholeConceptId'],coverageLimitation=coverage,cardiacReferenceSite=True,distinctTissueGeometryEstablished=False);assert cid not in concepts;a['concepts'].append(concept);site_rows.append(dict(term=term,conceptId=cid,key=key,parent=parent,name=concept['name'],partIds=carriers,kind='exact',referenceOnly=True,representation=row['representation'],anchor=dict(partId=carrier,point=point),points=points,referenceCids=row['referenceCids'],tissueGeometryEstablished=False));key_by_term[term]=key
 # Resolve parents after all sites have been created (e.g. the conducting
 # system and pericardium are themselves new collective references).
 byterm={g['term']:g['key']for g in m['termGroups']if g['kind']=='exact'}|{g['term']:g['key']for g in site_rows}
 for row in site_rows:
  row['parent']=byterm.get(authority[row['term']]['parent_id'],'heart');parent_cid=next((g['conceptId']for g in [*m['termGroups'],*site_rows]if g['key']==row['parent']),m['wholeConceptId']);next(c for c in a['concepts']if c['id']==row['conceptId'])['partOf']=parent_cid
 poses={pid:dict(geometryOrigin=pmap[pid].get('geometryOrigin',[0,0,0]),geometryInstance=pmap[pid].get('geometryInstance'),chunk=pmap[pid]['chunk'],positions=pmap[pid]['positions'],indices=pmap[pid]['indices'],indexCount=pmap[pid]['indexCount'])for pid in set(g['anchor']['partId']for g in site_rows)}
 a['cardiacLabelIntegration']=dict(schema='cardiac-qualified-label-sites/v1',sex=sex,bankSha256=m['coupledBankSha256'],referencePlanSha256=digest_file(out/'reference-plan.json'),termGroups=site_rows,carrierPoses=poses,geometryBackedTermCount=len(existing),referenceTermCount=len(site_rows),labeledTermCount=len(existing)+len(site_rows),distinctTissueCoverageUnchanged=True)
 target=out/sex;target.mkdir();atomic_json(target/'atlas.json',a);atomic_json(target/'cardiac-label-sites.json',dict(schema='cardiac-qualified-label-bundle/v1',integration=a['cardiacLabelIntegration'],concepts=[c for c in a['concepts']if c['id']in{g['conceptId']for g in site_rows}]));atomic_json(target/'label-proof.json',dict(source_catalogue_sha256=digest_file(paths[sex]),candidate_catalogue_sha256=digest_file(target/'atlas.json'),all_parts_literal=a['parts']==atlases[sex]['parts'],all_chunks_literal=a['chunks']==atlases[sex]['chunks'],fixed_skeleton=True,qualified_reference_terms=len(site_rows),geometry_backed_terms=len(existing),labeled_terms=len(existing)+len(site_rows),all_anchors_finite=all(np.isfinite(g['anchor']['point']).all()for g in site_rows),new_meshes=0,new_textures=0,anatomical_acceptance=False))
print(f'Prepared {len(authority)} labeled cardiac terms per model: {len(existing)} geometry-backed controls and {len(records)} explicitly provisional references; no geometry changed. Elapsed {time.monotonic()-started:.2f}s.',flush=True)
