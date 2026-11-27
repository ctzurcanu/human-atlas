import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
import {Quaternion,Vector3} from 'three';

const built=await build({stdin:{contents:"export * from './app/heart-wall';export * from './app/heart-system';export * from './app/organ-frame';export * from './app/guest-hierarchy';export * from './app/cardiac-relations';export * from './app/ta98-modeled-parent';export * from './app/cardiac-surface-features';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const authority=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
const read=sex=>JSON.parse(fs.readFileSync(`${process.env.ATLAS_HEART_CANDIDATE??'.local-models/ta98-review/bowel-wall-local-atlas'}/${sex}/atlas.json`));

for(const sex of ['male','female']){
 test(`${sex}: installed heart exposes real tissues in TA98 and Systems`,()=>{
  const a=read(sex),m=api.cardiacIntegration(a);assert.ok(m);
  const groups=api.cardiacTermParts(a);assert.equal(groups.size,42);
  for(const [side,term] of [['right','A12.1.00.027'],['left','A12.1.00.032']]){
   const id=`TA98HEART:${sex}:fibrous:${side}-annulus`,p=a.parts.find(p=>p.id===id);
   assert.ok(p?.pairedTissue&&p.sectionSolid&&p.cardiacAnnularPartition);
   assert.deepEqual(groups.get(term).map(q=>q.id),[id]);
   const contacts=m.coupledAttachmentSemantics.filter(r=>[r.first,r.second].includes(id));
   assert.ok(contacts.some(r=>r.kind==='annular-leaflet-insertion'&&r.fusion_inferred));
   assert.ok(contacts.some(r=>r.kind==='annular-myocardial-insertion'&&r.fusion_inferred));
   for(const g of p.tissueDrawGroups.filter(g=>g.neighbor)){
    const peer=a.parts.find(q=>q.id===g.neighbor);
    assert.ok(peer.tissueDrawGroups.some(q=>q.neighbor===id));
   }
  }
  if(m.chordalPartIds?.length){
   assert.equal(m.chordalPartIds.length,5);
   assert.deepEqual(new Set(groups.get('A12.1.00.023').map(p=>p.id)),new Set(m.chordalPartIds));
   const parent=a.concepts.find(c=>c.id===`TA98HEART:${sex}:chordae`);assert.equal(parent.partOf,m.wholeConceptId);
   for(const id of m.chordalPartIds){const p=a.parts.find(p=>p.id===id);assert.ok(p.sectionSolid&&p.pairedTissue&&p.normalSigns!==undefined);assert.ok(p.tissueDrawGroups.some(g=>m.valvePartIds.includes(g.neighbor)));assert.equal(a.concepts.find(c=>c.id===p.conceptId).partOf,parent.id);}
  }
  const queue=[...api.resolveGuestHierarchy(a,authority)],terms=new Map();
  while(queue.length){const n=queue.pop();if(terms.has(n.id))continue;terms.set(n.id,n);queue.push(...n.children);}
  for(const [term,parts] of groups)assert.ok(parts.every(p=>terms.get(term)?.directParts.some(q=>q.id===p.id)),term);
  const hearts=api.withHeartSystem(a,[]);assert.equal(hearts.length,1);assert.equal(hearts[0].name,'Heart');
  const owned=new Set([...m.wallPartIds,...m.valvePartIds,...m.papillaryPartIds,...m.chordalPartIds??[],...m.rootPartIds]);
  assert.deepEqual(new Set(hearts[0].parts.map(p=>p.id)),owned);
  assert.ok(hearts[0].nodes.some(n=>n.name==='Vessel roots'));
  assert.ok([...owned].every(id=>!a.parts.find(p=>p.id===id).suppressed));
  for(const p of a.parts.filter(p=>p.suppressed&&p.replacementIds?.some(id=>owned.has(id))))assert.ok(!hearts[0].parts.some(q=>q.id===p.id));
 });
 test(`${sex}: peeling accounts for every selected heart part once`,()=>{
  const a=read(sex),m=api.cardiacIntegration(a),whole=a.concepts.find(c=>c.id===m.wholeConceptId),steps=api.cardiacWallPeelSteps(a,whole);
  assert.equal(steps.length,5);assert.deepEqual(steps.slice(0,3).map(s=>s.ids.length),[1,10,1]);
  assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(whole.elements));
  assert.equal(steps.flatMap(s=>s.ids).length,whole.elements.length);
  assert.equal(api.cardiacWallPeelSteps(a,{...whole,id:'selection-set'}),null);
  const broken=structuredClone(a);broken.parts.find(p=>p.id===m.wallPartIds[1]).geometryOrigin=[2,3,4];
  assert.equal(api.cardiacWallPeelSteps(broken,whole),null);
 });
 test(`${sex}: chambers retain their own valve and venous inflow relationships`,()=>{
  const a=read(sex),graph=api.cardiacRelations(a.parts),m=a.finiteCardiacIntegration;
  for(const g of m.termGroups.filter(g=>['mitral-valve','tricuspid-valve','aortic-valve','pulmonary-valve'].includes(g.key)))for(const id of g.partIds)assert.ok(graph.get(id)?.some(r=>r.kind==='adjacent'&&r.target.id.startsWith(`TA98HEART:${sex}:myocardium:`)),g.name);
  for(const id of m.rootPartIds){const p=a.parts.find(p=>p.id===id);if(/pulmonary vein|vena cava/i.test(p.name))assert.ok(graph.get(id)?.some(r=>r.kind==='after'&&/atrium$/.test(r.target.name.toLowerCase())),p.name);}
  for(const id of m.chordalPartIds??[]){const p=a.parts.find(p=>p.id===id),links=graph.get(id)??[],c=p.cardiacChordalPartition;assert.ok(links.some(r=>r.kind==='continuous'&&r.target.id===c.leafletPartId));assert.ok(links.some(r=>r.kind===(c.papillaryInterfaceVerified?'continuous':'adjacent')&&r.target.id===c.papillaryPartId));assert.ok(!links.some(r=>r.kind===(c.papillaryInterfaceVerified?'adjacent':'continuous')&&r.target.id===c.papillaryPartId));}
  for(const contact of m.apparatusInterfaceSemantics??[]){const links=graph.get(contact.first)??[],kind=contact.kind.startsWith('native-')?'continuous':'adjacent';assert.ok(links.some(r=>r.kind===kind&&r.target.id===contact.second),contact.kind);if(kind==='adjacent')assert.ok(!links.some(r=>r.kind==='continuous'&&r.target.id===contact.second));}
  for(const contact of m.coupledAttachmentSemantics??[])if(['papillary-wall-insertion','leaflet-host-attachment'].includes(contact.kind))assert.ok(graph.get(contact.first)?.some(r=>r.kind==='continuous'&&r.target.id===contact.second),contact.kind);
 });
 test(`${sex}: named exterior features bind real stored facets and labels`,()=>{
  const a=read(sex);if(!a.cardiacFeatureIntegration)return;
  const marker=api.cardiacFeatureRegistry(a);assert.ok(marker);assert.equal(marker.termGroups.length,6);
  const p=a.parts.find(p=>p.id===marker.hostPartId),url=a.chunks[p.chunk].url;
  const bank=fs.readFileSync(url.startsWith('/local-models/')?'.local-models/'+url.slice('/local-models/'.length):'public'+url);
  for(const g of marker.termGroups){
   const choice=a.concepts.find(c=>c.id===g.conceptId);assert.equal(api.selectedCardiacSurfaceFeature(a,choice).conceptId,g.conceptId);
   const actual=new Vector3();for(let k=0;k<3;k++){const i=bank.readUInt32LE(p.indices+(g.anchor.face*3+k)*4);const local=new Vector3(...[0,1,2].map(axis=>bank.readFloatLE(p.positions+i*12+axis*4)));actual.addScaledVector(api.organAnatomicalPoint(p,local),g.anchor.barycentric[k]);}
   assert.ok(actual.distanceTo(new Vector3(...g.anchor.point))<1e-12,g.term);assert.deepEqual(api.cardiacSurfaceLabelReference(a,choice).point,g.anchor.point);
  }
  for(const mutate of [b=>b.cardiacFeatureIntegration.bankSha256='0'.repeat(64),b=>b.cardiacFeatureIntegration.termGroups[0].surface.faces=[-1]]){const broken=structuredClone(a);mutate(broken);assert.equal(api.cardiacFeatureRegistry(broken),null);}
 });
}

test('placed female geometry, labels and exploded poses use the same positive similarity',()=>{
 const a=read('female'),m=a.finiteCardiacIntegration,p=a.parts.find(p=>p.id===m.wallPartIds[0]);assert.ok(p.geometryInstance);
 const local=new Vector3(.013,-.008,.027),world=api.organAnatomicalPoint(p,local);
 assert.ok(api.organGeometryPoint(p,world).distanceTo(local)<1e-12);
 const assembled=api.organPose(p,new Vector3());assert.ok(local.clone().applyMatrix4(assembled).distanceTo(world)<1e-12);
 const box=api.organBounds(p),posed=api.organBounds(p,assembled);assert.ok(box.min.distanceTo(posed.min)<1e-12);assert.ok(box.max.distanceTo(posed.max)<1e-12);
 const turn=new Quaternion().setFromAxisAngle(new Vector3(0,0,1),.42),offset=new Vector3(.08,.01,-.02),pivot=new Vector3(0,1.1,0);
 const expected=world.clone().add(offset).sub(pivot).applyQuaternion(turn).add(pivot);
 assert.ok(local.clone().applyMatrix4(api.organPose(p,offset,turn,pivot)).distanceTo(expected)<1e-12);
});

test('posterior source identity needs its source-qualified evidence',()=>{
 const p={id:'ZA:Inferior leaflet of right atrioventricular valve',conceptId:'ZA:Inferior leaflet of right atrioventricular valve',name:'Inferior leaflet of right atrioventricular valve',system:'cardiac'};
 assert.match(api.cardiacSourceLimitation(p),/correspondence.*unresolved/);
 assert.match(api.cardiacSourceLimitation({...p,ta98Term:'A12.1.02.005'}),/correspondence.*unresolved/);
 assert.ok(!/correspondence.*unresolved/.test(api.cardiacSourceLimitation({...p,ta98Term:'A12.1.02.005',cardiacIdentityEvidence:{named_source:'FJ2433'}})??''));
});

test('the recovered papillary head needs its declared native registration',()=>{
 const p=read('male').parts.find(p=>p.id==='TA98CARDIAC:FJ2418');assert.equal(api.nativePartFrame(p),'ZA:');
 for(const change of [{registration:undefined},{sourceId:'another donor'},{reviewSource:'female'},{provenance:{...p.provenance,url:'https://example.com/unregistered'}},{registration:{...p.registration,pairs:[]}}, {registration:{...p.registration,maxLandmarkResidualMeters:-1}}])assert.equal(api.nativePartFrame({...p,...change}),undefined);
});
