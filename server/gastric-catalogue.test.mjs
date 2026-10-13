import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const read=path=>readFile(path,'utf8').then(JSON.parse);
test('Installed female finite tissues resolve named TA98 unions without duplicate wall coats',async t=>{
 let atlas;try{atlas=await read('.local-models/ta98-female-runtime.json');}catch(e){if(e.code==='ENOENT'){t.skip('Local review assets unavailable');return;}throw e;}
 if(!atlas.parts.some(p=>p.pairedTissue)){t.skip('Finite wall not installed');return;}
 const parts=atlas.parts.filter(p=>p.pairedTissue&&p.id.startsWith('GASTRIC:female:'));
 assert.equal(parts.length,33);assert.equal(new Set(parts.map(p=>`${p.chunk}:${p.positions}`)).size,1,'one immutable position bank after runtime repack');assert.equal(new Set(parts.map(p=>`${p.chunk}:${p.normals}`)).size,1);
 const ids=new Set(parts.map(p=>p.id));const concepts=new Map(atlas.concepts.map(c=>[c.id,c]));const layer=c=>concepts.get('GASTRIC:female:'+c);
 for(const c of ['019','020','021','022','023','024','025','026','027','029']){assert.ok(layer(c)?.elements.length);assert.equal(layer(c).ta98Term,'A05.5.01.'+c);assert.ok(layer(c).elements.every(id=>ids.has(id)));}
 assert.deepEqual(new Set(layer('021').elements),new Set(['022','023','025'].flatMap(c=>layer(c).elements)));assert.ok(layer('023').elements.includes(layer('024').elements[0]));assert.equal(layer('024').elements.length,1,'sphincter is one local thickening');assert.equal(layer('024').partOf,layer('023').id);
 assert.ok(layer('027').elements.every(id=>['027','029'].some(c=>id.startsWith('GASTRIC:female:'+c+':'))),'mucosa parent selects real tissue union');
 for(const p of parts){assert.ok(p.geometryOrigin&&p.normalSigns!==undefined);assert.equal(p.tissueDrawGroups.reduce((n,g)=>n+g.count,0),p.indexCount);for(const g of p.tissueDrawGroups)if(g.neighbor){assert.ok(ids.has(g.neighbor));assert.ok(parts.find(q=>q.id===g.neighbor).tissueDrawGroups.some(q=>q.neighbor===p.id),'joined boundary has reciprocal owner');}}
 for(const old of atlas.parts)assert.ok(!old.id.startsWith('TA98FEMALE:REGION:A05.5.01.'),'original coarse wall is not drawn over layered wall');
 // Check the actual outer cardia contour exposed by hiding the fundus.
 // A closed shell alone cannot detect the previous angular source partition.
 const cardia=parts.find(p=>p.id==='GASTRIC:female:019:007');
 const cardiaBytes=await readFile('.local-models/'+atlas.chunks[cardia.chunk].url.replace(/^\/local-models\//,''));
 const sheetCount=cardia.vertexCount/9,edges=new Map();
 for(let k=0;k<cardia.indexCount;k+=3){
  const tri=[0,1,2].map(j=>cardiaBytes.readUInt32LE(cardia.indices+(k+j)*4));
  if(tri.some(i=>i>=sheetCount))continue;
  for(let j=0;j<3;j++){const edge=[tri[j],tri[(j+1)%3]].sort((a,b)=>a-b),key=edge.join(':');edges.set(key,{edge,count:(edges.get(key)?.count??0)+1});}
 }
 const adjacency=new Map();
 for(const {edge:[a,b],count} of edges.values())if(count===1){for(const [i,j] of [[a,b],[b,a]])adjacency.set(i,[...(adjacency.get(i)??[]),j]);}
 assert.ok([...adjacency.values()].every(n=>n.length===2),'opening and collar remain closed curves');
 const remaining=new Set(adjacency.keys()),loops=[];
 while(remaining.size){const stack=[remaining.values().next().value],loop=[];remaining.delete(stack[0]);while(stack.length){const i=stack.pop();loop.push(i);for(const j of adjacency.get(i))if(remaining.delete(j))stack.push(j);}loops.push(loop);}
 assert.equal(loops.length,2,'cardia has its opening and one continuous adjoining collar');
 const collar=loops.sort((a,b)=>b.length-a.length)[0];
 const position=i=>[0,1,2].map(j=>cardiaBytes.readFloatLE(cardia.positions+(i*3+j)*4));
 for(const i of collar){const p=position(i),[x,y]=adjacency.get(i).map(j=>position(j).map((v,k)=>v-p[k]));const cosine=x.reduce((sum,v,k)=>sum+v*y[k],0)/Math.hypot(...x)/Math.hypot(...y);assert.ok(Math.acos(Math.max(-1,Math.min(1,cosine)))*180/Math.PI>140,'cardia collar must not regain acute corners');}
 const parents=(await read('app/data/ta98-parents.json')).parents;assert.equal(parents['A05.5.01.024'],'A05.5.01.023');assert.equal(parents['A05.5.01.025'],'A05.5.01.021');
 // A generated fold field may decorate only the innermost mucosal sheet.
 // It must not move the geometry that the cap counts or leak onto the exterior.
 assert.ok(atlas.gastricWallRelief,'runtime retains the installed generated fold field');
 {
  const mucosal=parts.filter(p=>p.layerTerm==='A05.5.01.027');assert.equal(mucosal.length,4);
  for(const p of mucosal){
   const url=atlas.chunks[p.chunk].url;const raw=await readFile('.local-models/'+url.replace(/^\/local-models\//,''));
   const mask=new Uint8Array(raw.buffer,raw.byteOffset+p.gastricLumen,p.vertexCount);
   const sheet=p.vertexCount/9;assert.ok(Number.isInteger(sheet));
   assert.ok(mask.subarray(0,8*sheet).every(v=>v===0),'outer tissues have no fold mask');
   assert.ok(mask.subarray(8*sheet).every(v=>v===255),'the whole luminal interface receives folds');
   const material=atlas.materials[p.material];assert.ok(material.gastricRelief?.scale>0);assert.ok(new URL(material.gastricRelief.heightMap,'http://localhost').pathname.endsWith('mucosal-fold-height.png'));
  }
  assert.ok(parts.filter(p=>p.layerTerm!=='A05.5.01.027').every(p=>p.gastricLumen===undefined),'other layers remain independent of mucosal relief');
  const serosa=parts.find(p=>p.layerTerm==='A05.5.01.019');assert.equal(serosa.tissueAppearance.pattern,'fibrous');
  assert.ok(serosa.tissueAppearance.color[0]>serosa.tissueAppearance.color[1]*1.5,'serosa remains pink');
  assert.equal(new Set(parts.filter(p=>p.layerTerm!=='A05.5.01.024').map(p=>p.tissueAppearance.color.join(','))).size,8,'eight distinguishable tissue colors');
 }
});
