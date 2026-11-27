import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';

const bytes=readFileSync(new URL('../public/assets/ta98.json',import.meta.url));
const hierarchy=JSON.parse(bytes),nodes=new Map(hierarchy.nodes.map(n=>[n.id,n]));
const index=JSON.parse(readFileSync(new URL('../app/data/ta98-parents.json',import.meta.url)));

test('Female and male opening states have distinct, sex-qualified canonical parents',()=>{
 for(const [parent,children]of [
  ['A08.3.01.028F',['A09.2.03.001','A09.2.03.002']],
  ['A08.3.01.027M',['A09.4.02.001','A09.4.02.002']],
 ])for(const code of children){
  assert.equal(index.parents[code],parent);
  assert.ok(nodes.get(parent).childIds.includes(code));
  const wrong=parent.endsWith('F')?'A08.3.01.027M':'A08.3.01.028F';
  assert.ok(!nodes.get(wrong).childIds.includes(code));
 }
});

test('Runtime ancestry is bound to the exact published hierarchy',()=>{
 assert.equal(index.hierarchySha256,createHash('sha256').update(bytes).digest('hex'));
});

test('Rectovesicalis appears in bladder context with its canonical muscle ancestry intact',()=>{
 const occurrence=nodes.get('TA98CONTEXT:A08.3.01.001:A04.5.03.017');
 assert.equal(occurrence.matchFrom,'A04.5.03.017');
 assert.ok(nodes.get('A08.3.01.001').childIds.includes(occurrence.id));
 assert.equal(index.parents['A04.5.03.017'],'A04.5.03.016');
 assert.ok(nodes.get('A04.5.03.016').childIds.includes('A04.5.03.017'));
 assert.ok(!nodes.get('A08.3.01.001').childIds.includes('A04.5.03.017'));
});
