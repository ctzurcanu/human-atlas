import test from 'node:test';import assert from 'node:assert/strict';import {spawnSync} from 'node:child_process';
test('Eye source correction accepts rebuilt broad placement and is idempotent, while rejecting conflicting anatomy',()=>{
const script=String.raw`import importlib.util
s=importlib.util.spec_from_file_location('eye','scripts/stage-ta98-female-eye-identities.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def atlas(kind):
 a={'parts':[],'concepts':[]}
 for side in ['r','l']:
  for source,term in [('iris','A15.2.03.020'),('pupil','A15.2.03.028')]:
   cid='LOCAL:female:'+side+'_'+source;pid=cid+':1';a['parts'].append({'id':pid,'name':source,'bounds':[[0,0,0],[1,1,1]],'positions':123});c={'id':cid,'name':source,'elements':[pid]}
   if kind!='unmapped':c.update(ta98Term=term if kind=='exact' else m.WHOLE,ta98Kind='exact' if kind=='exact' else 'parent')
   a['concepts'].append(c)
 return a
for kind in ['exact','parent','unmapped']:
 original=atlas(kind);a,changes=m.corrected(original);m.verify(original,a,changes);again,againChanges=m.corrected(a);assert again==a;m.verify(a,again,againChanges)
bad=atlas('exact');bad['concepts'][0]['ta98Term']='A15.2.02.001'
try:m.corrected(bad)
except AssertionError:pass
else:raise AssertionError('Conflicting anatomy accepted')
print('Rebuild and idempotence checks passed')`;
const result=spawnSync('python3',['-c',script],{encoding:'utf8'});assert.equal(result.status,0,result.stderr);assert.match(result.stdout,/checks passed/);
});
