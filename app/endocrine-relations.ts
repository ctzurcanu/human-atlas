import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

function combinedParathyroids(p:Part):boolean{
 const sex=/^LOCAL:(male|female):parathyroid_glands:\d+$/.exec(p.id)?.[1];
 return !!sex&&!p.suppressed&&p.system==='endocrine'&&p.conceptId===`LOCAL:${sex}:parathyroid_glands`&&p.sourceId==='parathyroid_glands'&&p.name==='Parathyroid glands';
}

export function endocrineDescription(p:Part):{text:string;url:string}|undefined{
 const type=pituitaryKind(p);
 if(type==='anterior')return {text:'The adenohypophysis is the glandular part of the pituitary. Hypothalamic releasing and inhibiting hormones regulate its secretion through the hypophysial portal circulation.',url:'https://www.ncbi.nlm.nih.gov/books/NBK279126/'};
 if(type==='posterior')return {text:'The neurohypophysis is the neural part of the pituitary. Axons from hypothalamic neurons deliver vasopressin and oxytocin for release into the circulation.',url:'https://www.ncbi.nlm.nih.gov/books/NBK279126/'};
 if(type==='whole')return {text:'The pituitary gland includes glandular and neural regions with different connections to the hypothalamus: portal blood vessels for anterior-lobe regulation and axons for posterior-lobe neurosecretion.',url:'https://www.ncbi.nlm.nih.gov/books/NBK279126/'};
 if(combinedParathyroids(p))return {text:'The parathyroid glands are separate endocrine glands, usually near the posterior surfaces of the thyroid lobes. Their number and location vary; this combined source does not identify each superior or inferior gland individually.',url:'https://www.ncbi.nlm.nih.gov/books/NBK244/'};
}

function pituitaryKind(p:Part):'whole'|'anterior'|'posterior'|undefined{
 if(p.suppressed||p.system!=='endocrine')return;
 if(p.id==='LOCAL:female:pituitary_gland:64'&&p.conceptId==='LOCAL:female:pituitary_gland'&&p.sourceId==='pituitary_gland'&&p.name==='Pituitary gland')return 'whole';
 if(p.id==='ZA:Adenohypophysis'&&p.conceptId===p.id&&p.sourceId==='Adenohypophysis'&&p.name==='Adenohypophysis')return 'anterior';
 if(p.id==='ZA:Neurohypophysis'&&p.conceptId===p.id&&p.sourceId==='Neurohypophysis'&&p.name==='Neurohypophysis')return 'posterior';
}
export function endocrineSourceLimitation(p:Part):string|undefined{
 const type=pituitaryKind(p);
 if(type)return `${type==='whole'?'Whole native pituitary source; separately selectable anterior/posterior lobes are missing.':'Named native '+type+' pituitary source.'} Cell populations, portal capillaries, axonal terminals, stalk interfaces and tissue extent remain unverified. A hormonal association is not a verified mesh connection.`;
 if(combinedParathyroids(p))return 'Combined named parathyroid source. Separate superior/inferior and left/right gland identities, capsule interfaces and microscopic cell populations are not established by the combined mesh; gland number and location vary.';
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Endotext NBK279126 distinguishes portal control of anterior pituitary
 * from hypothalamic axonal projections to posterior pituitary. NBK244:
 * parathyroids are distinct glands, usually posterolateral to thyroid lobes. */
export function endocrineRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>();
 const add=(a:Part,b:Part,kind:'after'|'adjacent',note:string,via?:string[])=>{
  for(const [from,to,k] of [[a,b,kind],[b,a,kind==='after'?'before':'adjacent']] as const){const links=graph.get(from.id)??[];links.push({kind:k,target:to,note,...(via?{via}: {})});graph.set(from.id,links);}
 };
 const hypothalami=parts.filter(p=>!p.suppressed&&p.system==='nervous'&&(p.id==='ZA:Hypothalamus'&&p.conceptId===p.id&&p.sourceId==='Hypothalamus'&&p.name==='Hypothalamus'||/^LOCAL:female:(?:l_hypothalamus:29|r_hypothalamus:15)$/.test(p.id)&&p.conceptId==='LOCAL:female:'+p.sourceId&&/^Hypothalamus \((left|right)\)$/.test(p.name)));
 for(const gland of parts){
  const type=pituitaryKind(gland);if(!type)continue;
  for(const hypothalamus of hypothalami.filter(p=>nativePartFrame(p)===nativePartFrame(gland))){
   const note=type==='anterior'?'Hypothalamic releasing and inhibiting hormones reach the adenohypophysis through the median eminence and hypophysial portal vessels. This hormonal route is distinct from the axonal route to the neurohypophysis. Portal capillary beds and source stalk/lumen interfaces remain unverified.':type==='posterior'?'Hypothalamic neurons project axons to the neurohypophysis, where vasopressin and oxytocin are released into blood. This neurosecretory route is distinct from anterior-lobe portal signaling; individual tracts, terminals and source stalk continuity remain unverified.':'Hypothalamic regulation reaches the anterior pituitary through portal vessels and the posterior pituitary through axonal projections. This whole pituitary source does not separately represent those lobes, vessels or tracts; the stalk interface remains unverified.';
   add(hypothalamus,gland,'after',note,type==='anterior'?['Median eminence','Hypophysial portal vessels']:type==='posterior'?['Hypothalamohypophysial tract']:['Hypophysial stalk (infundibulum)']);
  }
 }
 const thyroid=parts.find(p=>!p.suppressed&&p.system==='endocrine'&&p.id==='LOCAL:female:thyroid_gland:1081'&&p.conceptId==='LOCAL:female:thyroid_gland'&&p.sourceId==='thyroid_gland'&&p.name==='Thyroid gland');
 const para=parts.find(p=>!p.suppressed&&p.system==='endocrine'&&p.id==='LOCAL:female:parathyroid_glands:1083'&&p.conceptId==='LOCAL:female:parathyroid_glands'&&p.sourceId==='parathyroid_glands'&&p.name==='Parathyroid glands');
 if(thyroid&&para)add(thyroid,para,'adjacent','Parathyroids are distinct endocrine glands usually near the posterior thyroid lobes; number and location vary. They are not subdivisions of thyroid tissue. This combined native parathyroid source does not establish individual gland positions, shared capsules or verified packing.');
 cache.set(parts,graph);return graph;
}
