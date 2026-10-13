import type {Part} from './anatomy';
import type {RelationKind,ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const normalized=(p:Part)=>p.name.toLowerCase().replace(/\s*\((?:left|right)\)$/,'');
type BowelRegion='duodenum'|'jejunum'|'ileum';
const aliases:Record<string,BowelRegion>={
 'duodenum':'duodenum','intestine duodenum':'duodenum',
 'jejunum':'jejunum','small intestine jejunum':'jejunum',
 'ileum':'ileum','small intestine ileum':'ileum','small intestine illium':'ileum',
};
const region=(p:Part)=>p.system==='digestive'?aliases[normalized(p)]:undefined;
const scoped=(p:Part)=>p.system==='digestive'&&!!(p.bowelWallSource||aliases[normalized(p)]||/^(?:proximal|middle|distal) part of (?:jejunum|ileum)$/.test(normalized(p)));
const routeFrame=(p:Part)=>nativePartFrame(p)??(p.system==='digestive'&&p.sourceId==='Duodenum'&&/^TA98PART:ZA:Duodenum:A05\.6\.02\.(?:002|005|007|008)$/.test(p.id)?'ZA:':undefined);
/** Keep an independent reference mesh out of a native bowel continuation. */
export function smallBowelRouteFrameCompatible(a:Part,b:Part):boolean{
 if(!scoped(a)&&!scoped(b))return true;
 const frame=routeFrame(a),other=routeFrame(b);
 return !a.suppressed&&!b.suppressed&&(!(frame||other)||!!frame&&frame===other);
}
/** A named donor bowel region does not establish functional mesh registration. */
export function smallBowelFunctionalFrameCompatible(kind:RelationKind,a:Part,b:Part):boolean{
 if(!['arterial','venous','innervation','supplies','drains','innervates','lymphaticDrainage','lymphaticTributaries'].includes(kind))return true;
 if(!scoped(a)&&!scoped(b))return true;
 const frame=nativePartFrame(a);
 return !a.suppressed&&!b.suppressed&&!!frame&&nativePartFrame(b)===frame;
}
function nativeWhole(p:Part):BowelRegion|undefined{
 if(p.suppressed||!region(p))return;
 const identities:Record<string,[string,string]>={
  'ZA:Jejunum':['Jejunum','jejunum'],
  'ZA:Duodenum':['Duodenum','duodenum'],
  'LOCAL:female:intestine_duodenum:183':['intestine_duodenum','duodenum'],
  'LOCAL:male:small_intestine_illium:1263':['small_intestine_illium','ileum'],
  'LOCAL:female:small_intestine_jejunum:182':['small_intestine_jejunum','jejunum'],
  'LOCAL:female:small_intestine_ileum:181':['small_intestine_ileum','ileum'],
 };
 const identity=identities[p.id];
 return identity&&p.sourceId===identity[0]&&region(p)===identity[1]?region(p):undefined;
}
export function smallBowelSourceLimitation(p:Part):string|undefined{
 if(!nativeWhole(p))return;
 return 'Named whole-bowel source exterior. Finite wall layers, mucosal relief, lumen continuity, mesenteric support, vascular junctions and tissue packing remain unverified. Collector and nerve associations describe anatomical routes through missing intermediates; they do not establish a local attachment or certify this surface.';
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Only an explicitly generated replacement in the same native frame may
 * inherit the whole organ's route. The source remains available as evidence,
 * even though its coarse exterior is suppressed in the review catalogue.
 */
export function smallBowelWallSource(p:Part,parts:Part[]):Part|undefined{
 if(p.suppressed||!p.bowelWallSource||p.system!=='digestive'||!p.pairedTissue)return;
 const source=parts.find(s=>s.id===p.bowelWallSource);
 if(!source?.suppressed||!source.replacementParts?.includes(p.id)||p.sourceId!==source.sourceId||nativePartFrame(p)!==nativePartFrame(source))return;
 const which=nativeWhole({...source,suppressed:false});
 const tail=p.id.split(`:BOWELWALL:${which}:`)[1];
 const fold=tail==='circular-folds'||/^circular-folds:(submucosa|muscularis-mucosae|lamina-propria)$/.test(tail??'');
 const regular=['serosa','subserosa','longitudinal-muscle','circular-muscle','submucosa','muscularis-mucosae','lamina-propria','epithelium'].includes(tail??'');
 const concept=fold?`BOWEL-FOLD:${which}:${tail==='circular-folds'?'epithelium':tail?.split(':').at(-1)}`:which==='duodenum'&&['serosa','subserosa'].includes(tail??'')?`BOWEL-PROVISIONAL:${which}:${tail}`:`SIW:${which}:${tail}`;
 if(!which||!fold&&!regular||p.conceptId!==concept||!p.id.includes(`:BOWELWALL:${which}:`))return;
 return source;
}
/** OpenStax A&P 23.5 and NCBI NBK532263: extrinsic small-bowel supply.
 * Jejunal/ileal lymph: https://pmc.ncbi.nlm.nih.gov/articles/PMC9748447/
 * Named collectors and extrinsic nerves are not direct tissue junctions.
 */
export function smallBowelRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const previous=cache.get(parts);if(previous)return previous;
 const graph=new Map<string,ResolvedRelation[]>();
 const pair=(bowel:Part,target:Part,kind:'arterial'|'venous'|'innervation'|'lymphaticDrainage',via:string[],note:string)=>{
  const reverse=kind==='arterial'?'supplies':kind==='venous'?'drains':kind==='innervation'?'innervates':'lymphaticTributaries';
  for(const [from,to,k] of [[bowel,target,kind],[target,bowel,reverse]] as const){
   const links=graph.get(from.id)??[];links.push({kind:k,target:to,via,note});graph.set(from.id,links);
  }
 };
 for(const bowel of parts){
  const source=smallBowelWallSource(bowel,parts),which=nativeWhole(source?{...source,suppressed:false}:bowel),frame=nativePartFrame(bowel);if(!which||!frame)continue;
  const routeNote=source?' Organ-level association inherited from the native whole bowel; it does not establish a direct attachment to this tissue interval.':'';
  for(const target of parts){
   if(target.suppressed||nativePartFrame(target)!==frame)continue;
   const name=normalized(target),title=which==='jejunum'?'Jejunal':'Ileal';
   if(which==='duodenum'){
    const note='Duodenal organ-level route; regional distribution, individual branches, wall/lumen continuity and tissue attachment remain unverified.'+routeNote;
    if(target.system==='arterial'){
     if(/^(?:(?:anterior|posterior) )?(?:superior|inferior) pancreaticoduodenal artery$/.test(name))pair(bowel,target,'arterial',['Duodenal branches'],note);
     if(name==='gastroduodenal artery')pair(bowel,target,'arterial',['Superior pancreaticoduodenal arteries','Duodenal branches'],note);
     if(name==='superior mesenteric artery')pair(bowel,target,'arterial',['Inferior pancreaticoduodenal artery','Duodenal branches'],note);
    }
    if(target.system==='venous'){
     if(/^(?:(?:anterior|posterior) )?(?:superior|inferior) pancreaticoduodenal vein$/.test(name))pair(bowel,target,'venous',['Duodenal venules'],note);
     if(['superior mesenteric vein','hepatic portal vein'].includes(name))pair(bowel,target,'venous',['Pancreaticoduodenal veins'],note+' Collector routes vary; this is not a claim that every region drains directly into both collectors.');
    }
    if(target.system==='nervous'&&['vagus nerve','vagus nerve (x)'].includes(name))pair(bowel,target,'innervation',['Celiac / superior mesenteric autonomic plexuses','Enteric plexuses'],note+' Vagal influence is mediated through autonomic and enteric pathways, not a direct nerve-to-layer junction.');
    continue;
   }
   if(target.system==='arterial'&&name==='superior mesenteric artery')pair(bowel,target,'arterial',[title+' arteries'],'Supply through jejunal/ileal branches and their intestinal arcades. The named source collector is not a direct bowel-wall junction; individual branches, walls, lumens and source packing remain unverified.'+routeNote);
   if(target.system==='venous'&&name==='superior mesenteric vein')pair(bowel,target,'venous',[title+' veins'],'Drainage through intestinal tributaries to the superior mesenteric vein. Tributaries, continuous vessel lumens and source bowel/vascular attachments remain unverified.'+routeNote);
   if(target.system==='lymphatic'){
    const note='Organ-level jejunal/ileal lymph drainage through wall lymphatic plexuses and mesenteric collectors. This does not attach a node directly to this tissue surface; individual lacteals, collectors, regional node membership and source junctions remain unverified.'+routeNote;
    if(name==='juxta-intestinal mesenteric nodes')pair(bowel,target,'lymphaticDrainage',['Intestinal wall lymphatic plexuses','Mesenteric lymphatic collectors'],note);
    if(['superior mesenteric lymph nodes','central superior mesenteric nodes'].includes(name))pair(bowel,target,'lymphaticDrainage',['Intestinal wall lymphatic plexuses','Juxta-intestinal mesenteric nodes','Intermediate mesenteric nodes'],note+' The grouped collector is not a claim that every node receives this organ directly.');
   }
   if(target.system==='nervous'&&['vagus nerve','vagus nerve (x)'].includes(name))pair(bowel,target,'innervation',['Superior mesenteric plexus','Enteric plexuses'],'Extrinsic vagal parasympathetic influence reaches the small bowel through autonomic and enteric pathways. This association is not a direct cranial-nerve-to-wall attachment; intermediate plexuses and source tissue junctions remain unverified.'+routeNote);
  }
 }
 cache.set(parts,graph);return graph;
}
