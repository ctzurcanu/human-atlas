import {readFile} from 'node:fs/promises';
import {join} from 'node:path';
import type {Atlas} from '../app/anatomy';
import {availableModels,combinedTa98,maleTa98Path,modelCataloguePath,type ItemCount} from '../app/model-catalogues';
import {BOWEL_WALL_REVIEW_OPTIONS} from '../app/bowel-wall-review-navigation';
import {withPreferredSourceRepresentations,withHiddenReviewAlternatives} from '../app/preferred-source-representations';
import {withSystemGeometry} from '../app/system-geometry';

/** Build a small menu index without downloading other models in the browser. */
export async function menuItemCounts(root:string, dev:boolean, maleDirectory:string|undefined, embryo?:string) {
 const entries = availableModels(dev,!!maleDirectory).map(({value})=>({
  id:String(value), model:value, url:modelCataloguePath(value,new URLSearchParams(),dev,embryo),
 }));
 if (dev) for (const option of BOWEL_WALL_REVIEW_OPTIONS) {
  const model = option.sex==='male' ? 'local-ta98' : 'local-female-ta98';
  entries.push({id:option.value,model,url:modelCataloguePath(model,new URLSearchParams({bowel:'prototype'}),true)});
 }
 const counts:Record<string,ItemCount> = {};
 const cache = new Map<string,ItemCount>();
 for (const {id,model,url} of entries) {
  const file = url===maleTa98Path && maleDirectory ? join(maleDirectory,'atlas.json') :
   url.startsWith('/local-models/') ? join(root,'.local-models',url.slice(14)) : join(root,'public',url);
  let count = cache.get(file);
  if (!count) {
   let bytes;
   try { bytes = await readFile(file,'utf8'); }
   catch (error) { if ((error as NodeJS.ErrnoException).code==='ENOENT') continue; throw error; }
   let atlas:Atlas = withPreferredSourceRepresentations(withSystemGeometry(JSON.parse(bytes)));
   if (combinedTa98(model)) atlas=withHiddenReviewAlternatives(atlas);
   count={total:atlas.parts.length,selectable:atlas.parts.filter(p=>!p.suppressed).length};
   cache.set(file,count);
  }
  counts[id]=count;
 }
 return counts;
}
