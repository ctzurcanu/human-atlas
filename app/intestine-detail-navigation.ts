import type {HierarchyChoice} from './hierarchy-choice';
/** A separate magnified illustration never supplies a whole-organ mesh. */
export function intestinalVillusDetailUrl(choice:HierarchyChoice):string|null{
 const id=choice.id.replace(/^hierarchy:guest:ta98:/,'');
 const match=/^SIW:(duodenum|jejunum|ileum):villi$/.exec(id);
 if(match)return `./intestinal-villus-review.html?organ=${match[1]}`;
 return choice.terminology?.ta98==='A05.6.01.011'||choice.ta98Term==='A05.6.01.011'?'./intestinal-villus-review.html':null;
}
export function intestinalCryptDetailUrl(choice:HierarchyChoice):string|null{
 const id=choice.id.replace(/^hierarchy:guest:ta98:/,'');
 const match=/^SIW:(duodenum|jejunum|ileum):intestinal-glands$/.exec(id);
 if(match)return `./intestinal-crypt-review.html?organ=${match[1]}`;
 return choice.terminology?.ta98==='A05.6.01.012'||choice.ta98Term==='A05.6.01.012'?'./intestinal-crypt-review.html':null;
}
export function intestinalLymphoidDetailUrl(choice:HierarchyChoice):string|null{
 const id=choice.id.replace(/^hierarchy:guest:ta98:/,'');
 const match=/^SIW:(duodenum|jejunum|ileum):(solitary-nodules|aggregated-nodules)$/.exec(id);
 if(match)return `./intestinal-lymphoid-review.html?mode=${match[2]==='solitary-nodules'?'solitary':'aggregate'}&organ=${match[1]}`;
 const term=choice.terminology?.ta98??choice.ta98Term;
 return term==='A05.6.01.013'?'./intestinal-lymphoid-review.html?mode=solitary':term==='A05.6.01.014'?'./intestinal-lymphoid-review.html?mode=aggregate':null;
}
