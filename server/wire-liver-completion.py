"""Narrow, checked edits for the existing local viewer integration."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def change(file,old,new):
    p=ROOT/file;s=p.read_text();assert old in s,(file,old[:100]);p.write_text(s.replace(old,new,1))
def imp(file,text):
    p=ROOT/file;s=p.read_text();assert text not in s;p.write_text(text+'\n'+s)
imp('app/liver-reference-selection.ts',"import {liverCompletionRegistry} from './liver-completion-controls';")
change('app/liver-reference-selection.ts',"...(liverLandmarkRegistry(atlas)?.marker.termGroups??[])];","...(liverLandmarkRegistry(atlas)?.marker.termGroups??[]),...(liverCompletionRegistry(atlas)?.termGroups.filter(g=>g.curves)??[])];")
change('app/liver-surface-features.ts',"m.termGroups.length!==terms.length||new Set(m.termGroups.map(g=>g.term)).size!==terms.length","![terms.length,terms.length+1].includes(m.termGroups.length)||new Set(m.termGroups.map(g=>g.term)).size!==m.termGroups.length")
change('app/liver-surface-features.ts',"!terms.includes(g.term)","![...terms,'A05.8.01.017'].includes(g.term)")
imp('app/anatomy-label-plan.ts',"import {liverCompletionLabelReference} from './liver-completion-controls';")
change('app/anatomy-label-plan.ts',"??liverLandmarkLabelReference(atlas,choice));","??liverLandmarkLabelReference(atlas,choice)??liverCompletionLabelReference(atlas,choice));")
imp('app/liver-system.ts',"import {liverCompletionTermParts,liverWholePresentationParts} from './liver-completion-controls';")
change('app/liver-system.ts',"...liverLandmarkTermParts(atlas)]","...liverLandmarkTermParts(atlas),...liverCompletionTermParts(atlas)]")
change('app/liver-system.ts',"const parts=unique([...(direct.get(term.id)??[]),...children.flatMap(c=>c.parts)]),","const parts=liverWholePresentationParts(atlas,term.id,unique([...(direct.get(term.id)??[]),...children.flatMap(c=>c.parts)])),")
imp('app/guest-hierarchy.ts',"import {liverCompletionTermParts,liverWholePresentationParts} from './liver-completion-controls';")
change('app/guest-hierarchy.ts',"const cricoidTerms=", "for(const [term,parts]of hierarchy.id==='ta98'?liverCompletionTermParts(atlas):[])replacementTermParts.set(term,parts);\n const cricoidTerms=")
change('app/guest-hierarchy.ts',"const parts=children.every(child=>child.parts===directParts)?directParts:[...new Map([...directParts,...children.flatMap(child=>child.parts)].map(part=>[part.id,part])).values()];", "const parts=liverWholePresentationParts(atlas,selectors.id??'',children.every(child=>child.parts===directParts)?directParts:[...new Map([...directParts,...children.flatMap(child=>child.parts)].map(part=>[part.id,part])).values()]);")
imp('app/preferred-source-representations.ts',"import {liverCompletionRegistry} from './liver-completion-controls';")
change('app/preferred-source-representations.ts',"new Set(hepaticConnectionRegistry(atlas)?.marker.regionPartIds??[])","new Set([...(hepaticConnectionRegistry(atlas)?.marker.regionPartIds??[]),...(liverCompletionRegistry(atlas)?.partIds??[])])")
change('app/viewer-state.ts',"if(p.suppressed)return false;","if(p.suppressed)return false;\n if(liverAlternativeActive(p,s))return false;")
change('app/viewer-state.ts',"if(p.suppressed||s.depth!==undefined&&sets.depthHidden.has(depthLayerFor(p)))return false;","if(p.suppressed||s.depth!==undefined&&sets.depthHidden.has(depthLayerFor(p)))return false;\n if(liverAlternativeActive(p,s))return false;")
change('app/viewer-state.ts',"export function partVisible(","""function liverAlternativeActive(p:Part,s:SceneState){
 const ids=(p as Part&{liverAlternativeParts?:string[]}).liverAlternativeParts,sets=visibilitySets(s);
 return p.id.startsWith('TA98PART:ZA:Liver:')&&!!ids?.length&&ids.every(id=>/^TA98LIVER:male:segment:(?:I|II|III|IV|V|VI|VII|VIII)$/.test(id))
  &&ids.some(id=>sets.selected.has(id)||!s.isolate&&!sets.hidden.has(id)&&sets.visible.has('digestive'));
}
export function partVisible(""")
change('app/local-ta98-layer-catalogue.ts'," if(model==='local-ta98'&&query.get('hepatic')==='features-v143')", " if(model==='local-ta98'&&query.get('hepatic')==='completion-v154')return '/local-models/geometry-workbench/liver-completion/completion-v154/atlas.json';\n if(model==='local-ta98'&&query.get('hepatic')==='features-v143')")
change('app/page.tsx',"|'local-female-ta98';","|'local-female-ta98'|'liver-unit';")
change('app/page.tsx'," {value:'cell',label:'Human eukaryotic cell'},"," {value:'cell',label:'Human eukaryotic cell'},\n ...(import.meta.env.DEV?[{value:'liver-unit' as const,label:'Liver · microscopic unit (1000×)'}]:[]),")
change('app/page.tsx',"const cataloguePath=(model:Model)=>", "const cataloguePath=(model:Model)=>model==='liver-unit'?'/local-models/geometry-workbench/liver-completion/completion-v154/microscopic-atlas.json':")
imp('app/selection-inspector.tsx',"import {liverMicroscopicDetailUrl} from './liver-completion-controls';")
change('app/selection-inspector.tsx',"const microscopicDetail=", "const liverMicroDetail=import.meta.env.DEV?liverMicroscopicDetailUrl(atlas,choice):null;\n const microscopicDetail=")
change('app/selection-inspector.tsx',"{microscopicDetail&&", "{liverMicroDetail&&<p className=\"context-note\"><a href={liverMicroDetail} target=\"_blank\" rel=\"noopener noreferrer\">Open magnified liver unit</a></p>}\n   {microscopicDetail&&")
print('Liver completion controls wired')
