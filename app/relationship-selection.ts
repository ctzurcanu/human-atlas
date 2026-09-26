export interface RelationshipExpansion {
 anchorId:string;
 selectedKey:string;
 frontier:string[];
 exhausted:boolean;
 depth:number;
}

// Expand one relationship step from the pieces added by the preceding click.
// Keeping the frontier separate prevents cycles from selecting the same piece twice.
export function expandRelationshipSelection(
 selected:string[],
 previous:RelationshipExpansion|null,
 related:(id:string)=>Iterable<string>,
):{selected:string[];expansion:RelationshipExpansion|null}{
 if(!selected.length)return {selected,expansion:null};
 const continuing=previous?.selectedKey===selected.join('|');
 const anchorId=continuing?previous.anchorId:selected[0];
 const frontier=continuing?previous.frontier:[anchorId];
 const seen=new Set(selected);
 const added:string[]=[];
 for(const sourceId of frontier){
  for(const targetId of related(sourceId)){
   if(seen.has(targetId))continue;
   seen.add(targetId);
   added.push(targetId);
  }
 }
 const expanded=added.length?[...selected,...added]:selected;
 return {selected:expanded,expansion:{anchorId,selectedKey:expanded.join('|'),frontier:added,exhausted:added.length===0,depth:(continuing?previous.depth:0)+1}};
}
