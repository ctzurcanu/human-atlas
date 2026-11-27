/** Saved presentation preferences control general group labels. The current
 * anatomical selection always retains its identifying callout. */
export function selectedLabelsVisible(state:{labels?:boolean;selected:string[];dermatomeSelected?:string[]}):boolean{
 return state.labels!==false||state.selected.length>0||!!state.dermatomeSelected?.length;
}
