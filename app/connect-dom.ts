export type ExpandedControl={key:string;expanded:boolean};
const controls=()=>[...document.querySelectorAll<HTMLButtonElement>('#anatomy-browser button[data-connect-key],.detail-sheet button[data-slot="accordion-trigger"]')];
const key=(button:HTMLButtonElement)=>button.dataset.connectKey??`details:${button.textContent?.trim().replace(/\s+/g,' ')}`;
export function expandedControls():ExpandedControl[]{return controls().map(button=>({key:key(button),expanded:button.getAttribute('aria-expanded')==='true'}));}
/** Parents must open before their children exist, so callers retry after a frame. */
export function applyExpandedControls(values:ExpandedControl[]){
 const byKey=new Map(controls().map(button=>[key(button),button]));let pending=false;
 for(const value of values){const button=byKey.get(value.key);if(!button){pending=true;continue;}if((button.getAttribute('aria-expanded')==='true')!==value.expanded){button.click();pending=true;}}
 return pending;
}
