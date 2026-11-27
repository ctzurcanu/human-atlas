import type {Part,SceneState} from './anatomy';
import {depthLayerFor} from './depth-layers';
import {partLayerOpacity} from './depth-control';

/** Checkbox visibility changes retain every camera and placement field. */
export function toggleSystemParts(current:SceneState,parts:Part[],availableParts:Part[]):SceneState{
 const ids=new Set(parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]);
 const anyOn=parts.some(part=>current.visible.includes(part.system)&&!hiddenNow.has(part.id)&&!current.depthHidden?.includes(depthLayerFor(part))&&partLayerOpacity(part,current)>0);
 const newlyEnabled=new Set(parts.filter(part=>!current.visible.includes(part.system)).map(part=>part.system));
 const restored=new Set(parts.map(depthLayerFor).filter(id=>current.depthHidden?.includes(id)));
 if(!anyOn)for(const part of availableParts)if((newlyEnabled.has(part.system)||restored.has(depthLayerFor(part)))&&!ids.has(part.id))hiddenNow.add(part.id);
 if(anyOn)for(const id of ids)hiddenNow.add(id);else for(const id of ids)hiddenNow.delete(id);
 const systems=[...new Set(parts.map(part=>part.system))];
 return {...current,visible:anyOn?current.visible:[...new Set([...current.visible,...systems])],hidden:[...hiddenNow],depthHidden:anyOn?current.depthHidden:current.depthHidden?.filter(id=>!restored.has(id as ReturnType<typeof depthLayerFor>)),skinOpacity:!anyOn&&parts.some(part=>depthLayerFor(part)==='skin')&&!(current.skinOpacity??0)?.1:current.skinOpacity,selected:anyOn?current.selected.filter(id=>!ids.has(id)):current.selected,isolate:false};
}
