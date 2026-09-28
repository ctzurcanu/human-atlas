import type {SceneState} from './anatomy';
import type {ToolDocument} from './advanced-tools';
import type {RelationshipExpansion} from './relationship-selection';
import type {ExpandedControl} from './connect-dom';
import type {Terminology} from './anatomical-terminology';

export type ConnectRole='host'|'guest';
export type ConnectPeer={id:string;ip:string;name:string;online:boolean};
export type ScenePose={model:string;camera:number[];up:number[];rotation:number[]};
export type AdvancedPresentation={document:ToolDocument;slideIndex:number|null};
type ChoiceSummary={id:string;name:string;elements:string[];group:boolean;terminology?:Terminology};
export type ConnectSnapshot={
 version:1;model:string;hierarchy:string;guestSources:string[];state:SceneState;
 choice:(ChoiceSummary&{children?:ChoiceSummary[]})|null;
 ui:{panel:'search'|'advanced'|'sections'|null;frontPanel:'layers'|'search'|'advanced'|'sections'|'details';details:boolean;layersVisible:boolean;mobileLayersOpen:boolean;addSelection:boolean;query:string;expanded:ExpandedControl[]};
 covering:string[];relationshipExpansion:RelationshipExpansion|null;pose:ScenePose|null;advanced:AdvancedPresentation|null;
};
