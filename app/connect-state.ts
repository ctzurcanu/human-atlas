import type {SceneState} from './anatomy';
import type {ToolDocument} from './advanced-tools';
import type {RelationshipExpansion} from './relationship-selection';
import type {ExpandedControl} from './connect-dom';
import type {Terminology} from './anatomical-terminology';

export type ConnectRole='host'|'guest';
export type ConnectPeer={id:string;ip:string;name:string;online:boolean;inSession?:boolean};
export type SceneWorldPose={position:number[];target:number[];up:number[];quaternion:number[];fov:number;zoom:number;near:number;far:number;aspect?:number;filmGauge?:number;filmOffset?:number;view?:{enabled:boolean;fullWidth:number;fullHeight:number;offsetX:number;offsetY:number;width:number;height:number}|null};
export type ScenePose={model:string;camera:number[];up:number[];rotation:number[];world?:SceneWorldPose};
export type AdvancedPresentation={document:ToolDocument;slideIndex:number|null};
type ChoiceSummary={id:string;name:string;elements:string[];group:boolean;terminology?:Terminology};
export type ConnectSnapshot={
 version:1;model:string;hierarchy:string;guestSources:string[];state:SceneState;
 choice:(ChoiceSummary&{children?:ChoiceSummary[]})|null;
 ui:{panel:'search'|'advanced'|'sections'|'visit'|null;frontPanel:'layers'|'search'|'advanced'|'sections'|'visit'|'details';details:boolean;layersVisible:boolean;mobileLayersOpen:boolean;addSelection:boolean;query:string;expanded:ExpandedControl[]};
 covering:string[];relationshipExpansion:RelationshipExpansion|null;pose:ScenePose|null;advanced:AdvancedPresentation|null;
};
