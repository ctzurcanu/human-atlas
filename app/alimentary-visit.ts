import type {Atlas} from './anatomy';
export type AlimentaryVisit={points:[number,number,number][];segments:{name:string;startIndex:number;points:[number,number,number][]}[];bindings:{id:string;positions:number;indices:number;vertexCount:number;chunkUrl:string;bounds:number[][]}[]};
export function matchesAlimentaryVisit(atlas:Atlas,value:unknown):value is AlimentaryVisit{
 if(!value||typeof value!=='object')return false;const route=value as AlimentaryVisit;if(!Array.isArray(route.points)||!Array.isArray(route.bindings)||!Array.isArray(route.segments))return false;
 return route.points.length>10&&route.points.every(p=>p.length===3&&p.every(Number.isFinite))&&route.bindings.every(binding=>{const p=atlas.parts.find(p=>p.id===binding.id);return p&&!p.suppressed&&p.positions===binding.positions&&p.indices===binding.indices&&p.vertexCount===binding.vertexCount&&atlas.chunks[p.chunk]?.url===binding.chunkUrl&&JSON.stringify(p.bounds)===JSON.stringify(binding.bounds);});
}
