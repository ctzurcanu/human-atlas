import {useCallback,useReducer} from 'react';
import type {SceneState} from './anatomy';
type Update=SceneState|((s:SceneState)=>SceneState);
type History={present:SceneState;past:SceneState[];future:SceneState[];time:number;key:string};
type Action=({type:'set';update:Update;time:number}|{type:'load';state:SceneState}|{type:'undo'|'redo'})&{camera?:number[]};
export function historyReducer(h:History,a:Action):History{
 if(a.type==='load')return {present:a.state,past:[],future:[],time:0,key:''};
 if(a.type==='undo'){const present=h.past.at(-1);return present?{present,past:h.past.slice(0,-1),future:[a.camera?{...h.present,camera:a.camera}:h.present,...h.future],time:0,key:''}:h;}
 if(a.type==='redo'){const present=h.future[0];return present?{present,past:[...h.past,a.camera?{...h.present,camera:a.camera}:h.present],future:h.future.slice(1),time:0,key:''}:h;}
 if(a.type!=='set')return h;const next=typeof a.update==='function'?a.update(h.present):a.update;if(JSON.stringify(next)===JSON.stringify(h.present))return h;
 const keys=Object.keys(next).filter(k=>JSON.stringify(next[k as keyof SceneState])!==JSON.stringify(h.present[k as keyof SceneState]));const key=keys.length===1&&['contextOpacity','skinOpacity','section','explode'].includes(keys[0])?keys[0]:'';
 return {present:next,past:key&&key===h.key&&a.time-h.time<800?h.past:[...h.past,a.camera?{...h.present,camera:a.camera}:h.present].slice(-80),future:[],time:a.time,key};
}
export function useViewHistory(initial:SceneState,getCamera?:()=>number[]|undefined){
 const [h,dispatch]=useReducer(historyReducer,{present:initial,past:[],future:[],time:0,key:''});
 return {state:h.present,setState:useCallback((update:Update)=>dispatch({type:'set',update,time:Date.now(),camera:getCamera?.()}),[getCamera]),loadState:useCallback((state:SceneState)=>dispatch({type:'load',state}),[]),undo:()=>dispatch({type:'undo',camera:getCamera?.()}),redo:()=>dispatch({type:'redo',camera:getCamera?.()}),canUndo:!!h.past.length,canRedo:!!h.future.length};
}
