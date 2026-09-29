/** Integrate controller pointing changes while a separate trigger is held. */
export function orbitAngleDelta(current:number,previous:number){
 return Math.atan2(Math.sin(current-previous),Math.cos(current-previous));
}

export function advanceControllerOrbit(previous:{yaw:number;pitch:number}|null,current:{yaw:number;pitch:number},orbit:{yaw:number;pitch:number}){
 if(!previous)return orbit;
 const yaw=orbitAngleDelta(current.yaw,previous.yaw),pitch=current.pitch-previous.pitch;
 // Ignore tracking jumps; clamp ordinary fast wrist turns to a stable step.
 if(Math.abs(yaw)>.6||Math.abs(pitch)>.6)return orbit;
 return {yaw:orbit.yaw-Math.max(-.12,Math.min(.12,yaw)),pitch:Math.max(-1.25,Math.min(1.25,orbit.pitch-Math.max(-.12,Math.min(.12,pitch))))};
}
