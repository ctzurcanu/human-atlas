const clamp=(value:number)=>Math.max(0,Math.min(1,value));

export function smoothStep(value:number){const t=clamp(value);return t*t*t*(10+t*(-15+6*t));}

/** Exact critically damped motion, with velocity retained when the target changes. */
export function dampMotion(value:number,velocity:number,target:number,response:number,dt:number){
 const omega=2/Math.max(.001,response),offset=value-target,j=velocity+omega*offset,decay=Math.exp(-omega*Math.max(0,dt));
 return {value:target+(offset+j*dt)*decay,velocity:(velocity-omega*j*dt)*decay};
}

/** Keep a fixed track, giving the first two intervals twice their old width. */
export function explosionMarks(steps:number){
 const count=Math.max(1,Math.floor(steps));
 if(count<=2)return Array.from({length:count+1},(_,i)=>i/count);
 const first=count>4?2/count:2/(count+2),rest=(1-2*first)/(count-2);
 return Array.from({length:count+1},(_,i)=>i<=2?i*first:2*first+(i-2)*rest);
}
export function explosionToSlider(amount:number,steps:number){
 const marks=explosionMarks(steps),progress=clamp(amount)*(marks.length-1),from=Math.min(marks.length-2,Math.floor(progress));
 return marks[from]+(marks[from+1]-marks[from])*(progress-from);
}
export function sliderToExplosion(position:number,steps:number,snap=false){
 const marks=explosionMarks(steps);let value=clamp(position);
 if(snap){
  let nearest=0;for(let i=1;i<marks.length;i++)if(Math.abs(marks[i]-value)<Math.abs(marks[nearest]-value))nearest=i;
  const width=nearest===0?marks[1]:nearest===marks.length-1?1-marks[nearest-1]:Math.min(marks[nearest]-marks[nearest-1],marks[nearest+1]-marks[nearest]);
  if(Math.abs(value-marks[nearest])<=Math.min(.04,width*.28))value=marks[nearest];
 }
 let from=0;while(from<marks.length-2&&value>marks[from+1])from++;
 return (from+(value-marks[from])/(marks[from+1]-marks[from]))/(marks.length-1);
}
