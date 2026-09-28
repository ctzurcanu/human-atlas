export function validCamera(value){
 if(!Array.isArray(value)||![6,8,12].includes(value.length)||!value.every(item=>typeof item==='number'&&Number.isFinite(item)&&Math.abs(item)<10000))return false;
 if(value.length!==12)return Math.hypot(value[0]-value[3],value[1]-value[4],value[2]-value[5])>.000001;
 if([...value.slice(0,3),...value.slice(6,12)].some(item=>item<0||item>1)||value[8]<=0||value[8]>=1)return false;
 return Math.hypot(...value.slice(0,3).map(item=>2*item-1))>.001&&Math.hypot(...value.slice(9,12).map(item=>2*item-1))>.001;
}
