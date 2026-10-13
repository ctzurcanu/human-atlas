/** Original idealized villus fragment, in metres. Not a registered organ source.
 * One outer/inner carrier; the basal annulus is an exposed specimen cut.
 * Capillaries, lacteal, cells, crypts and brush border are not represented.
 */
export function villusDetailGeometry(){
 const sides=16,profile=[[0,.07],[.10,.073],[.25,.067],[.45,.065],[.60,.065],[.66,.065]],thickness=.025;
 const outerProfile=profile.map(([y,r])=>[y,r+thickness]);
 for(const angle of [Math.PI/6,Math.PI/3,Math.PI*5/12]){profile.push([.66+.065*Math.sin(angle),.065*Math.cos(angle)]);outerProfile.push([.66+.09*Math.sin(angle),.09*Math.cos(angle)]);}
 const count=profile.length*sides+1,tip=count-1,outer=new Float32Array(count*3),inner=new Float32Array(count*3),faces:number[]=[],rims:number[]=[];
 for(let k=0;k<profile.length;k++)for(let j=0;j<sides;j++){
  const angle=j/sides*Math.PI*2,[y,r]=profile[k],[Y,R]=outerProfile[k],i=(k*sides+j)*3;
  inner.set([r*Math.cos(angle)*.001,y*.001,r*Math.sin(angle)*.001],i);
  outer.set([R*Math.cos(angle)*.001,Y*.001,R*Math.sin(angle)*.001],i);
  if(k<profile.length-1){const a=k*sides+j,b=k*sides+(j+1)%sides,A=a+sides,B=b+sides;faces.push(a,A,b,b,A,B);}
 }
 inner.set([0,.725*.001,0],tip*3);outer.set([0,.750*.001,0],tip*3);
 for(let j=0;j<sides;j++){const a=(profile.length-1)*sides+j,b=(profile.length-1)*sides+(j+1)%sides;faces.push(a,tip,b);rims.push((j+1)%sides,j);}
 const coreFaces=[...faces];for(let j=0;j<sides;j++)coreFaces.push((j+1)%sides,count,j);
 return {outer,inner,faces:new Uint32Array(faces),rims:new Uint32Array(rims),coreFaces:new Uint32Array(coreFaces),sides,count,epitheliumRadialThicknessMeters:thickness*.001,heightMeters:.000750};
}
/** Exact containment for the untwisted regular polygon rings of this fragment.
 * Used for picking stencil caps; a large support quad is not itself tissue.
 */
export function insideVillusCarrier(point:readonly number[],positions:Float32Array,sides:number){
 const rings=(positions.length/3-1)/sides,tip=positions.length-3,direction=positions[tip+1]>positions[1]?1:-1,y=point[1]*direction;
 if(y<positions[1]*direction||y>positions[tip+1]*direction)return false;
 for(let k=0;k<rings;k++){
  const a=k*sides*3,b=k===rings-1?tip:(k+1)*sides*3,lo=positions[a+1]*direction,hi=positions[b+1]*direction;if(y>hi)continue;
  const t=(y-lo)/(hi-lo),r=positions[a]*(1-t)+positions[b]*t,sector=2*Math.PI/sides;
  const angle=((Math.atan2(point[2],point[0])%sector)+sector)%sector-sector/2;
  return Math.hypot(point[0],point[2])<=r*Math.cos(sector/2)/Math.cos(angle)+1e-12;
 }return false;
}
