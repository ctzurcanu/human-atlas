/** Original illustrative intestinal crypt specimen, metres. Three envelopes
 * bound epithelial and stromal cells; the lumen is empty and the mouth open.
 * The outer sleeve and mouth annuli are artificial specimen boundaries.
 */
export function cryptDetailGeometry(){
 const sides=16,profile=[[0,.04],[.08,.043],[.15,.042],[.20,.04]];
 const profiles=[0,.015,.055].map(offset=>profile.map(([depth,r])=>[-depth,r+offset]));
 for(const angle of [Math.PI/6,Math.PI/3,Math.PI*5/12])for(let k=0;k<3;k++){
  const r=.04+[0,.015,.055][k];profiles[k].push([-.20-r*Math.sin(angle),r*Math.cos(angle)]);
 }
 const count=profiles[0].length*sides+1,tip=count-1,interfaces=profiles.map((p,k)=>{
  const bank=new Float32Array(count*3);
  for(let i=0;i<p.length;i++)for(let j=0;j<sides;j++){const a=j/sides*2*Math.PI;bank.set([p[i][1]*Math.cos(a)*.001,p[i][0]*.001,p[i][1]*Math.sin(a)*.001],(i*sides+j)*3);}
  bank.set([0,(-.24-[0,.015,.055][k])*.001,0],tip*3);return bank;
 });
 const faces:number[]=[],rims:number[]=[];
 for(let k=0;k<profiles[0].length-1;k++)for(let j=0;j<sides;j++){const a=k*sides+j,b=k*sides+(j+1)%sides,A=a+sides,B=b+sides;faces.push(a,b,A,b,B,A);}
 for(let j=0;j<sides;j++){const a=(profiles[0].length-1)*sides+j,b=(profiles[0].length-1)*sides+(j+1)%sides;faces.push(a,b,tip);rims.push(j,(j+1)%sides);}
 return {lumen:interfaces[0],epithelium:interfaces[1],stroma:interfaces[2],faces:new Uint32Array(faces),rims:new Uint32Array(rims),sides,count,depthMeters:.000295,epithelialThicknessMeters:.000015,stromaThicknessMeters:.00004};
}
