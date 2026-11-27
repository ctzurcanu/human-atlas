/** A lower-opacity neighbor must not erase a more visible tissue boundary.
 * Equal-opacity assembled tissues omit their internal shared interface. */
export function tissueNeighborOccludes(sourceOpacity:number,neighborOpacity:number):boolean{
 return Number.isFinite(sourceOpacity)&&Number.isFinite(neighborOpacity)&&sourceOpacity>=.001&&neighborOpacity>=.001&&neighborOpacity+1e-6>=sourceOpacity;
}

const zeroOrigin:readonly number[]=[0,0,0];
/** Coordinate-bank origins do not separate anatomical neighbours. Compare the
 * actual affine tissue poses after removing each bank's intrinsic origin. */
export function tissuePosesCoincide(first:{elements:ArrayLike<number>},firstOrigin:readonly number[]|undefined,
 second:{elements:ArrayLike<number>},secondOrigin:readonly number[]|undefined,tolerance=1e-7):boolean{
 const a=first.elements,b=second.elements,oa=firstOrigin??zeroOrigin,ob=secondOrigin??zeroOrigin;
 if(a.length!==16||b.length!==16||oa.length!==3||ob.length!==3||!(tolerance>=0)||!Number.isFinite(tolerance))return false;
 for(let i=0;i<16;i++){
  if(!Number.isFinite(a[i])||!Number.isFinite(b[i]))return false;
  if((i<12||i===15)&&Math.abs(a[i]-b[i])>tolerance)return false;
 }
 for(let axis=0;axis<3;axis++){
  if(!Number.isFinite(oa[axis])||!Number.isFinite(ob[axis]))return false;
  const ta=a[12+axis]-a[axis]*oa[0]-a[4+axis]*oa[1]-a[8+axis]*oa[2];
  const tb=b[12+axis]-b[axis]*ob[0]-b[4+axis]*ob[1]-b[8+axis]*ob[2];
  if(Math.abs(ta-tb)>tolerance)return false;
 }
 return true;
}
