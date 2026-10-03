/** Detect baked annotation glyphs, preserving source territory pixels verbatim.
 * Returns sampling metadata, not an edited color image. Large black UV voids
 * are excluded. Five-pixel bounds cover the glyph's antialiased perimeter.
 */
export function dermatomeInkMask(rgba:Uint8ClampedArray,width:number,height:number){
 const visited=new Uint8Array(width*height),mask=new Uint8Array(width*height),queue=new Int32Array(width*height);let annotations=0;
 const darkPixels=new Uint8Array(width*height);
 for(let i=0;i<darkPixels.length;i++){const r=rgba[i*4],g=rgba[i*4+1],b=rgba[i*4+2],max=Math.max(r,g,b),min=Math.min(r,g,b);darkPixels[i]=max<95||max<190&&max-min<30?1:0;}
 const dark=(i:number)=>!!darkPixels[i];
 const pigment=(x:number,y:number)=>x>=0&&x<width&&y>=0&&y<height&&!dark(y*width+x);
 const mark=(x:number,y:number,padding=4)=>{for(let yy=Math.max(0,y-padding);yy<=Math.min(height-1,y+padding);yy++)mask.fill(255,yy*width+Math.max(0,x-padding),yy*width+Math.min(width,x+padding+1));};
 for(let seed=0;seed<visited.length;seed++){
  if(visited[seed]||!dark(seed))continue;
  let start=0,end=1,minX=seed%width,maxX=minX,minY=Math.floor(seed/width),maxY=minY;queue[0]=seed;visited[seed]=1;
  while(start<end){const i=queue[start++],x=i%width,y=Math.floor(i/width);minX=Math.min(minX,x);maxX=Math.max(maxX,x);minY=Math.min(minY,y);maxY=Math.max(maxY,y);
   for(const [dx,dy] of [[-1,0],[1,0],[0,-1],[0,1]]){const nx=x+dx,ny=y+dy,j=ny*width+nx;if(nx<0||nx>=width||ny<0||ny>=height||visited[j]||!dark(j))continue;visited[j]=1;queue[end++]=j;}
  }
  // Unused UV islands have large solid black backgrounds. Tiny dark source
  // pigment is not discarded; this threshold only isolates annotation ink.
  if(end>12000||maxX-minX>180||maxY-minY>180){
   // Leaders can touch the unused black UV background. Recover only narrow
   // strokes with actual pigment on opposing sides, rather than masking the void.
   for(let j=0;j<end;j++){const i=queue[j],x=i%width,y=Math.floor(i/width);let thin=false;
    for(const radius of [2,4,6,8]){for(const [dx,dy] of [[1,0],[0,1],[1,1],[1,-1]])if(pigment(x+dx*radius,y+dy*radius)&&pigment(x-dx*radius,y-dy*radius)){thin=true;break;}if(thin)break;}
    if(thin)mark(x,y);
   }continue;
  }
  annotations++;
  for(let y=Math.max(0,minY-5);y<=Math.min(height-1,maxY+5);y++)mask.fill(255,y*width+Math.max(0,minX-5),y*width+Math.min(width,maxX+6));
 }
 return {mask,annotations};
}
