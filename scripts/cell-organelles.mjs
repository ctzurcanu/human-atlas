/**
 * Schematic organelles added to the UMCG cell cutaway.
 *
 * The source model has no peroxisomes, smooth ER, microtubules, intermediate filaments,
 * lipid droplets, glycogen granules, secretory vesicles, autophagosomes or proteasomes.
 * These are generated here as simple, deterministic shapes placed in free space between
 * the source meshes (never inside or through them). They are labelled as schematic
 * additions in the atlas; they are not part of the CC BY-NC-SA source scan.
 */

function prng(seed){
 let state=seed>>>0;
 return ()=>{
  state=(state+0x6d2b79f5)>>>0;
  let t=state;
  t=Math.imul(t^(t>>>15),t|1);
  t^=t+Math.imul(t^(t>>>7),t|61);
  return ((t^(t>>>14))>>>0)/4294967296;
 };
}

const sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]];
const add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]];
const mul=(a,s)=>[a[0]*s,a[1]*s,a[2]*s];
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const length=a=>Math.hypot(a[0],a[1],a[2]);
const unit=a=>{const l=length(a)||1;return [a[0]/l,a[1]/l,a[2]/l];};

/** Accumulates triangles with per-vertex smooth normals. */
class Mesh{
 constructor(){this.points=[];this.normals=[];this.indices=[];}
 get count(){return this.points.length/3;}
 vertex(p,n){this.points.push(p[0],p[1],p[2]);this.normals.push(n[0],n[1],n[2]);return this.count-1;}
 /** UV ellipsoid, optionally rotated by an orthonormal basis [x,y,z]. */
 ellipsoid(center,radii,basis,segments=12,rings=8){
  const base=this.count,[bx,by,bz]=basis??[[1,0,0],[0,1,0],[0,0,1]];
  for(let r=0;r<=rings;r++){
   const phi=Math.PI*r/rings;
   for(let s=0;s<=segments;s++){
    const theta=2*Math.PI*s/segments;
    const local=[Math.sin(phi)*Math.cos(theta),Math.cos(phi),Math.sin(phi)*Math.sin(theta)];
    const world=add(add(mul(bx,local[0]*radii[0]),mul(by,local[1]*radii[1])),mul(bz,local[2]*radii[2]));
    const normal=unit(add(add(mul(bx,local[0]/radii[0]),mul(by,local[1]/radii[1])),mul(bz,local[2]/radii[2])));
    this.vertex(add(center,world),normal);
   }
  }
  for(let r=0;r<rings;r++)for(let s=0;s<segments;s++){
   const a=base+r*(segments+1)+s,b=a+segments+1;
   if(r>0)this.indices.push(a,b,a+1);
   if(r<rings-1)this.indices.push(a+1,b,b+1);
  }
 }
 /** Closed tube along a polyline using rotation-minimising frames; end caps are fans. */
 tube(path,radius,sides=6,taper=null){
  if(path.length<2)return;
  const tangents=path.map((p,i)=>unit(sub(path[Math.min(i+1,path.length-1)],path[Math.max(i-1,0)])));
  let normal=cross(tangents[0],Math.abs(tangents[0][1])<.9?[0,1,0]:[1,0,0]);normal=unit(normal);
  const base=this.count;
  for(let i=0;i<path.length;i++){
   if(i>0){normal=sub(normal,mul(tangents[i],dot(normal,tangents[i])));normal=unit(normal);}
   const binormal=cross(tangents[i],normal),r=radius*(taper?taper(i/(path.length-1)):1);
   for(let s=0;s<sides;s++){
    const angle=2*Math.PI*s/sides,out=add(mul(normal,Math.cos(angle)),mul(binormal,Math.sin(angle)));
    this.vertex(add(path[i],mul(out,r)),out);
   }
  }
  for(let i=0;i<path.length-1;i++)for(let s=0;s<sides;s++){
   const a=base+i*sides+s,b=base+i*sides+(s+1)%sides,c=a+sides,d=b+sides;
   this.indices.push(a,c,b,b,c,d);
  }
  for(const [index,flip] of [[0,true],[path.length-1,false]]){
   const apex=this.vertex(path[index],mul(tangents[index],flip?-1:1)),ring=base+index*sides;
   for(let s=0;s<sides;s++){
    const a=ring+s,b=ring+(s+1)%sides;
    this.indices.push(...(flip?[apex,b,a]:[apex,a,b]));
   }
  }
 }
 /** Capped cylinder along an axis (proteasome barrel); the end rings are slightly wider. */
 barrel(center,axis,radius,height,sides=10){
  const a=unit(axis),half=mul(a,height/2),end0=sub(center,half),end1=add(center,half);
  this.tube([end0,add(end0,mul(a,height*.18)),sub(end1,mul(a,height*.18)),end1],radius,sides,t=>t<.15||t>.85?1.15:1);
 }
}

/** Spatial hash of every vertex of the source meshes, for clearance checks. */
class Obstacles{
 constructor(cell=.03){this.cell=cell;this.grid=new Map();}
 key(x,y,z){const c=this.cell;return `${Math.floor(x/c)},${Math.floor(y/c)},${Math.floor(z/c)}`;}
 addPoints(points){
  for(let i=0;i<points.length;i+=3){
   const key=this.key(points[i],points[i+1],points[i+2]);
   let bucket=this.grid.get(key);if(!bucket){bucket=[];this.grid.set(key,bucket);}
   bucket.push(points[i],points[i+1],points[i+2]);
  }
 }
 /** True when no stored vertex lies within radius of p. */
 clear(p,radius){
  const c=this.cell,span=Math.ceil(radius/c),x=Math.floor(p[0]/c),y=Math.floor(p[1]/c),z=Math.floor(p[2]/c),r2=radius*radius;
  for(let i=-span;i<=span;i++)for(let j=-span;j<=span;j++)for(let k=-span;k<=span;k++){
   const bucket=this.grid.get(`${x+i},${y+j},${z+k}`);
   if(!bucket)continue;
   for(let n=0;n<bucket.length;n+=3){
    const dx=bucket[n]-p[0],dy=bucket[n+1]-p[1],dz=bucket[n+2]-p[2];
    if(dx*dx+dy*dy+dz*dz<r2)return false;
   }
  }
  return true;
 }
}

/** Cytoplasm region of the cutaway: inside the plasma membrane slab, below the microvilli. */
const REGION={min:[-.56,.08,-.06],max:[.56,1.42,.36]};

export function buildOrganelles({obstacles:obstaclePoints,centrosome,golgi,nucleus}){
 const random=prng(20240917);
 const range=(lo,hi)=>lo+(hi-lo)*random();
 const obstacles=new Obstacles();
 for(const points of obstaclePoints)obstacles.addPoints(points);
 const inRegion=(p,margin=0)=>p.every((v,i)=>v>=REGION.min[i]+margin&&v<=REGION.max[i]-margin);
 const reserved=[];  // placed spheres [center,radius], so new items do not overlap each other
 const free=(p,radius,gap=.012,inside=null)=>(inside?inside(p,radius):inRegion(p,radius))&&obstacles.clear(p,radius+gap)&&reserved.every(([c,r])=>length(sub(c,p))>r+radius+gap);
 const randomPoint=(box=REGION)=>[range(box.min[0],box.max[0]),range(box.min[1],box.max[1]),range(box.min[2],box.max[2])];
 function place(radius,box,tries=4000,inside=null){
  for(let i=0;i<tries;i++){
   const p=randomPoint(box);
   if(free(p,radius,.012,inside)){reserved.push([p,radius]);return p;}
  }
  return null;
 }
 const randomBasis=()=>{
  const x=unit([range(-1,1),range(-1,1),range(-1,1)]),helper=Math.abs(x[1])<.9?[0,1,0]:[1,0,0],y=unit(cross(x,helper));
  return [x,y,cross(x,y)];
 };
 const parts=[];
 const stat={};
 const add_=(key,name,mesh,color,opacity=1)=>{if(mesh.count)parts.push({key,name,mesh,color,opacity});};

 // Peroxisomes: single-membrane bodies with a dense crystalloid matrix (TH H1.00.01.3.01074-01076).
 {
  const membrane=new Mesh(),matrix=new Mesh();let placed=0;
  for(let i=0;i<10;i++){
   const radius=range(.03,.042),p=place(radius*1.2);
   if(!p)continue;
   const basis=randomBasis(),radii=[radius*1.15,radius*.9,radius*.9];
   membrane.ellipsoid(p,radii,basis,14,10);
   matrix.ellipsoid(p,[radius*.55,radius*.5,radius*.5],basis,10,7);
   placed++;
  }
  stat.peroxisomes=placed;
  add_('peroxisomal-membrane','Peroxisomal membrane',membrane,[233,196,96],.55);
  add_('peroxisomal-matrix','Peroxisomal matrix',matrix,[176,104,48],1);
 }

 // Fat droplets (TH H1.00.01.3.01125-01126).
 {
  const mesh=new Mesh();let placed=0;
  for(let i=0;i<6;i++){
   const radius=range(.034,.056),p=place(radius);
   if(!p)continue;
   mesh.ellipsoid(p,[radius,radius,radius],null,16,10);placed++;
  }
  stat.fatDroplets=placed;
  add_('fat-droplet','Fat droplet',mesh,[247,226,140],.92);
 }

 // Glycogen granules, in clusters (TH H1.00.01.3.01121).
 {
  const mesh=new Mesh();let placed=0;
  for(let cluster=0;cluster<7;cluster++){
   const centre=place(.05,undefined,2000);
   if(!centre)continue;
   reserved.pop();
   for(let i=0;i<24;i++){
    const offset=mul(unit([range(-1,1),range(-1,1),range(-1,1)]),range(0,.05)),p=add(centre,offset),radius=range(.006,.0085);
    if(!free(p,radius,.004))continue;
    reserved.push([p,radius]);mesh.ellipsoid(p,[radius,radius,radius],null,7,5);placed++;
   }
  }
  stat.glycogenGranules=placed;
  add_('glycogen-granule','Glycogen granule',mesh,[176,126,205],1);
 }

 // Secretory vesicles gathered near the Golgi trans face (TH H1.00.01.3.01080).
 {
  const mesh=new Mesh(),box={min:sub(golgi,[.22,.28,.2]),max:add(golgi,[.22,.28,.2])};let placed=0;
  for(let i=0;i<16;i++){
   const radius=range(.016,.026),p=place(radius,box);
   if(!p)continue;
   mesh.ellipsoid(p,[radius,radius,radius],null,12,8);placed++;
  }
  stat.secretoryVesicles=placed;
  add_('secretory-vesicle','Secretory vesicle',mesh,[128,196,120],.85);
 }

 // Autophagosomes: double-membrane bodies around cargo (TH H1.00.01.3.01110).
 {
  const mesh=new Mesh();let placed=0;
  for(let i=0;i<2;i++){
   const radius=range(.05,.062),p=place(radius*1.1);
   if(!p)continue;
   const basis=randomBasis();
   mesh.ellipsoid(p,[radius,radius*.85,radius*.85],basis,16,10);
   mesh.ellipsoid(p,[radius*.78,radius*.64,radius*.64],basis,14,9);
   mesh.ellipsoid(p,[radius*.36,radius*.3,radius*.3],basis,10,7);placed++;
  }
  stat.autophagosomes=placed;
  add_('autophagosome','Autophagosome',mesh,[170,176,226],.7);
 }

 // Proteasomes: small capped barrels (TH H1.00.01.3.01113).
 {
  const mesh=new Mesh();let placed=0;
  for(let i=0;i<24;i++){
   const p=place(.02);
   if(!p)continue;
   reserved.pop();reserved.push([p,.017]);
   mesh.barrel(p,unit([range(-1,1),range(-1,1),range(-1,1)]),.008,.026,10);placed++;
  }
  stat.proteasomes=placed;
  add_('proteasome','Proteasome',mesh,[112,152,222],1);
 }

 // Walks that stay in free space: used by smooth ER and the filament systems.
 function walk(start,direction,steps,stepLength,clearance,wobble,bias){
  const path=[start];let heading=unit(direction),p=start;
  for(let i=0;i<steps;i++){
   let moved=false;
   for(let attempt=0;attempt<10&&!moved;attempt++){
    const wander=mul(unit([range(-1,1),range(-1,1),range(-1,1)]),wobble*(1+attempt*.35)),goal=bias?mul(unit(sub(bias(p),p)),.12):[0,0,0];
    const candidate=unit(add(add(heading,wander),goal)),q=add(p,mul(candidate,stepLength));
    if(inRegion(q,clearance)&&obstacles.clear(q,clearance)){heading=candidate;p=q;path.push(q);moved=true;}
   }
   if(!moved)break;
  }
  return path;
 }
 const smooth=path=>path.length<4?path:path.map((p,i)=>i===0||i===path.length-1?p:mul(add(add(path[i-1],mul(p,2)),path[i+1]),.25));

 // Smooth ER: branching tubules continuing from the rough ER around the nucleus (TH H1.00.01.3.01038).
 {
  const mesh=new Mesh();let tubes=0;
  const around={min:[-.5,.14,-.32],max:[.5,.95,.2]};
  for(let i=0;i<28&&tubes<22;i++){
   const start=place(.01,around,1500);
   if(!start)continue;
   reserved.pop();
   const path=smooth(walk(start,[range(-1,1),range(-1,1),range(-1,1)],26,.026,.016,.55));
   if(path.length<8)continue;
   mesh.tube(path,.0075,6);tubes++;
   for(const index of [Math.floor(path.length*.35),Math.floor(path.length*.7)]){
    const branch=smooth(walk(path[index],cross(sub(path[index+1],path[index-1]),[range(-1,1),range(-1,1),range(-1,1)]),10,.024,.016,.5));
    if(branch.length>=5)mesh.tube(branch,.006,6);
   }
  }
  stat.smoothErTubes=tubes;
  add_('smooth-er','Smooth endoplasmic reticulum',mesh,[140,205,186],.95);
 }

 // Microtubules: long, nearly straight filaments leaving the centrosome (TH H1.00.01.3.02002).
 {
  const mesh=new Mesh();let tubes=0;
  for(let i=0;i<60&&tubes<30;i++){
   const direction=unit([range(-1,1),range(-1,.4),range(-.5,.5)]);
   const path=smooth(walk(centrosome,direction,48,.026,.012,.12));
   if(path.length<14)continue;
   mesh.tube(path,.0045,5);tubes++;
  }
  stat.microtubules=tubes;
  add_('microtubule','Microtubule',mesh,[96,170,150],1);
 }

 // Intermediate filaments: wavy ropes ringing the nucleus and reaching the membrane (TH H1.00.01.3.02019).
 {
  const mesh=new Mesh();let tubes=0;
  const toward=nucleus;
  for(let i=0;i<40&&tubes<16;i++){
   const start=place(.008,undefined,500);
   if(!start)continue;
   reserved.pop();
   const path=smooth(walk(start,[range(-1,1),range(-1,1),range(-.3,.3)],34,.03,.011,.3,p=>{
    const away=sub(p,toward),distance=length(away);return distance<.3?add(toward,mul(unit(away),.34)):p;
   }));
   if(path.length<14)continue;
   mesh.tube(path,.0035,5);tubes++;
  }
  stat.intermediateFilaments=tubes;
  add_('intermediate-filament','Intermediate filament',mesh,[150,200,120],1);
 }

 // ---- Second set: vesicle traffic, residual material, polyribosomes and nuclear bodies ----
 const nucleusBox={min:sub(nucleus,[.17,.17,.17]),max:add(nucleus,[.17,.17,.17])};
 // inside the nuclear envelope, on the visible (front) side of the cutaway
 const insideNucleus=(p,r)=>length(sub(p,nucleus))<.165-r&&p[2]>nucleus[2]-.02;
 const lumpy=(mesh,centre,radius,lumps,basis)=>{
  mesh.ellipsoid(centre,[radius,radius*.85,radius*.85],basis,12,8);
  for(let i=0;i<lumps;i++){
   const o=mul(unit([range(-1,1),range(-1,1),range(-1,1)]),radius*.55),r=radius*range(.4,.6);
   mesh.ellipsoid(add(centre,o),[r,r,r],null,8,6);
  }
 };

 // Multivesicular bodies: a vesicle holding smaller internal vesicles (TH H1.00.01.3.01105).
 {
  const shell=new Mesh(),inner=new Mesh();let placed=0;
  for(let i=0;i<5;i++){
   const radius=range(.03,.042),p=place(radius*1.1);
   if(!p)continue;
   shell.ellipsoid(p,[radius,radius,radius],null,14,10);
   for(let k=0;k<4;k++){
    const o=mul(unit([range(-1,1),range(-1,1),range(-1,1)]),radius*.45),r=radius*.26;
    inner.ellipsoid(add(p,o),[r,r,r],null,8,6);
   }
   placed++;
  }
  stat.multivesicularBodies=placed;
  add_('multivesicular-body','Multivesicular body',shell,[176,166,224],.5);
  add_('internal-vesicle','Internal vesicle',inner,[132,120,200],1);
 }

 // Exosomes: very small released vesicles (TH H1.00.01.3.01107).
 {
  const mesh=new Mesh();let placed=0;
  for(let i=0;i<14;i++){
   const radius=range(.007,.011),p=place(radius);
   if(!p)continue;
   mesh.ellipsoid(p,[radius,radius,radius],null,8,6);placed++;
  }
  stat.exosomes=placed;
  add_('exosome','Exosome',mesh,[214,190,236],.95);
 }

 // Clathrin-coated vesicles: a vesicle inside an icosahedral clathrin cage (TH H1.00.01.3.01097).
 {
  const cage=new Mesh(),vesicle=new Mesh();let placed=0;
  const phi=(1+Math.sqrt(5))/2,verts=[[-1,phi,0],[1,phi,0],[-1,-phi,0],[1,-phi,0],[0,-1,phi],[0,1,phi],[0,-1,-phi],[0,1,-phi],[phi,0,-1],[phi,0,1],[-phi,0,-1],[-phi,0,1]].map(unit);
  const edges=[];for(let i=0;i<12;i++)for(let j=i+1;j<12;j++)if(length(sub(verts[i],verts[j]))<1.1)edges.push([i,j]);
  for(let i=0;i<9;i++){
   const radius=range(.02,.026),p=place(radius*1.25);
   if(!p)continue;
   const basis=randomBasis(),at=v=>add(p,add(add(mul(basis[0],v[0]*radius*1.2),mul(basis[1],v[1]*radius*1.2)),mul(basis[2],v[2]*radius*1.2)));
   vesicle.ellipsoid(p,[radius*.75,radius*.75,radius*.75],null,10,7);
   for(const [a,b] of edges)cage.tube([at(verts[a]),at(verts[b])],.0024,4);
   placed++;
  }
  stat.clathrinCoatedVesicles=placed;
  add_('clathrin-coated-vesicle','Clathrin-coated vesicle',vesicle,[236,180,120],1);
  add_('clathrin-coat','Clathrin coat',cage,[225,225,225],1);
 }

 // Transport vesicles between the rough ER and the Golgi (TH H1.00.01.3.01115).
 {
  const mesh=new Mesh(),box={min:[-.12,.72,0],max:[.36,1.08,.36]};let placed=0;
  for(let i=0;i<14;i++){
   const radius=range(.01,.015),p=place(radius,box);
   if(!p)continue;
   mesh.ellipsoid(p,[radius,radius,radius],null,10,7);placed++;
  }
  stat.transportVesicles=placed;
  add_('transport-vesicle','Transport vesicle',mesh,[120,190,210],.9);
 }

 // Endolysosomes: hybrid of endosome and lysosome with dense contents (TH H1.00.01.3.01104).
 {
  const shell=new Mesh(),dense=new Mesh();let placed=0;
  for(let i=0;i<3;i++){
   const radius=range(.03,.04),p=place(radius*1.1);
   if(!p)continue;
   shell.ellipsoid(p,[radius,radius,radius],null,14,10);
   lumpy(dense,p,radius*.6,3,randomBasis());placed++;
  }
  stat.endolysosomes=placed;
  add_('endolysosome','Endolysosome',shell,[190,150,210],.5);
  add_('endolysosome-contents','Endolysosome contents',dense,[110,76,140],1);
 }

 // Residual bodies: dense, irregular remnants of digestion (TH H1.00.01.3.01112).
 {
  const mesh=new Mesh();let placed=0;
  for(let i=0;i<3;i++){
   const radius=range(.03,.042),p=place(radius*1.3);
   if(!p)continue;
   lumpy(mesh,p,radius,5,randomBasis());placed++;
  }
  stat.residualBodies=placed;
  add_('residual-body','Residual body',mesh,[120,96,70],1);
 }

 // Pigment granules (lipofuscin-like) in small clusters (TH H1.00.01.3.01117).
 {
  const mesh=new Mesh();let placed=0;
  for(let cluster=0;cluster<5;cluster++){
   const centre=place(.04,undefined,2000);
   if(!centre)continue;
   reserved.pop();
   for(let i=0;i<6;i++){
    const p=add(centre,mul(unit([range(-1,1),range(-1,1),range(-1,1)]),range(0,.035))),radius=range(.007,.01);
    if(!free(p,radius,.004))continue;
    reserved.push([p,radius]);mesh.ellipsoid(p,[radius,radius,radius],null,8,6);placed++;
   }
  }
  stat.pigmentGranules=placed;
  add_('pigment-granule','Pigment granule',mesh,[96,70,40],1);
 }

 // Polyribosomes: ribosomes strung along a messenger RNA strand (TH H1.00.01.3.01045).
 {
  const beads=new Mesh(),strand=new Mesh();let chains=0;
  for(let i=0;i<30&&chains<9;i++){
   const start=place(.01,undefined,800);
   if(!start)continue;
   reserved.pop();
   const path=smooth(walk(start,[range(-1,1),range(-1,1),range(-1,1)],16,.02,.016,.5));
   if(path.length<9)continue;
   strand.tube(path,.0022,4);
   for(let k=1;k<path.length-1;k+=2){beads.ellipsoid(path[k],[.0085,.0085,.0085],null,8,6);}
   chains++;
  }
  stat.polyribosomes=chains;
  add_('polyribosome','Polyribosome',beads,[214,170,86],1);
  add_('messenger-rna-strand','Messenger RNA strand',strand,[226,226,226],1);
 }

 // Nuclear bodies on the visible side of the nucleus (TH H1.00.01.2.02014, 02015, 02017).
 {
  const coiled=new Mesh(),pml=new Mesh(),speckle=new Mesh();let a=0,b=0,c=0;
  for(let i=0;i<3;i++){
   const radius=range(.012,.017),p=place(radius*1.2,nucleusBox,4000,insideNucleus);
   if(!p)continue;
   lumpy(coiled,p,radius,3,randomBasis());a++;
  }
  for(let i=0;i<3;i++){
   const radius=range(.011,.015),p=place(radius*1.3,nucleusBox,4000,insideNucleus);
   if(!p)continue;
   pml.ellipsoid(p,[radius,radius,radius],null,12,8);pml.ellipsoid(p,[radius*.5,radius*.5,radius*.5],null,8,6);b++;
  }
  for(let i=0;i<6;i++){
   const centre=place(.02,nucleusBox,4000,insideNucleus);
   if(!centre)continue;
   reserved.pop();
   for(let k=0;k<7;k++){
    const p=add(centre,mul(unit([range(-1,1),range(-1,1),range(-1,1)]),range(0,.018))),r=range(.004,.0065);
    if(!free(p,r,.003,insideNucleus))continue;
    reserved.push([p,r]);speckle.ellipsoid(p,[r,r,r],null,7,5);c++;
   }
  }
  stat.nuclearBodies={coiled:a,pml:b,speckleGranules:c};
  add_('coiled-body','Coiled body',coiled,[236,120,120],1);
  add_('pml-nuclear-body','PML nuclear body',pml,[120,200,236],.85);
  add_('nuclear-speckle','Nuclear speckle',speckle,[240,220,100],1);
 }

 // ---- Plasma-membrane structures: desmosomes, hemidesmosomes and caveolae ----
 // The source membrane is a thin four-sided sheet: lateral walls near x = -0.62 / +0.60, base near y = 0.
 const WALL={left:-.598,right:.58},BASE=.012;
 const onWall=(p,r)=>p[1]>.12+r&&p[1]<1.38-r&&p[2]>.0+r&&p[2]<.29-r;
 const onBase=(p,r)=>Math.abs(p[0])<.52-r&&p[2]>.0+r&&p[2]<.29-r;
 {
  const junction=new Mesh(),plaque=new Mesh();let placed=0;
  for(let i=0;i<18;i++){
   const side=i%2?'right':'left',inward=side==='left'?1:-1;
   const probe=[WALL[side],range(.15,1.35),range(.03,.27)];
   if(!free(probe,.026,.004,onWall))continue;
   reserved.push([probe,.026]);
   junction.ellipsoid([WALL[side]-inward*.004,probe[1],probe[2]],[.0035,.03,.022],null,12,8);
   plaque.ellipsoid([WALL[side]+inward*.008,probe[1],probe[2]],[.006,.027,.019],null,12,8);
   placed++;
  }
  stat.desmosomes=placed;
  add_('desmosome','Desmosome',junction,[214,72,120],1);
  add_('desmosomal-plaque','Desmosomal plaque',plaque,[160,110,200],1);
 }
 {
  const mesh=new Mesh();let placed=0;
  for(let i=0;i<14;i++){
   const probe=[range(-.5,.5),BASE,range(.03,.27)];
   if(!free(probe,.026,.004,onBase))continue;
   reserved.push([probe,.026]);
   mesh.ellipsoid([probe[0],BASE+.004,probe[2]],[.03,.005,.02],null,12,8);placed++;
  }
  stat.hemidesmosomes=placed;
  add_('hemidesmosome','Hemidesmosome',mesh,[240,150,90],1);
 }
 {
  const mesh=new Mesh();let placed=0;
  for(let i=0;i<60&&placed<34;i++){
   const mode=random(),inwardAxis=mode<.4?[1,0,0]:mode<.8?[-1,0,0]:[0,1,0];
   const surface=mode<.4?[WALL.left,range(.15,1.35),range(.03,.27)]:mode<.8?[WALL.right,range(.15,1.35),range(.03,.27)]:[range(-.5,.5),BASE,range(.03,.27)];
   const body=add(surface,mul(inwardAxis,.018));
   if(!free(body,.012,.004,mode<.8?onWall:onBase))continue;
   reserved.push([body,.012]);
   mesh.ellipsoid(body,[.0105,.0105,.0105],null,10,7);
   mesh.tube([surface,add(surface,mul(inwardAxis,.01))],.0035,5);placed++;
  }
  stat.caveolae=placed;
  add_('caveola','Caveola',mesh,[120,170,230],1);
 }
 return {parts,stat};
}
