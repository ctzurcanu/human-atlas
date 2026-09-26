// The publisher's 27 spinal cord tissue meshes are straight extrusions
// from y=1.348 to 1.853. That places their upper half above the neck, while
// the independently modeled spinal dura follows the vertebral canal from
// approximately y=1.055 to 1.545. Bend and fit these tissues into that canal.
export const spinalCordTissues=new Set([
 'Anterior corticospinal tract',
 'Anterior fasciculus proprius',
 'White matter of spinal cord',
 'Anterior horn of spinal cord',
 'Posterior horn of spinal cord',
 'Anterior spinocerebellar tract',
 'Anterior spinothalamic tract',
 'Central canal',
 'Gracile fasciculus',
 'Intermediolateral nucleus',
 'Intermediomedial nucleus',
 'Lateral corticospinal tract',
 'Lateral fasciculus proprius',
 'Lateral intermediate substance',
 'Lateral reticulospinal tract',
 'Lateral spinothalamic tract',
 'Lateral vestibulospinal tract',
 'Medial reticulospinal tract',
 'Medial vestibulospinal tract',
 'Nucleus proprius',
 'Posterior fasciculus proprius',
 'Posterior spinocerebellar tract',
 'Posterolateral tract',
 'Rubrospinal tract',
 'Spinal reticular process',
 'Spinotectal tract',
 'Tectospinal tract',
]);
const originalMin=1.3475721642932;
const targetMin=1.105,targetMax=1.525;
const yScale=(targetMax-targetMin)/.505;
const centerline=[
 [1.105,-.031],[1.125,-.0365],[1.15,-.044],[1.175,-.052],
 [1.20,-.0595],[1.225,-.064],[1.25,-.068],[1.275,-.071],
 [1.30,-.0725],[1.325,-.0715],[1.35,-.069],[1.375,-.0635],
 [1.40,-.055],[1.425,-.0425],[1.45,-.032],[1.475,-.0245],
 [1.50,-.0215],[1.525,-.0265],
];
function canalAt(y){
 let index=centerline.length-2;
 for(let i=0;i<centerline.length-1;i++)if(y<=centerline[i+1][0]){index=i;break;}
 const [y0,z0]=centerline[index],[y1,z1]=centerline[index+1];
 const t=Math.max(0,Math.min(1,(y-y0)/(y1-y0)));
 const before=centerline[Math.max(0,index-1)],after=centerline[Math.min(centerline.length-1,index+2)];
 const m0=(z1-before[1])/(y1-before[0]),m1=(after[1]-z0)/(after[0]-y0);
 const t2=t*t,t3=t2*t,h=y1-y0;
 const z=(2*t3-3*t2+1)*z0+(t3-2*t2+t)*h*m0+(-2*t3+3*t2)*z1+(t3-t2)*h*m1;
 const slope=((6*t2-6*t)*z0+(3*t2-4*t+1)*h*m0+(-6*t2+6*t)*z1+(3*t2-2*t)*h*m1)/h;
 return {z,slope};
}
export function placeSpinalCord(point,normal){
 const y=Math.max(targetMin,Math.min(targetMax,targetMin+(point.y-originalMin)*yScale));
 const {z,slope}=canalAt(y),oldZ=point.z;
 const xScale=y<1.13?.45+(y-targetMin)/(.025)*.15:.6,zScale=.4;
 point.set(point.x*xScale,y,z+(oldZ+.0459)*zScale);
 normal.set(normal.x/xScale,normal.y/yScale-slope*normal.z/zScale,normal.z/zScale).normalize();
}
