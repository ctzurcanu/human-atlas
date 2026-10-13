/** Move one shared interface towards its inward neighbor without creating
 * per-part geometry. The scalar field must have passed the geometry gates. */
export function relativeTissueInterface(outer:Float32Array,inner:Float32Array,weights:Float32Array) {
  if(outer.length!==inner.length||outer.length!==weights.length*3
    ||weights.some(w=>!Number.isFinite(w)||w<0||w>.5))throw Error('Invalid relative tissue field');
  return Float32Array.from(outer,(v,i)=>v+(inner[i]-v)*weights[Math.floor(i/3)]);
}
