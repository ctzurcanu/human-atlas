/** Keep exact atlas selection IDs in the smallest common shader attribute type. */
export function partIndexBuffer(partIndex:number,vertexCount:number,totalParts:number):Uint16Array|Float32Array {
 if(!Number.isInteger(partIndex)||partIndex<0||partIndex>=totalParts)throw new RangeError('Invalid atlas part index');
 if(!Number.isInteger(vertexCount)||vertexCount<0)throw new RangeError('Invalid anatomy vertex count');
 // Every geometry in an atlas must use the same type for Three's merged batches.
 // Unsigned-short attributes are converted to exact floats by the float shader input.
 const values=totalParts<=65536?new Uint16Array(vertexCount):new Float32Array(vertexCount);
 return values.fill(partIndex);
}
