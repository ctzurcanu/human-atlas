/** A lower-opacity neighbor must not erase a more visible tissue boundary.
 * Equal-opacity assembled tissues omit their internal shared interface. */
export function tissueNeighborOccludes(sourceOpacity:number,neighborOpacity:number):boolean{
 return Number.isFinite(sourceOpacity)&&Number.isFinite(neighborOpacity)&&sourceOpacity>=.001&&neighborOpacity>=.001&&neighborOpacity+1e-6>=sourceOpacity;
}
