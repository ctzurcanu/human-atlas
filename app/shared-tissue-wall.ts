import {BufferAttribute, BufferGeometry, Box3, Vector3, Sphere} from 'three';

export interface SharedWallInput {
  /** Ordered interfaces of one repaired, paired wall, in organ-local metres. */
  interfaces: readonly Float32Array[];
  faces: Uint32Array;
  /** Directed edges of the reference surface's real openings. */
  rims: Uint32Array;
  layerIds: readonly string[];
}

/** One position/normal pool, independent tissue draw lists, no texture storage.
 * Coordinates are immutable. Moving a peeled shell changes its object matrix,
 * never its interface positions. The caller must supply validated anatomy.
 */
export function createSharedTissueWall(input: SharedWallInput) {
  const {interfaces, faces, rims, layerIds} = input;
  const count = interfaces[0]?.length / 3;
  if (!Number.isInteger(count) || count < 3 || interfaces.length !== layerIds.length + 1
    || new Set(layerIds).size !== layerIds.length || faces.length % 3 || rims.length % 2
    || interfaces.some(p => p.length !== count * 3 || p.some(x => !Number.isFinite(x)))
    || faces.some(i => i >= count) || rims.some(i => i >= count)) throw new Error('Invalid shared tissue wall');
  // One outward normal field per shared interface. Material slot 1 reverses
  // that field for an inner face; slot 2 shades annular rims flat. Keeping the
  // orientation in the material avoids duplicating all position/normal data.
  const bank = count * interfaces.length;
  const values = new Float32Array(bank * 3), normals = new Int16Array(bank * 3);
  for (let k = 0; k < interfaces.length; k++) {
    const p = interfaces[k], sums = new Float64Array(p.length);
    values.set(p, k * p.length);
    for (let j = 0; j < faces.length; j += 3) {
      const a = faces[j] * 3, b = faces[j+1] * 3, c = faces[j+2] * 3;
      const ux = p[b]-p[a], uy = p[b+1]-p[a+1], uz = p[b+2]-p[a+2];
      const vx = p[c]-p[a], vy = p[c+1]-p[a+1], vz = p[c+2]-p[a+2];
      const n = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx];
      for (const id of [a,b,c]) for (let axis = 0; axis < 3; axis++) sums[id+axis] += n[axis];
    }
    for (let j = 0; j < count; j++) {
      const length = Math.hypot(sums[j*3], sums[j*3+1], sums[j*3+2]);
      if (!length) throw new Error('Undefined interface normal');
      for (let axis = 0; axis < 3; axis++) {
        const n = Math.round(sums[j*3+axis]/length*32767), id = k*count*3+j*3+axis;
        normals[id] = n;
      }
    }
  }
  const position = new BufferAttribute(values,3), normal = new BufferAttribute(normals,3,true);
  const full: BufferGeometry[] = [], displayed: BufferGeometry[] = [];
  const Index = bank <= 65536 ? Uint16Array : Uint32Array;
  const fill = (k: number, outer: boolean, inner: boolean) => {
    const indices = new Index((Number(outer)+Number(inner))*faces.length+rims.length*3);
    let cursor = 0;
    if (outer) for (const f of faces) indices[cursor++] = k*count+f;
    if (inner) for (let j = 0; j < faces.length; j += 3)
      for (const f of [faces[j],faces[j+2],faces[j+1]]) indices[cursor++] = (k+1)*count+f;
    const surfaceCount = cursor;
    for (let j = 0; j < rims.length; j += 2) {
      const a = k*count+rims[j], b = k*count+rims[j+1];
      const A = (k+1)*count+rims[j], B = (k+1)*count+rims[j+1];
      for (const id of [b,a,A,b,A,B]) indices[cursor++] = id;
    }
    return {indices,surfaceCount};
  };
  const setDraw = (g: BufferGeometry,outer: boolean,inner: boolean) => {
    g.clearGroups();
    if (outer) g.addGroup(0,faces.length,0);
    if (inner) g.addGroup(faces.length,faces.length,1);
    g.addGroup(faces.length*2,rims.length*3,2); // flat annular rim material
  };
  for (let k = 0; k < layerIds.length; k++) {
    const box = new Box3(), point = new Vector3();
    for (const p of [interfaces[k],interfaces[k+1]]) for (let j = 0; j < p.length; j += 3)
      box.expandByPoint(point.fromArray(p,j));
    const make = () => {
      const g = new BufferGeometry(); g.setAttribute('position',position); g.setAttribute('normal',normal);
      g.boundingBox = box.clone(); g.boundingSphere = box.getBoundingSphere(new Sphere());
      g.userData.tissueLayer = layerIds[k]; g.userData.pairedTissueWall = true;
      return g;
    };
    const closed = make(); closed.setIndex(new BufferAttribute(fill(k,true,true).indices,1));
    setDraw(closed,true,true); full.push(closed);
    const view = make(); view.setIndex(closed.index); setDraw(view,true,true); displayed.push(view);
  }
  let signature = '';
  return {
    full, displayed, position, normal,
    /** Hide internal coincident faces only where adjacent shells still touch.
     * A separated junction exposes both sides. No camera or scene mutation.
     */
    update(visible: readonly boolean[], separated: readonly boolean[] = []) {
      if (visible.length !== layerIds.length || separated.length && separated.length !== layerIds.length-1)
        throw new Error('Invalid tissue visibility');
      const key = visible.map(Number).join('')+':'+separated.map(Number).join('');
      if (key === signature) return false; signature = key;
      for (let k = 0; k < layerIds.length; k++) {
        if (!visible[k]) { displayed[k].setDrawRange(0,0); continue; }
        setDraw(displayed[k],k===0 || !visible[k-1] || !!separated[k-1],
          k===layerIds.length-1 || !visible[k+1] || !!separated[k]);
        displayed[k].setDrawRange(0,Infinity);
      }
      return true;
    },
    bytes: {position: values.byteLength, normal: normals.byteLength,
      closedIndices: full.reduce((sum,g) => sum+g.index!.array.byteLength,0)},
    dispose() { for (const g of [...full,...displayed]) g.dispose(); },
  };
}
