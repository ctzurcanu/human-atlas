import type {Atlas, Part} from './anatomy';

type Offsets = Pick<Part, 'chunk'|'positions'|'normals'|'indices'|'uvs'|'colors'|'normalSigns'|'tissueChart'|'gastricLumen'>;
export type SystemGeometry = {chunks:Atlas['chunks']; parts:Map<string, Part>};
type PackagedAtlas = Atlas & {sourceGeometry?:{chunks:Atlas['chunks']; parts:Record<string,Offsets>}};

/** Verified source poses identify anatomy. The renderer reads the same bytes
 * from the system banks, whose offsets can differ from the original poses. */
export function withSystemGeometry(atlas:PackagedAtlas):Atlas {
 const {sourceGeometry, ...model} = atlas;
 if (!sourceGeometry) return atlas;
 const parts = model.parts.map(part => {
  const pose = sourceGeometry.parts[part.id];
  if (!pose || !sourceGeometry.chunks[pose.chunk] || !model.chunks[part.chunk]) throw Error('Incomplete system geometry for ' + part.id);
  return {...part, ...pose};
 });
 return {...model, parts, chunks:sourceGeometry.chunks,
  systemGeometry:{chunks:model.chunks, parts:new Map(model.parts.map(p => [p.id, p]))}};
}
