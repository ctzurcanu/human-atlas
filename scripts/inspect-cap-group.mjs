/** Inspect one publisher concept's axial cut without changing its geometry. */
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
import * as T from 'three';
import {mergeGeometries} from 'three/addons/utils/BufferGeometryUtils.js';

const [concept='REF:Rectus abdominis muscle.l',height='1.05257']=process.argv.slice(2);
const atlas=JSON.parse(await readFile('.local-models/reference.json','utf8'));
const parts=atlas.parts.filter(part=>part.conceptId===concept&&part.bounds[0][1]<=Number(height)&&part.bounds[1][1]>=Number(height));
if(!parts.length)throw Error(`No primitives intersect ${height} for ${concept}`);
async function load(entry){const output=await build({entryPoints:[entry],bundle:true,platform:'node',format:'esm',write:false});return import(`data:text/javascript;base64,${Buffer.from(output.outputFiles[0].text).toString('base64')}`);}
const [{sectionCapGeometry},{sectionCapProfile},{analyzeMeshTopology}]=await Promise.all([load('app/section-cap.ts'),load('app/section-cap-profile.ts'),load('app/mesh-topology.ts')]);
const cache=new Map();
async function geometry(part){let chunk=cache.get(part.chunk);if(!chunk){const url=atlas.chunks[part.chunk].url;chunk=await readFile(url.startsWith('/local-models/')?`.local-models/${url.slice('/local-models/'.length)}`:`public${url}`);cache.set(part.chunk,chunk);}
 const result=new T.BufferGeometry();result.setAttribute('position',new T.BufferAttribute(new Float32Array(chunk.buffer,chunk.byteOffset+part.positions,part.vertexCount*3),3));result.setIndex(new T.BufferAttribute(new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount),1));return result;}
const plane=new T.Plane(new T.Vector3(0,-1,0),Number(height)),fragments=await Promise.all(parts.map(geometry));
function inspect(part,source){const topology=analyzeMeshTopology(source),profile=sectionCapProfile(part,topology),cap=sectionCapGeometry(source,plane,new T.Matrix4(),profile);return {name:part.name,triangles:topology.triangles,boundaryRatio:Number(topology.boundaryRatio.toFixed(4)),profile,cap:cap?.userData??null};}
console.log(JSON.stringify({concept,height:Number(height),individual:fragments.map((g,i)=>inspect(parts[i],g)),combined:inspect(parts[0],mergeGeometries(fragments)??fragments[0])},null,2));
