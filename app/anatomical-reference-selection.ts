import type {Atlas} from './anatomy';
import {liverReferenceChoices} from './liver-reference-selection';
import {tongueReferenceChoices} from './tongue-controls';
import {carotidReferenceChoices} from './carotid-surface-features';
export {liverReferenceActive as anatomicalReferenceActive,toggleLiverReferenceChoice as toggleAnatomicalReferenceChoice} from './liver-reference-selection';
const cache=new WeakMap<Atlas,ReturnType<typeof liverReferenceChoices>>();
export function anatomicalReferenceChoices(atlas:Atlas){
 let references=cache.get(atlas);if(references)return references;
 references=new Map([...liverReferenceChoices(atlas),...tongueReferenceChoices(atlas),...carotidReferenceChoices(atlas)]);cache.set(atlas,references);return references;
}
