"""Bound existing read-only brain inventories into a practical source decision."""
import argparse,hashlib,json,time
from pathlib import Path
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,required=True);args=ap.parse_args();out=args.directory.resolve();began=time.monotonic()
 destination=out/'comparison.json';assert not destination.exists()
 npath=out/'native-brain.json';dpath=out/'biodigital-brain.json';n=json.loads(npath.read_text());d=json.loads(dpath.read_text())
 for inventory in [n,d]:
  for path,digest in inventory['inputs'].items():assert sha(path)==digest,path
 assert n['nativeBrainNeuralTissue']['placedPartFragments']==sum(n['categories'][k]['placedPartFragments'] for k in ['brainstem-neural-tissue','cerebellar-neural-tissue','cortical-neural-surface-region','deep-cerebral-diencephalic-white-matter-neural-tissue'])
 assert n['nativeBrainNeuralTissue']['uniqueSurfaceTriangleTotal']==sum(n['categories'][k]['uniqueSurfaceTriangleTotal'] for k in ['brainstem-neural-tissue','cerebellar-neural-tissue','cortical-neural-surface-region','deep-cerebral-diencephalic-white-matter-neural-tissue'])
 def native(k):
  r=n['categories'][k];return dict(namedRegionsWithLaterality=r['genuineNamedSourceRegionsWithLaterality'],unsidedNamedTypes=r['unsidedSourceRegionNames'],physicalMaterialFragments=r['placedPartFragments'],uniqueExactWorldSurfaceMeshes=r['uniqueExactWorldSurfaceMeshes'],placedTriangles=r['placedFragmentTriangleTotal'])
 def donor(k):
  r=d['neuralRegionFamilies'][k];return dict(namedRegionsWithLaterality=r['physicalNamedInstances'],unsidedNamedTypes=r['distinctUnsidedNames'],uniquePrimarySharedAssets=r['uniquePrimaryGeometryAssets'],uniqueDecodedContentAssets=r['uniqueExactDecodedAssetContents'],uniqueAssetAndMatrixPlacements=r['uniqueExactAssetAndMatrixPlacements'],placedTriangles=r['instanceTriangles'],sharedAssetTriangles=r['uniqueAssetTriangles'])
 comparisons=[dict(region='Cerebral cortex',native=native('cortical-neural-surface-region'),donor=donor('cerebral-cortex'),assessment='Native has more named cortical subdivisions at a similar total triangle count; counts do not establish the better visible cortical contour.'),
  dict(region='Cerebellum',native=native('cerebellar-neural-tissue'),donor=donor('cerebellum'),assessment='Native includes35separate named neural meshes including lobular/vermal subdivisions; donor has one coarse whole-cerebellum mesh.'),
  dict(region='Brainstem',native=native('brainstem-neural-tissue'),donor=donor('brainstem'),assessment='Native has substantially more distinct named brainstem geometry and tessellation; donor still has separate midbrain/pons/medulla and collicular/peduncular meshes.')]
 result=dict(schema='native-versus-biodigital-brain-source-comparison/v1',status='complete-read-only-source-decision',recommendation='retain-current-native-brain',
  rationale='The donor is not clearly better. The current native source has more genuine named brain regions and far richer cerebellar, brainstem and deep-brain geometry. Preserve it while eventual peripheral nervous or lymphatic sources are considered separately.',
  evidenceInputs={str(npath):sha(npath),str(dpath):sha(dpath),str(out/'future-import-reuse-notes.json'):sha(out/'future-import-reuse-notes.json')},
  nativeActiveBrain=n['nativeBrainNeuralTissue'],donorBrainNeuralTissue=d['brainNeuralTissue'],comparisons=comparisons,
  nativeSeparateNonParenchymalCategories={k:v for k,v in n['categories'].items() if k not in ['brainstem-neural-tissue','cerebellar-neural-tissue','cortical-neural-surface-region','deep-cerebral-diencephalic-white-matter-neural-tissue']},
  donorSeparateNonParenchymalCategories={k:v for k,v in d['categories'].items() if k!='brain-neural-tissue'},
  countingPolicy='Named source regions are distinct from nongeometric concepts and from native material fragments. Source shared assets are counted once; bilateral placed tissue cost is reported separately. Native pre-bake sharing is not recoverable from its packed world-coordinate catalogue, so native unique world meshes are not equated with original donor asset count.',
  sourceScopeWarning='Original donor Brain organizer includes40eye/lacrimal meshes,12cranial nerve meshes,7cavity/boundary meshes,3associated endocrine meshes and2choroid-plexus meshes. These are separately counted and are not brain-neural detail. Septal partition is separate for both sources.',
  qualityQualification='No rendered surface, cross-section, continuity or self-intersection comparison was performed. Named subdivisions and actual triangle density establish resolution evidence, not anatomical acceptance or superior surface quality. Donor primary source assets contain24exactzero-area neural triangles and no decoder-rejected triangles.',
  reusedPipelineNotes=str(out/'future-import-reuse-notes.json'),brainOrNervousOrLymphaticImportPerformed=False,brainDeformationPerformed=False,profileOrDefaultOrSourceGeometryChanged=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
 destination.write_text(json.dumps(result,indent=2)+'\n');(out/'checkpoint.json').write_text(json.dumps(dict(status=result['status'],recommendation=result['recommendation'],comparisonSha256=sha(destination),elapsedSeconds=result['elapsedSeconds']))+'\n')
 print(json.dumps(dict(output=str(destination),sha256=sha(destination),recommendation=result['recommendation'],native=n['nativeBrainNeuralTissue'],donor=d['brainNeuralTissue'],elapsedSeconds=result['elapsedSeconds'])))
if __name__=='__main__':main()
