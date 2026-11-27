from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'app/scene.tsx';s=p.read_text()
def change(old,new):
 global s
 assert old in s,old[:120];s=s.replace(old,new,1)
s="import {selectedLiverBoundary} from './liver-completion-controls';\nimport {createLiverBoundaryOverlay} from './liver-boundary-overlay';\n"+s
change('const liverSurfaceRef=useRef(',"const liverBoundaryRef=useRef(selectedLiverBoundary(atlas,labelChoice));liverBoundaryRef.current=selectedLiverBoundary(atlas,labelChoice);\n const liverSurfaceRef=useRef(")
change('const surfaceFeatureOverlay=createSurfaceFeatureOverlay(desktopAnatomy,clipPlanes);','const surfaceFeatureOverlay=createSurfaceFeatureOverlay(desktopAnatomy,clipPlanes);\n  const liverBoundaryOverlay=createLiverBoundaryOverlay(desktopAnatomy,clipPlanes);let lastLiverBoundary=\'\',liverBoundaryIndex=-1;')
change("const surfaceFeature=liverSurfaceRef.current,featureKey=surfaceFeature?.conceptId??'';", "const boundary=liverBoundaryRef.current,boundaryKey=boundary?.conceptId??'';\n   if(boundaryKey!==lastLiverBoundary){lastLiverBoundary=boundaryKey;liverBoundaryIndex=boundary?atlas.parts.findIndex(p=>p.id===boundary.partIds[0]):-1;dirty=true;}\n   const surfaceFeature=liverSurfaceRef.current,featureKey=surfaceFeature?.conceptId??'';")
change('surfaceFeatureOverlay.update(surfaceFeature,',"liverBoundaryOverlay.update(boundary,liverBoundaryIndex>=0&&data[liverBoundaryIndex*4+3]>.5?pickers[liverBoundaryIndex]:undefined,liverBoundaryIndex>=0?atlas.parts[liverBoundaryIndex].geometryOrigin:undefined);\n    surfaceFeatureOverlay.update(surfaceFeature,")
change('surfaceFeatureOverlay.dispose();','surfaceFeatureOverlay.dispose();liverBoundaryOverlay.dispose();')
p.write_text(s);print('Boundary rendering wired')
