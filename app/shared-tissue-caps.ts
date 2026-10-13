import * as T from 'three';
import {tissueWallMaterial, type TissueAppearance} from './tissue-wall-material';

/** Live stencil caps for closed tissue shells. Updating a plane or peel matrix
 * changes uniforms/transforms only: no contour jobs and no release debounce.
 * The counting surface includes inner faces and annular rims, so the lumen
 * subtracts itself. Never use this with a single unpaired outer surface.
 */
export function createSharedTissueCaps(scene: T.Scene, shells: readonly T.BufferGeometry[],
  appearances: readonly TissueAppearance[], plane: T.Plane) {
  if (shells.length!==appearances.length) throw new Error('Missing tissue appearance');
  const planeGeometry=new T.PlaneGeometry(1,1);
  const counters: T.Mesh[] = [], caps: T.Mesh[] = [], materials: T.Material[] = [];
  const frames: T.Matrix4[] = [];
  for(let k=0;k<shells.length;k++) {
    if(!shells[k].userData.pairedTissueWall) throw new Error('Caps require a validated paired tissue shell');
    for(const side of [T.BackSide,T.FrontSide]) {
      const m=new T.MeshBasicMaterial({side,colorWrite:false,depthWrite:false,depthTest:false,clippingPlanes:[plane]});
      m.stencilWrite=true;m.stencilFunc=T.AlwaysStencilFunc;
      m.stencilFail=m.stencilZFail=m.stencilZPass=side===T.BackSide?T.IncrementWrapStencilOp:T.DecrementWrapStencilOp;
      const mesh=new T.Mesh(shells[k],m);mesh.matrixAutoUpdate=false;mesh.frustumCulled=false;
      mesh.visible=false;mesh.renderOrder=-1000+k*3+(side===T.BackSide?0:1);scene.add(mesh);counters.push(mesh);materials.push(m);
    }
    const frame=new T.Matrix4(),{material}=tissueWallMaterial(appearances[k],frame);frames.push(frame);
    material.side=T.DoubleSide;material.stencilWrite=true;material.stencilFunc=T.NotEqualStencilFunc;
    material.stencilRef=0;material.stencilFail=material.stencilZFail=material.stencilZPass=T.ZeroStencilOp;
    const cap=new T.Mesh(planeGeometry,material);cap.visible=false;cap.frustumCulled=false;
    cap.renderOrder=-998+k*3;scene.add(cap);caps.push(cap);materials.push(material);
  }
  const sphere=new T.Sphere(),center=new T.Vector3(),scale=new T.Vector3();
  const forward=new T.Vector3(0,0,1),facing=new T.Vector3();
  return {
    caps,counters,
    update(visible: readonly boolean[], matrices: readonly T.Matrix4[], enabled: boolean) {
      if(visible.length!==shells.length || matrices.length!==shells.length) throw new Error('Missing tissue cap state');
      for(let k=0;k<shells.length;k++) {
        const show=enabled&&visible[k];caps[k].visible=show;
        for(let j=k*2;j<k*2+2;j++){const counter=counters[j];counter.visible=show;counter.matrix.copy(matrices[k]);counter.matrixWorldNeedsUpdate=true;}
        if(!show)continue;
        frames[k].copy(matrices[k]).invert();
        const bounds=shells[k].boundingSphere;if(!bounds)throw new Error('Tissue shell needs bounds');
        sphere.copy(bounds).applyMatrix4(matrices[k]);plane.projectPoint(sphere.center,center);
        caps[k].position.copy(center);caps[k].quaternion.setFromUnitVectors(forward,facing.copy(plane.normal).negate());
        scale.setScalar(sphere.radius*2.01);scale.z=1;caps[k].scale.copy(scale);
      }
    },
    dispose() {for(const mesh of [...counters,...caps])scene.remove(mesh);for(const material of materials)material.dispose();planeGeometry.dispose();},
  };
}
