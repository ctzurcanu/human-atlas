import * as T from 'three';
import {fetchAsset} from './asset-cache';

type Coordinates = [number, number, number];
interface Volume {
  url: string;
  min: Coordinates;
  max: Coordinates;
  size: Coordinates;
}

export interface FingerSurfacePreview {
  skinId: string;
  affectedParts: string[];
  volumes: [Volume, Volume];
}

// This local review overlay treats the measured skin as the boundary of the
// assembled view. The source meshes remain available for dissection.
export function createFingerSurfacePreview(
  config: FingerSurfacePreview | undefined,
  indices: Map<string, number>,
  width: number,
  signal: AbortSignal,
  requestRender: () => void,
) {
  if (!config || !import.meta.env.DEV || config.volumes.length !== 2 ||
    config.volumes.some(volume =>
      !/^\/local-models\/female-ta98-finger-surface-preview\/[a-z]+-[0-9a-f]{20}\.bin$/.test(volume.url) ||
      volume.size.length !== 3 || volume.size.some(n => !Number.isInteger(n) || n < 2 || n > 256) ||
      !volume.min.every((v, i) => Number.isFinite(v) && Number.isFinite(volume.max[i]) && v < volume.max[i])
    )) return;
  const skinIndex = indices.get(config.skinId);
  if (skinIndex === undefined) return;
  const mask = new Uint8Array(width);
  for (const id of config.affectedParts) {
    const index = indices.get(id);
    if (index !== undefined && index !== skinIndex) mask[index] = 255;
  }
  const parts = new T.DataTexture(mask, width, 1, T.RedFormat);
  parts.needsUpdate = true;
  const enabled = {value: 0};
  const volumes = config.volumes.map(volume => {
    const texture = new T.Data3DTexture(new Uint8Array([255]), 1, 1, 1);
    texture.format = T.RedFormat;
    texture.minFilter = texture.magFilter = T.LinearFilter;
    texture.needsUpdate = true;
    return {
      texture,
      min: new T.Vector3(...volume.min),
      span: new T.Vector3(...volume.max).sub(new T.Vector3(...volume.min)),
    };
  });
  let ready = false;
  void Promise.all(config.volumes.map(async (volume, i) => {
    const response = await fetchAsset(volume.url, {signal});
    if (!response.ok) throw new Error('Finger preview field unavailable');
    const bytes = new Uint8Array(await response.arrayBuffer());
    const [width, height, depth] = volume.size;
    if (bytes.length !== width * height * depth) throw new Error('Invalid finger preview field size');
    if (signal.aborted) return;
    volumes[i].texture.image = {data: bytes, width, height, depth};
    volumes[i].texture.needsUpdate = true;
  })).then(() => {
    if (signal.aborted) return;
    ready = true;
    requestRender();
  }).catch(error => {
    if (!signal.aborted) console.error(error);
  });
  signal.addEventListener('abort', () => {
    parts.dispose();
    for (const volume of volumes) volume.texture.dispose();
  }, {once: true});
  return {
    skinIndex,
    update(skinOpaque: boolean) {
      const next = ready && skinOpaque ? 1 : 0;
      if (next !== enabled.value) {
        enabled.value = next;
        requestRender();
      }
      return !ready ? 'loading' : next ? 'active' : 'inactive';
    },
    apply(shader: Parameters<T.Material['onBeforeCompile']>[0]) {
      shader.uniforms.fingerPreviewEnabled = enabled;
      shader.uniforms.fingerPreviewParts = {value: parts};
      for (const [i, volume] of volumes.entries()) {
        shader.uniforms[`fingerPreviewVolume${i}`] = {value: volume.texture};
        shader.uniforms[`fingerPreviewMin${i}`] = {value: volume.min};
        shader.uniforms[`fingerPreviewSpan${i}`] = {value: volume.span};
      }
      shader.vertexShader = 'varying vec3 fingerPreviewPoint; varying float fingerPreviewPart;\n' + shader.vertexShader;
      shader.vertexShader = shader.vertexShader.replace('#include <project_vertex>', `
        fingerPreviewPoint = (modelMatrix * vec4(transformed, 1.0)).xyz;
        fingerPreviewPart = partIndex;
        #include <project_vertex>
      `);
      shader.fragmentShader = `
        uniform float fingerPreviewEnabled;
        uniform sampler2D fingerPreviewParts;
        uniform highp sampler3D fingerPreviewVolume0;
        uniform highp sampler3D fingerPreviewVolume1;
        uniform vec3 fingerPreviewMin0, fingerPreviewMin1, fingerPreviewSpan0, fingerPreviewSpan1;
        uniform float stateWidth;
        varying vec3 fingerPreviewPoint;
        varying float fingerPreviewPart;
        bool withinFingerPreview(vec3 uvw) {
          return all(greaterThanEqual(uvw, vec3(0.0))) && all(lessThanEqual(uvw, vec3(1.0)));
        }
      ` + shader.fragmentShader;
      shader.fragmentShader = shader.fragmentShader.replace('#include <clipping_planes_fragment>', `
        #include <clipping_planes_fragment>
        if (fingerPreviewEnabled > 0.5 && texture2D(fingerPreviewParts, vec2((fingerPreviewPart + 0.5) / stateWidth, 0.5)).r > 0.5) {
          vec3 uvw0 = (fingerPreviewPoint - fingerPreviewMin0) / fingerPreviewSpan0;
          vec3 uvw1 = (fingerPreviewPoint - fingerPreviewMin1) / fingerPreviewSpan1;
          if (withinFingerPreview(uvw0) && texture(fingerPreviewVolume0, uvw0).r < 0.60) discard;
          if (withinFingerPreview(uvw1) && texture(fingerPreviewVolume1, uvw1).r < 0.60) discard;
        }
      `);
    },
  };
}
