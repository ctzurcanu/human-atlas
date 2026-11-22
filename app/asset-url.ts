/// <reference types="vite/client" />

/** Resolve public assets beneath the deployment path, including catalogue chunk URLs. */
export function assetUrl(path:string) {
 const file=path.replace(/^\/+/, '');
 const repositoryAssets=import.meta.env.VITE_FEMALE_TA98_ASSETS;
 if(repositoryAssets && file==='models/female-ta98/atlas.json')return import.meta.env.BASE_URL + file + '.gz';
 if(repositoryAssets && /^models\/female-ta98\/(geometry|textures)\//.test(file)) {
  return repositoryAssets + file;
 }
 return import.meta.env.BASE_URL + file;
}
