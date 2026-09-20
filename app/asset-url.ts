/// <reference types="vite/client" />

/** Resolve public assets beneath the deployment path, including catalogue chunk URLs. */
export const assetUrl=(path:string)=>`${import.meta.env.BASE_URL}${path.replace(/^\/+/,'')}`;
