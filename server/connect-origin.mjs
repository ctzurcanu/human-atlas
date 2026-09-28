// Browsers on our static viewer may use an HTTPS relay on another server.
// Keep the allowlist narrow; an invitation does not authorize arbitrary sites.
export const CONNECT_VIEWER_ORIGINS=['https://ctzurcanu.github.io'];
export function allowedConnectOrigin(request){
 try{
  const value=String(request.headers.origin),origin=new URL(value);
  if(value!==origin.origin)return false;
  const host=String(request.headers['x-forwarded-host']??request.headers.host??'').split(',')[0].trim();
  const secure=!!request.socket.encrypted||String(request.headers['x-forwarded-proto']??'').split(',')[0].trim()==='https';
  const loopback=['localhost','127.0.0.1','[::1]'].includes(origin.hostname);
  return secure&&origin.protocol==='https:'&&(origin.host===host||CONNECT_VIEWER_ORIGINS.includes(origin.origin))||!secure&&loopback&&origin.protocol==='http:'&&origin.host===host;
 }catch{return false;}
}
export function connectIdentityCors(request,response){
 if(!request.headers.origin)return true; // Same-origin browser GETs omit Origin.
 if(!allowedConnectOrigin(request))return false;
 response.setHeader('Access-Control-Allow-Origin',String(request.headers.origin));
 response.setHeader('Vary','Origin');return true;
}
