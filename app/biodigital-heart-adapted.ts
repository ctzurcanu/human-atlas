import type {SourceCardiacLabelSite,SourceCardiacRegistry} from './biodigital-heart';

// Local review references are not accepted tissue. Their exact geometry and
// registration are fingerprinted by the reproducible offline preparation job.
// BEGIN REGISTERED CARDIAC DETAIL FINGERPRINTS
const approved:Record<string,string>={"A12.1.00.006":"2ba96dcb6f1748e3","A12.1.00.008":"11f9730c73452890","A12.1.00.021":"5e89c842d5be57fa","A12.1.00.024":"73ddb003d3efcd2b","A12.1.00.025":"6e8801a16d81a711","A12.1.00.026":"ff1073f9e3b35735","A12.1.00.028":"ecb3874b54b2fafb","A12.1.00.029":"d9175cee5528ca6a","A12.1.00.030":"31e9d527d363a767","A12.1.01.003":"e80b61d830dbda50","A12.1.01.004":"f0193528233cc96c","A12.1.01.005":"ecd62f2d3821e1b9","A12.1.01.006":"6669d26294a76596","A12.1.01.007":"619572a24c65dc7e","A12.1.01.008":"ba2190827cf08e76","A12.1.01.012":"0804e4c503fe9ed1","A12.1.01.013":"8a4276e4011a5b08","A12.1.01.014":"ef17e9266936d576","A12.1.02.004":"1e12d0f502845a25","A12.1.02.005":"a5440f0f715072bf","A12.1.02.006":"9fe72d879d999fa7","A12.1.02.007":"f6b66e67f12c17c3","A12.1.02.014":"ad2f4c8b7d5dbb33","A12.1.02.015":"510225ece968a4e8","A12.1.03.003":"d91a243ce5938674","A12.1.03.005":"780b9da58fce5ab5","A12.1.04.004":"cb8a4370f492d138","A12.1.04.005":"3d02fd8c6d229658","A12.1.04.006":"549673fe36479f12","A12.1.04.009":"045a737584b4f841","A12.1.04.011":"100051646e5d6a58","A12.1.04.016":"0ebb5dcb3132b727","A12.1.04.017":"2c086890ac99c064","A12.1.06.002":"aa3a7c64fd8aaa00","A12.1.06.003":"dbecdd4a8dd8e022","A12.1.06.004":"647fe67e5f0dc862","A12.1.06.005":"59f787ff32fe1eaf","A12.1.06.006":"bd5033ef16fcea2f","A12.1.06.007":"9a6749399780237d","A12.1.06.008":"1bf66628d8183b28","A12.1.07.001":"6d82eafea69e1d96","A12.1.07.002":"f61e6a5ef69bac52"};
// END REGISTERED CARDIAC DETAIL FINGERPRINTS
const canonical=(v:unknown):unknown=>Array.isArray(v)?v.map(canonical):v&&typeof v==='object'?Object.fromEntries(Object.entries(v).sort(([a],[b])=>a.localeCompare(b)).map(([k,x])=>[k,canonical(x)])):v;
export function cardiacDetailFingerprint(value:unknown){
 const text=JSON.stringify(canonical(value));let h1=0xdeadbeef,h2=0x41c6ce57;
 for(let i=0;i<text.length;i++){const ch=text.charCodeAt(i);h1=Math.imul(h1^ch,2654435761);h2=Math.imul(h2^ch,1597334677);}
 return (h1>>>0).toString(16).padStart(8,'0')+(h2>>>0).toString(16).padStart(8,'0');
}
const validated=new WeakSet<object>();
function freeze(value:unknown){if(value&&typeof value==='object'&&!Object.isFrozen(value)){Object.values(value).forEach(freeze);Object.freeze(value);}}
/** Cache only deeply frozen, already fingerprinted reference records. A changed
 * record must be a new object and is checked again; no stale mutable cache. */
export function sourceAdaptedCardiacReference(site:SourceCardiacLabelSite,m:SourceCardiacRegistry){
 const p=site as SourceCardiacLabelSite&{adaptation?:{schema?:string;status?:string;targetBankSha256?:string;anatomicalExtentVerified?:boolean}};
 if(!approved[p.term]||p.adaptation?.schema!=='registered-cardiac-detail/v1'||p.adaptation.status!=='provisional'||p.adaptation.targetBankSha256!==m.bankSha256||p.adaptation.anatomicalExtentVerified!==false||!p.coverageLimitation.startsWith('Provisional adapted ')||p.representation!==p.geometry?.type)return false;
 if(validated.has(p))return true;
 if(cardiacDetailFingerprint(p)!==approved[p.term])return false;
 freeze(p);validated.add(p);return true;
}
