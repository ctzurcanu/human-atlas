// User-chosen embryo view, in the same normalized format as shared URLs.
const EMBRYO_FRAME=[.944333,.274867,.543392,.422751,.446347,.784766,.5,.450126,.369345,.054656,.724182,.462439];
export function defaultModelCamera(model:string){
 return model==='embryo'?[...EMBRYO_FRAME]:undefined;
}
