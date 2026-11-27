import {mkdir,readFile,writeFile} from 'node:fs/promises';
import {resolve} from 'node:path';

async function writeWhenChanged(file,content){
 let existing;
 try{existing=await readFile(file,'utf8');}catch(error){if(error.code!=='ENOENT')throw error;}
 if(existing!==content)await writeFile(file,content);
}

/** Generate Vite's ignored web entries from the tracked comparison sources. */
export async function prepareComparisonPage(projectRoot){
 const root=resolve(projectRoot),directory=resolve(root,'web');
 const sourceHtml=resolve(root,'app/comparison-page.html');
 const sourceEntry=resolve(root,'app/comparison-page.ts');
 const html=resolve(directory,'compare.html'),entry=resolve(directory,'compare.ts');
 const template=await readFile(sourceHtml,'utf8');
 await mkdir(directory,{recursive:true});
 await writeWhenChanged(html,template);
 await writeWhenChanged(entry,"import '../app/comparison-page';\n");
 return {html,entry,sourceHtml,sourceEntry};
}
