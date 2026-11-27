from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def change(file,old,new):
 p=ROOT/file;s=p.read_text();assert old in s,(file,old[:90]);p.write_text(s.replace(old,new,1))
change('app/page.tsx','completion-v154/microscopic-atlas.json','completion-v155/microscopic-atlas.json')
change('app/local-ta98-layer-catalogue.ts',"completion-v155/atlas.json","completion-v155/atlas-ready.json")
change('app/scene.tsx',"const surfaceHitVisible=(index:number,hit:T.Intersection)=>{", "const surfaceHitVisible=(index:number,hit:T.Intersection)=>{\n   if((atlas.parts[index] as Part&{liverComputationalFaces?:number[]}).liverComputationalFaces?.includes(hit.faceIndex??-1))return false;")
change('app/scene.tsx','for(const group of p.tissueDrawGroups)displayed.addGroup(group.start,group.count,0);','for(const group of p.tissueDrawGroups)if(!(group as {computationalOnly?:boolean}).computationalOnly)displayed.addGroup(group.start,group.count,0);')
change('app/scene.tsx','for(const group of atlas.parts[index].tissueDrawGroups??[]){const replacement=', 'for(const group of atlas.parts[index].tissueDrawGroups??[]){if((group as {computationalOnly?:boolean}).computationalOnly)continue;const replacement=')
print('Functional section and mouth display integrated')
