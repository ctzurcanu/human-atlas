"""One explicit exposed-tongue mucosal carrier, separate from basal envelopes.

The posterior/inferior boundary represents reflection into floor-of-mouth
mucosa. It is not a cut around individual muscle envelopes or an opening solid.
"""
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'

def load_sources():
    whole=WORK/'basal-whole-protected-v1/whole-tongue.npz'
    items={'whole':whole,'GG-left':WORK/'genioglossus-v2-encoding-safe-final/genioglossus-left.npz',
           'GG-right':WORK/'genioglossus-v2-encoding-safe-final/genioglossus-right.npz',
           'HG-left':WORK/'source-v1/meshes/native-03.npz','HG-right':WORK/'source-v1/meshes/native-04.npz',
           'SG-left':WORK/'styloglossus-v2/styloglossus-left.npz','SG-right':WORK/'styloglossus-v2/styloglossus-right.npz',
           'mandible':WORK/'source-v1/meshes/native-05.npz','hyoid':WORK/'source-v1/meshes/native-13.npz'}
    return {key:(path,Mesh.load(path))for key,path in items.items()}

def segments(mesh,axis,level):
    output=[]
    for tri in mesh.triangles:
        d=tri[:,axis]-level
        if d.min()>0 or d.max()<0:continue
        hits=[]
        for a,b in [(0,1),(1,2),(2,0)]:
            if d[a]*d[b]<0:hits.append(tri[a]+d[a]/(d[a]-d[b])*(tri[b]-tri[a]))
        if len(hits)==2:output.append(hits)
    return np.asarray(output).reshape(-1,2,3)

def inspect(output):
    from PIL import Image,ImageDraw
    began=time.monotonic();output=output.resolve();output.mkdir(parents=True,exist_ok=False)
    inputs=load_sources();fig=Image.new('RGB',(1200,1000),'white');draw=ImageDraw.Draw(fig)
    colors={'whole':'black','GG-left':'#bf562e','GG-right':'#bf562e','HG-left':'#207aa2','HG-right':'#207aa2',
            'SG-left':'#7f42a3','SG-right':'#7f42a3','mandible':'#a8a596','hyoid':'#8e8c81'}
    views=[('Near-midline sagittal X=3 mm',0,.003,(2,1)),('Lateral sagittal X=11 mm',0,.011,(2,1)),
           ('Coronal Z=60 mm',2,.060,(0,1)),('Axial Y=1622 mm',1,1.622,(0,2))]
    for index,(title,axis,level,coords) in enumerate(views):
        origin=np.array([(index%2)*600+50,(index//2)*500+50]);extent=np.array([[20,1595],[90,1655]])if coords==(2,1)else(np.array([[-25,1595],[25,1655]])if coords==(0,1)else np.array([[-25,20],[25,90]]))
        scale=min(500/(extent[1,0]-extent[0,0]),400/(extent[1,1]-extent[0,1]))
        def xy(p):
            value=(p-extent[0])*scale
            return tuple(origin+[value[0],400-value[1]])
        draw.text(tuple(origin-[0,30]),title,fill='black');draw.rectangle(tuple(origin)+tuple(origin+[500,400]),outline='#dddddd')
        for key,(_,mesh)in inputs.items():
            lines=segments(mesh,axis,level)
            for line in lines[:,:,coords]*1000:draw.line([xy(p)for p in line],fill=colors[key],width=2)
    draw.text((50,970),'Black: corrected tongue | rust: GG | blue: HG | purple: SG | grey: fixed bone',fill='black');fig.save(output/'source-sections.png')
    atomic_json(output/'inspection.json',dict(schema='tongue-exposed-mucosa-inspection/v1',
        sources={key:dict(path=str(path),sha256=digest_file(path),geometrySha256=mesh.geometry_sha)for key,(path,mesh)in inputs.items()},
        figure=str(output/'source-sections.png'),elapsed_seconds=time.monotonic()-began,anatomical_acceptance=False))
    print(str(output/'source-sections.png'),flush=True)

def reflection_height(points):
    """Approximate floor-of-mouth reflection; no muscle-distance classifier."""
    p=np.asarray(points,float)
    def smooth(values):
        t=np.clip(values,0.,1.);return t*t*(3-2*t)
    return 1.6265-.0025*smooth((p[:,2]-.055)/.017)-.002*smooth(np.abs(p[:,0])/.016)

def construct(output):
    from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface,interpolate_clip_controls
    from geometry_tools.algorithms.tube_continuity import boundary_loops
    from geometry_tools.audit import topology
    from geometry_tools.tissue import layer_arrangement,normal_words,pack
    from geometry_tools.contact import pair_contacts
    began=time.monotonic();output=output.resolve()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new local tongue output directory')
    inputs=load_sources();source,whole=inputs['whole'];v,f=whole.vertices,whole.faces
    values=v[:,1]-reflection_height(v)
    clipped=clip_scalar_surface(v,f,values,0.);outer=clipped['positions'];faces=clipped['positive_faces'];parents=clipped['positive_parents']
    source_normals=whole.attributes['normals'].astype(float);source_normals/=np.linalg.norm(source_normals,axis=1)[:,None]
    thickness=.00035+.00007*np.clip(source_normals[:,1],0,1)-.00005*np.clip(-source_normals[:,1],0,1)
    inner=interpolate_clip_controls(v-thickness[:,None]*source_normals,clipped['new_edge_parents'],clipped['fractions'])
    new_normals=interpolate_clip_controls(source_normals,clipped['new_edge_parents'],clipped['fractions'])
    new_normals/=np.linalg.norm(new_normals,axis=1)[:,None]
    raw_words=np.rint(new_normals*32767).astype('<i2');raw_words[:len(v)]=whole.attributes['normals']
    carrier=Mesh(outer,faces,whole.frame,dict(normals=raw_words),dict(whole.source,derived_from_geometry=whole.geometry_sha))
    carrier_topology=topology(carrier);loops=boundary_loops(faces)
    if len(carrier_topology['components'])!=1 or len(loops)!=1:raise ValueError('The exposed mucosa must be one patch with one reflection boundary')
    protected_faces=np.all((v[f,1]>=1.633)|(v[f,2]>=.072),axis=1)
    full_literal=np.all(faces<len(v),axis=1)&np.all(faces==f[parents],axis=1)
    if not set(np.flatnonzero(protected_faces))<=set(parents[full_literal]):raise ValueError('Protected dorsal/free-blade source facets were cut or removed')
    plan=dict(schema='tongue-exposed-mucosa-boundary/v1',frame=whole.frame,source=str(source),sourceSha256=digest_file(source),
        sourceGeometrySha256=whole.geometry_sha,sourceCatalogueSha256=whole.source.get('catalogue_sha256'),
        references=[dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC6850428/',evidence='Adult fresh-cadaver anatomy: mucosa/floor-of-mouth fascia reflect onto anterior and lateral ventral tongue; GG is suspended below.'),
                    dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',evidence='Adult serial sections: body has a mushroom outline elevated by GG; free blade lies anterior to frenulum; extrinsic muscles insert within tongue.')],
        field=dict(kind='smooth-approximate-FOM-reflection-not-individual-muscle-holes',
                   formula='Y-1.6265+.0025*smoothstep((Z-.055)/.017)+.002*smoothstep(abs(X)/.016) >= 0',
                   referenceInformedModelingBoundary=True,sourceFrenulumNotEstablished=True,
                   kept='dorsum, lateral exposed body/root and genuine anterior ventral blade',
                   removed='computational posterior/inferior root boundary continuous with extrinsic-muscle entry'),
        protected=dict(dorsalYAtLeastM=1.633,freeBladeZAtLeastM=.072,literalWholeFaces=int(protected_faces.sum())),
        fixed=['All native bones and whole tongue source controls'],
        adjoiningInterfaces=[dict(first='mucosa',second='floor-of-mouth-mucosa',kind='reflection-boundary; neighboring FOM tissue is not yet constructed'),
                             dict(first='mucosa',second='extrinsic-muscle-insertions',kind='intralingual fit still required; crossings are not accepted')],
        anatomyAcceptance=False,installed=False)
    output.mkdir(parents=True);atomic_json(output/'reflection-boundary.json',plan)
    carrier.save(output/'carrier.npz')
    Mesh(outer,clipped['negative_faces'],whole.frame).save(output/'excluded-basal-carrier.npz')
    fields=np.array([outer,inner]);atomic_npz(output/'fields.npz',fields=fields,faces=faces,frame=np.array(whole.frame))
    rows=[dict(sourceFace=int(parent),carrierFace=int(i),cornerOrder=[0,1,2])for i,parent in enumerate(parents)if full_literal[i]]
    atomic_json(output/'source-face-provenance.json',dict(schema='tongue-exposed-carrier-ancestry/v1',
        sourceGeometrySha256=whole.geometry_sha,carrierGeometrySha256=carrier.geometry_sha,faces=rows,
        cutFaces=[dict(carrierFace=int(i),sourceFace=int(parent))for i,parent in enumerate(parents)if not full_literal[i]],
        edgeSubdivision=dict(sourceEdges=clipped['new_edge_parents'].tolist(),fractions=clipped['fractions'].tolist()),
        sourceControlsLiteralPrefix=len(v),originalNormalsLiteralPrefix=len(v)))
    atomic_json(output/'reflection-curve.json',dict(schema='anatomical-boundary-curve/v1',frame=whole.frame,
        source='floor-of-mouth reflection modeling boundary',loops=[outer[loop].tolist()for loop in loops],notSeparateSolid=True))
    positions,surfaces,cells,prisms=layer_arrangement(fields,faces,['mucosa'])
    words=normal_words(positions,surfaces);words['interface:0']=raw_words[faces]
    result=pack(positions,surfaces,cells,{'mucosa':'TA98TONGUE:male:mucosa'},whole.frame,output/'packed',
        origin=(0,1.625,0),corner_words=words,construction=prisms,max_pairs=100000,max_tests=1000000)
    emitted=Mesh.load(output/'packed/TA98TONGUE_male_mucosa.npz');pairs=[]
    for name,(path,other)in inputs.items():
        if name=='whole':continue
        if time.monotonic()-began>55:raise TimeoutError('Stopped before 60-second exposed-mucosa stage deadline')
        arrays,summary=pair_contacts(emitted,other,0.,max_pairs=100000,max_tests=1000000)
        witness=output/(name+'-contact.npz');atomic_npz(witness,**arrays)
        pairs.append(dict(other=name,otherPath=str(path),otherSha256=digest_file(path),summary=summary,witnesses=str(witness)))
    report=dict(**plan,carrier=str(output/'carrier.npz'),fields=str(output/'fields.npz'),packedManifest=str(output/'packed/tissue.json'),
        carrierTopology=carrier_topology,carrierTriangles=len(faces),sourceFacesRetainedLiteral=len(rows),cutTriangles=int((~full_literal).sum()),
        removedBasalTriangles=len(clipped['negative_faces']),finiteWallChecks=result['checks'],finiteWallPassed=result['passed'],
        actualEmittedPairs=pairs,elapsedSeconds=time.monotonic()-began,
        completeAssemblyPassed=False,remaining=['Fit intralingual GG/HG/SG to inner envelope above the reflection boundary',
            'Root must couple dorsal aponeurosis using the same shared inner interface',
            'Floor-of-mouth and palatoglossal mucosal continuity are pending; no open reflection edge is claimed complete'])
    atomic_json(output/'candidate.json',report)
    print(dict(finiteWallPassed=result['passed'],carrierTriangles=len(faces),remainingCrossings={r['other']:r['summary']['crossing_face_pairs']for r in pairs},elapsedSeconds=report['elapsedSeconds']),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inspect',action='store_true');p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.inspect:inspect(args.output)
    else:construct(args.output)
