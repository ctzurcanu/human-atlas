"""Bound source-only circle-of-Willis diagnostics; never fits anatomy."""
from pathlib import Path
import sys,time,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1'
OUT=BASE/'native-aca-source-review-v5'
FRAME='model:local-ta98:assembled-metres'

def section(m,center,axis,level):
    scalar=(m.vertices-center)@axis-level
    clipped=clip_scalar_surface(m.vertices,m.faces,scalar,0.)
    loops=boundary_loops(clipped['negative_faces']);assert len(loops)==1
    points=clipped['positions'][loops[0]];pivot=center+axis*level
    u=np.cross(axis,np.eye(3)[np.argmin(abs(axis))]);u/=np.linalg.norm(u);v=np.cross(axis,u);basis=np.c_[u,v]
    xy=(points-pivot)@basis;cross=xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1);signed=.5*cross.sum();assert signed!=0
    centroid=((xy+np.roll(xy,-1,axis=0))*cross[:,None]).sum(0)/(6*signed)
    diam=np.linalg.norm(points[:,None]-points[None,:],axis=2).max()
    return dict(level_m=float(level),sourceSectionPoints_m=points.tolist(),area_m2=float(abs(signed)),equivalentOuterDiameter_m=float(2*np.sqrt(abs(signed)/np.pi)),maximumFeretOuterDiameter_m=float(diam),centroid_m=(pivot+basis@centroid).tolist(),sourceCapsCreated=0)

def main():
    started=time.monotonic();assert digest_file(CAT)==EXPECTED;OUT.mkdir(exist_ok=False)
    atlas=json.loads(CAT.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();meshes={};sources=[];entries={};ancestry={}
    def read(pid):
        if pid not in meshes:meshes[pid]=reader.read_part(atlas,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED)
        return meshes[pid]
    def entry(pid,color,faces=None,suffix=''):
        m=read(pid);ids=np.arange(len(m.faces))if faces is None else np.asarray(faces,int);key=pid+suffix
        source=dict(m.source,diagnostic_source_face_ids=ids.tolist(),diagnostic_selection_only=True,no_emitted_cap=True,no_source_geometry_changed=True)
        q=Mesh(m.vertices,m.faces[ids],FRAME,m.attributes,source);path=OUT/('source-'+str(len(entries))+'.npz');q.save(path)
        assert np.array_equal(q.triangles,m.triangles[ids])
        entries[key]=dict(mesh=str(path),name=key,color=color)
        ancestry[key]=dict(sourcePartId=pid,sourceName=parts[pid]['name'],sourceTerm=parts[pid].get('ta98Term'),sourceLicence=parts[pid].get('license'),sourceBankSha256=m.source['bank_sha256'],sourceGeometrySha256=m.geometry_sha,sourceAttributeSha256=m.attributes_sha,sourceFaceIds=ids.tolist(),snapshot=str(path),snapshotSha256=digest_file(path),allIndexedSourceCoordinatesAndNormalWordsLiteral=True,anatomicalRegistrationVerified=False)
        return entries[key]
    arteries=['ZA:Anterior cerebral artery.l','ZA:Anterior cerebral artery.r','ZA:Anterior communicating artery','ZA:Internal carotid artery.l','ZA:Internal carotid artery.r']
    colors=[[.93,.24,.16,1],[.98,.56,.20,1],[.80,.10,.68,1],[.60,.10,.08,1],[.63,.32,.05,1]]
    vessel=[entry(p,c)for p,c in zip(arteries,colors)]
    chiasm=[entry('ZA:Optic chiasm.'+s,[.96,.88,.18,1])for s in ['l','r']]
    optic=[entry('ZA:Optic tract.'+s,[.75,.68,.16,1])for s in ['l','r']]
    sphenoid=entry('ZA:Sphenoid bone:Bone-4',[.65,.72,.79,1])
    brainids=[p['id']for p in atlas['parts']if p['id'].startswith('ZA:')and p['id'].rsplit(':',1)[-1]in ['Frontal lobe','Temporal lobe','Parietal lobe','Occipital lobe']]
    vv=[];ff=[];normals=[];groups=[];offset=0
    for pid in brainids:
        m=read(pid);vv.append(m.vertices);ff.append(m.faces+offset);normals.append(m.attributes['normals']);offset+=len(m.vertices)
        groups.append(dict(partId=pid,triangles=len(m.faces),sourceBankSha256=m.source['bank_sha256'],geometrySha256=m.geometry_sha,attributeSha256=m.attributes_sha))
    brain=Mesh(np.vstack(vv),np.vstack(ff),FRAME,dict(normals=np.vstack(normals)),dict(diagnostic_literal_context_aggregate=True,sourceCatalogueSha256=EXPECTED,sourceGroups=groups,no_source_geometry_changed=True))
    brainpath=OUT/'brain-context-literal.npz';brain.save(brainpath);brainentry=dict(mesh=str(brainpath),name='Fixed ZA cortex context; literal outer-lobe source faces',color=[.55,.54,.52,1])
    atomic_json(OUT/'brain-context-source-groups.json',groups)
    # These BP3D records are posterior segments, not an anterior junction.
    bp=[entry(pid,[.10,.78,.92,1])for pid in ['FJ1723','FJ1723M']if pid in parts]
    roots=[]
    for pid,c in zip(arteries[-2:],colors[-2:]):
        m=read(pid);ids=np.flatnonzero(np.all(m.triangles[:,:,1]>=1.704,axis=1));assert len(ids)
        roots.append(entry(pid,c,ids,' [literal intracranial-root face selection; uncapped]'))
    # Derived geometric centerline, explicitly not a supplied anatomical course.
    comm=read(arteries[2]);center=comm.vertices.mean(0);_,_,vt=np.linalg.svd(comm.vertices-center,full_matrices=False);axis=vt[0];axis*=1 if axis[0]>0 else-1
    projection=(comm.vertices-center)@axis;lo=float(projection.min());hi=float(projection.max())
    sections=[section(comm,center,axis,x)for x in np.linspace(lo,hi,21)[1:-1]]
    end0=comm.vertices[np.flatnonzero(projection==lo)].mean(0);end1=comm.vertices[np.flatnonzero(projection==hi)].mean(0)
    derived=np.vstack([end0,*[r['centroid_m']for r in sections],end1]);length=float(np.linalg.norm(np.diff(derived,axis=0),axis=1).sum())
    exact=json.loads((BASE/'native-aca-bounded-references-v3/exact-witnesses-v4/proof.json').read_text());gaps=[];planes=[]
    for r in exact['rows']:
        w=r['closestActualSurfaceWitness'];g=np.array(w['firstPoint_m'])-w['secondPoint_m'];outward=axis if r['side']=='l'else-axis;n=np.cross(axis,g);n/=np.linalg.norm(n)
        gaps.append(dict(side=r['side'],acomToAcaGapVector_m=g.tolist(),gap_m=float(np.linalg.norm(g)),angleToAComLongAxis_degrees=float(np.degrees(np.arccos(np.clip(g@outward/np.linalg.norm(g),-1,1)))),acaFace=w['firstSourceFace'],acomFace=w['secondSourceFace']))
        planes.append(dict(name=r['side']+'-closest-wall-longitudinal-cut',direction=n.tolist(),up=axis.tolist(),section=dict(normal=n.tolist(),point=w['secondPoint_m'],clear_outer=True)))
    measurements=dict(schema='native-acom-source-course/v5',sourceCatalogueSha256=EXPECTED,sourcePartId=arteries[2],sourceGeometrySha256=comm.geometry_sha,sourceBankSha256=comm.source['bank_sha256'],storedSourceCenterlineExists=False,derivedCenterlineMethod='19 literal closed-surface planar contour area centroids at 5–95 percent of PCA extent, joined to literal extreme source controls; geometric diagnostic, not anatomical ostia or an authoritative centerline',longAxis=axis.tolist(),axisExtent_m=hi-lo,axisEndpointSourceCoordinates_m=[end0.tolist(),end1.tolist()],derivedCenterline_m=derived.tolist(),derivedPolylineLength_m=length,sections=sections,sourceAxisAlignedDimensions_m=np.ptp(comm.vertices,axis=0).tolist(),gapWitnesses=gaps,outerCastNotFiniteWallOrLumen=True,noAnatomicalRepairJustifiedFromLengthAlone=True,sourceRegistrationUnverified=True)
    atomic_json(OUT/'source-course-measurements.json',measurements)
    whole=dict(frame=FRAME,meshes=vessel+chiasm+optic+[sphenoid,brainentry]+bp,views=[dict(name='native-source-inferior-context; BP3D posterior segments cyan',direction=[0,-1,.15],up=[0,0,1]),dict(name='native-source-anterior-oblique-context',direction=[.5,-.4,1]),dict(name='native-source-left-context',direction=[1,0,.2])],resolution=1100,wall_seconds=60)
    focus=dict(frame=FRAME,meshes=vessel[:3]+roots+chiasm+bp,views=[dict(name='complete-bilateral-ACA-courses-and-closed-ACom; inferior',direction=[0,-1,0],up=[0,0,1]),dict(name='complete-bilateral-ACA-courses-and-closed-ACom; anterior',direction=[0,0,1]),dict(name='complete-bilateral-ACA-courses-and-closed-ACom; superior-oblique',direction=[.2,1,.4],up=[0,0,1])],resolution=1100,wall_seconds=60)
    cuts=dict(frame=FRAME,meshes=vessel[:3]+roots+chiasm,views=planes+[dict(name='near-junction-axial-source-cut; uncapped',direction=[0,1,0],up=[0,0,1],section=dict(normal=[0,1,0],point=center.tolist(),clear_outer=True)),dict(name='near-junction-coronal-source-cut; uncapped',direction=[0,0,1],section=dict(normal=[0,0,1],point=center.tolist(),clear_outer=True))],resolution=1100,wall_seconds=60)
    queued=[]
    for name,spec in [('whole-context',whole),('focused-source-course',focus),('local-source-sections',cuts)]:
        atomic_json(OUT/(name+'-spec.json'),spec);result=prepare(None,spec);result['purpose']=name;queued.append(result)
    atomic_json(OUT/'source-face-ancestry.json',ancestry)
    proof=dict(schema='native-aca-source-review-queue/v5',sourceCatalogueSha256=EXPECTED,sourceCatalogue=str(CAT),allSourceCoursesAndContextFixed=True,noGeometryOrLiveMetadataChanged=True,noRegistrationAssumptionForBp3d=True,bp3dActualIdentity='FJ1723/FJ1723M: right/left precommunicating posterior cerebral artery A12.2.07.083; diagnostic alternative only, not a whole circle or ACom',brainContextParts=len(brainids),brainContextTriangles=len(brain.faces),queue=queued,measurements=str(OUT/'source-course-measurements.json'),elapsedSeconds=time.monotonic()-started)
    atomic_json(OUT/'proof.json',proof);assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(output=str(OUT/'proof.json'),queue=[dict(id=q['request_id'],purpose=q['purpose'])for q in queued],axisExtent_mm=1000*(hi-lo),derivedCourse_mm=1000*length,outerEquivalentDiameterRange_mm=[1000*min(r['equivalentOuterDiameter_m']for r in sections),1000*max(r['equivalentOuterDiameter_m']for r in sections)],gapAngles_degrees=[q['angleToAComLongAxis_degrees']for q in gaps],elapsedSeconds=proof['elapsedSeconds'])))

if __name__=='__main__':main()
