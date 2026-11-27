"""Source-only distal-overlap diagnostic; all overlays are measurement glyphs."""
from pathlib import Path
from collections import defaultdict
import sys,time,json,heapq
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1';P=BASE/'anterior-donor-counterpart-curves-v10/proof.json';OUT=BASE/'bp3d-distal-curve-diagnostic-v11';CAT=ROOT/'public/models/atlas.json';FRAME='source:BodyParts3D-4.0:packed-native-metres-y-up'
COLORS={56:[1.,.72,.03,1.],8:[.02,.9,.95,1.],6:[.95,.06,.83,1.]}

def glyph_segments(points,radius,closed=True):
    vv=[];ff=[];pairs=zip(points,np.roll(points,-1,axis=0))if closed else zip(points[:-1],points[1:])
    for a,b in pairs:
        axis=b-a;n=np.linalg.norm(axis)
        if n==0:continue
        axis/=n;u=np.cross(axis,np.eye(3)[np.argmin(abs(axis))]);u/=np.linalg.norm(u);v=np.cross(axis,u);r=min(radius,n*.12);theta=np.arange(6)*2*np.pi/6;ring=r*(np.cos(theta)[:,None]*u+np.sin(theta)[:,None]*v);off=len(vv);vv.extend(a+ring);vv.extend(b+ring)
        for j in range(6):
            k=(j+1)%6;ff.extend([[off+j,off+k,off+6+k],[off+j,off+6+k,off+6+j]])
    return np.array(vv),np.array(ff,int)

def locator(center,radius):
    v=center+radius*np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
    f=np.array([[0,2,4],[2,1,4],[1,3,4],[3,0,4],[2,0,5],[1,2,5],[3,1,5],[0,3,5]])
    return v,f

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);d=json.loads(P.read_text());a=json.loads(CAT.read_text());assert digest_file(CAT)==d['sourceCatalogueSha256'];parts={p['id']:p for p in a['parts']};reader=AtlasReader();entries=[];metrics=[];ancestry=[]
    def add(mesh,name,color):
        path=OUT/('view-'+str(len(entries))+'.npz');mesh.save(path);entries.append(dict(mesh=str(path),name=name,color=color));return str(path)
    for row in d['rows']:
        side=row['side'];aid=row['donorAca'];cid=row['donorDistalContinuation'];aca=reader.read_part(a,parts[aid],frame=FRAME,catalogue_path=CAT);calloso=reader.read_part(a,parts[cid],frame=FRAME,catalogue_path=CAT);components=row['sourceDistalContinuationIntersection']['components'];root=row['sourceRootIntersection']['components'][0];rootpoints=np.array(root['sourceCurvePoints_m']);rootcenter=np.array(root['arcLengthCentroid_m']);end=np.array(row['geodesicFarEndpointSourcePoints_m']).mean(0);direction=end-rootcenter;direction/=np.linalg.norm(direction)
        graph=defaultdict(dict)
        for tri in aca.faces:
            for x,y in zip(tri,np.roll(tri,-1)):
                length=float(np.linalg.norm(aca.vertices[x]-aca.vertices[y]));graph[int(x)][int(y)]=length;graph[int(y)][int(x)]=length
        rootfaces=np.unique(np.array(root['actualSourceFacePairs'])[:,0]);seed=np.unique(aca.faces[rootfaces]);dist=np.full(len(aca.vertices),np.inf);dist[seed]=0.;queue=[(0.,int(x))for x in seed]
        while queue:
            z,x=heapq.heappop(queue)
            if z!=dist[x]:continue
            for y,w in graph[x].items():
                if z+w<dist[y]:dist[y]=z+w;heapq.heappush(queue,(dist[y],y))
        supports=[np.array(c['actualSourceFacePairs'],int)for c in components]
        # Display colors partition literal source faces. Shared supports are
        # neutral; no face is duplicated or assigned an anatomical ownership.
        for owner,m in [(0,aca),(1,calloso)]:
            count=np.zeros(len(m.faces),int)
            for q in supports:count[np.unique(q[:,owner])]+=1
            groups=[('unaltered source faces',np.flatnonzero(count==0),[.76,.27,.20,1.]if owner==0 else[.58,.63,.68,1.]),('multiple-curve source face support',np.flatnonzero(count>1),[.9,.9,.9,1.])]
            for c,q in zip(components,supports):groups.append((str(c['nodes'])+'-node exact-curve source face support',np.unique(q[:,owner])[count[np.unique(q[:,owner])]==1],COLORS[c['nodes']]))
            allids=[]
            for label,ids,color in groups:
                if not len(ids):continue
                allids.extend(ids.tolist());source=dict(m.source,diagnostic_literal_face_partition=True,sourceFaceIds=ids.tolist(),anatomicalPartition=False);mesh=Mesh(m.vertices,m.faces[ids],FRAME,m.attributes,source)
                path=add(mesh,side+' '+parts[aid if owner==0 else cid]['name']+'; '+label,color);ancestry.append(dict(path=path,partId=aid if owner==0 else cid,sourceGeometrySha256=m.geometry_sha,sourceBankSha256=m.source['bank_sha256'],sourceFaceIds=ids.tolist(),sourceTrianglesAndNormalWordsLiteral=True,diagnosticColorsOnly=True))
            assert len(allids)==len(m.faces)and np.array_equal(np.sort(allids),np.arange(len(m.faces)))
        # White exact root curve establishes the upstream direction.
        v,f=glyph_segments(rootpoints,.000020);add(Mesh(v,f,FRAME,{},dict(measurementGlyphOnly=True,actualExactSourceIntersectionCurve=True,noAnatomicalTissue=True)),side+' upstream ACA–ICA exact root curve; WHITE measurement glyph',[.98,.98,.98,1.])
        for c,pairs in zip(components,supports):
            points=np.array(c['sourceCurvePoints_m']);center=np.array(c['arcLengthCentroid_m']);_,_,vt=np.linalg.svd(points-center,full_matrices=False);basis=vt[:2].T;xy=(points-center)@basis;cross=xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1);projected=.5*abs(cross.sum());normal=vt[-1];normal*=1 if normal[2]>=0 else-1;radius=.000012 if c['nodes']!=6 else.000001
            v,f=glyph_segments(points,radius);source=dict(measurementGlyphOnly=True,noAnatomicalTissue=True,centerlineIsLiteralExactSourceCurve=True,exactCurveSourceProof=str(P),sourceFacePairs=c['actualSourceFacePairs'],glyphMaximumRadius_m=radius)
            add(Mesh(v,f,FRAME,{},source),side+' '+str(c['nodes'])+'-node EXACT INTERSECTION CURVE; measurement line',COLORS[c['nodes']])
            # A labeled-by-color locator makes a 25–64 micrometre loop visible
            # in a whole-course image. It is not the curve's physical size.
            location=points[0]+normal*.0015;v,f=glyph_segments(np.array([points[0],location]),.000010,False);add(Mesh(v,f,FRAME,{},dict(measurementGlyphOnly=True,noAnatomicalTissue=True,notAnatomicalScale=True)),side+' '+str(c['nodes'])+'-node curve LOCATOR leader; not anatomy',COLORS[c['nodes']]);v,f=locator(location,.000080);add(Mesh(v,f,FRAME,{},dict(measurementGlyphOnly=True,noAnatomicalTissue=True,notAnatomicalScale=True)),side+' '+str(c['nodes'])+'-node curve LOCATOR marker; not curve size',COLORS[c['nodes']])
            sourceverts=np.unique(aca.faces[np.unique(pairs[:,0])]);metrics.append(dict(side=side,nodes=c['nodes'],color=COLORS[c['nodes']],perimeter_m=c['perimeter_m'],projectedBestFitPlaneArea_m2=float(projected),projectionIsNotAnAnatomicalAnnulus=True,axisAlignedExtent_m=np.ptp(points,axis=0).tolist(),principalExtents_m=np.ptp((points-center)@vt.T,axis=0).tolist(),planeMaximumResidual_m=c['planeMaximumResidual_m'],planeRmsResidual_m=c['planeRmsResidual_m'],curveArcLengthCentroid_m=c['arcLengthCentroid_m'],rootArcLengthCentroid_m=rootcenter.tolist(),rootToFarEndUnitDirection=direction.tolist(),projectedRootRelativePosition_m=float((center-rootcenter)@direction),sourceGraphDistanceRangeOnIncidentAcaControls_m=[float(dist[sourceverts].min()),float(dist[sourceverts].max())],graphDistanceIsNotAnatomicalCenterlineLength=True,actualSourceFacePairs=c['actualSourceFacePairs'],exactRationalKeys=c['exactRationalKeys'],allSourceGeometryFixed=True,notSelectedAsAnatomicalJunction=True,locatorOffset_m=.0015,locatorRadius_m=.000080))
    spec=dict(frame=FRAME,meshes=entries,views=[dict(name='BP3D complete upstream ACA and named callosomarginal course; anterior',direction=[0,0,1]),dict(name='Distal overlaps and exact 56 GOLD / 8 CYAN / 6 MAGENTA node curves; inferior',direction=[0,-1,.1],up=[0,0,1]),dict(name='Distal overlaps; left-oblique',direction=[1,-.25,1]),dict(name='Distal overlaps; right-oblique',direction=[-1,-.25,1])],resolution=1600,wall_seconds=60)
    atomic_json(OUT/'review-spec.json',spec);request=prepare(None,spec)
    proof=dict(schema='bp3d-distal-three-curve-source-diagnostic/v11',sourceCatalogueSha256=digest_file(CAT),counterpartProofSha256=digest_file(P),sourceFrame=FRAME,allUpstreamAcaAndDownstreamCallosomarginalSourceFacesShownExactlyOnce=True,allSourceGeometryAndNormalsLiteral=True,physicalSourceGeometryChanged=False,allAdditionalLinesAndLocatorsAreDiagnosticGlyphsOnly=True,noLargestLoopSelectedAsAnatomy=True,noFitDeformationWallReconstructionOrInstall=True,colorsByNodeCount=COLORS,metrics=metrics,sourceFaceAncestry=ancestry,queue=request,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',proof)
    print(json.dumps(dict(output=str(OUT/'proof.json'),reviewId=request['request_id'],meshes=len(entries),metrics=[dict(side=r['side'],nodes=r['nodes'],perimeter_mm=1000*r['perimeter_m'],projectedArea_mm2=1e6*r['projectedBestFitPlaneArea_m2'],principalExtents_mm=(1000*np.array(r['principalExtents_m'])).tolist())for r in metrics],elapsedSeconds=proof['elapsedSeconds'])))
if __name__=='__main__':main()
