"""Independently check source gauges and actual section clearance of a fit.

Usage: python -m geometry_tools.commands.audit_tube_course FIT_DIRECTORY --output NEW_REPORT.json
Fixed-plane calibre is a diagnostic, not a perpendicular lumen measurement.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import shapely
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from geometry_tools.course_io import sha,load_arrays,load_course_spec,write_json


def area(q):return float(np.linalg.norm(np.cross(q-q.mean(0),np.roll(q-q.mean(0),-1,axis=0)).sum(0))*.5)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('fit_directory');parser.add_argument('--output',required=True);args=parser.parse_args();start=time.monotonic();folder=Path(args.fit_directory).resolve();out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new report path')
    fit=json.loads((folder/'run.json').read_text());spec,g,source,plan,course,active,fixed,_,_=load_course_spec(fit['specPath'])
    if fit['bindings']['specSha256']!=sha(fit['specPath'])or fit.get('displacementsSha256')!=sha(folder/fit['displacements']):raise ValueError('Stale fit artifact')
    d=load_arrays(folder/fit['displacements'])['delta']
    if d.shape!=g['centres'].shape or not np.isfinite(d).all():raise ValueError('Invalid displacement field')
    protected=bool(np.array_equal(d[list(fixed)],np.zeros((len(fixed),3))));y_literal=bool(np.array_equal(d[:,1],np.zeros(len(d))));v=g['vertices']+d[g['cells']];gauge_errors=[];max_pair=0.;max_area=0.
    for p in source['profiles']:
        raw=np.array(p['rawIds'])+p['offset'];a=g['vertices'][raw];b=v[raw];pair=float(np.abs(np.linalg.norm(a[:,None]-a[None,:],axis=2)-np.linalg.norm(b[:,None]-b[None,:],axis=2)).max());before=area(a);after=area(b)
        relative=abs(after/before-1)if before>0 else None;max_pair=max(max_pair,pair)
        if relative is not None:max_area=max(max_area,relative)
        if relative is None or pair>5e-8 or relative>.005:gauge_errors.append(dict(partId=p['partId'],cell=p['cell'],pairDistanceError_m=pair,relativeAreaError=relative))
    sections=[];clear=True;fixed_plane_area_max=0.;radius_change_max=0.
    for row in course['rows']:
        context={k:shapely.from_wkb(val)for k,val in row['contextWKB'].items()};polygons=[];profiles=[]
        for p in row['profiles']:
            edges=np.array(p['sourceRawEdges']);weights=np.array(p['sourceEdgeWeights']);before=(g['vertices'][edges]*weights[:,:,None]).sum(1);after=(v[edges]*weights[:,:,None]).sum(1);polygon=shapely.polygons(after[:,[0,2]]);valid=bool(shapely.is_valid(polygon));oldarea=area(before);relative=abs(area(after)/oldarea-1)if oldarea>0 else None
            if relative is not None:fixed_plane_area_max=max(fixed_plane_area_max,relative)
            before_r=np.linalg.norm(before-before.mean(0),axis=1);after_r=np.linalg.norm(after-after.mean(0),axis=1);radius_change_max=max(radius_change_max,float(np.abs(after_r-before_r).max()))
            distances={};fits=False
            if valid:
                distances={k:float(shapely.distance(polygon,shape))for k,shape in context.items()if k!='container'};distances['container']=float(shapely.distance(polygon,shapely.boundary(context['container'])))
                fits=bool(shapely.contains(context['container'],polygon))and all(shapely.disjoint(polygon,shape)for k,shape in context.items()if k!='container')and min(distances.values())>=plan['physicalClearance_m']-1e-8
            clear&=fits;polygons.append(polygon);profiles.append(dict(partId=p['partId'],valid=valid,clear=fits,clearances_m=distances,fixedPlaneAreaRelativeError=relative))
        mutual=[]
        for i,a in enumerate(polygons):
            for j,b in enumerate(polygons[i+1:],i+1):
                gap=float(shapely.distance(a,b));ok=bool(shapely.disjoint(a,b))and gap>=plan['physicalClearance_m']-1e-8;clear&=ok;mutual.append(dict(profiles=[i,j],clearance_m=gap,clear=ok))
        sections.append(dict(index=row['index'],profiles=profiles,mutual=mutual))
        if time.monotonic()-start>15:raise TimeoutError('15s native independent course audit deadline')
    source_tri=g['vertices'][g['faces']];candidate_tri=v[g['faces']];oldnormal=np.cross(source_tri[:,1]-source_tri[:,0],source_tri[:,2]-source_tri[:,0]);normal=np.cross(candidate_tri[:,1]-candidate_tri[:,0],candidate_tri[:,2]-candidate_tri[:,0]);introduced_zero=int(np.count_nonzero((np.linalg.norm(oldnormal,axis=1)>1e-20)&(np.linalg.norm(normal,axis=1)<=1e-20)))
    report=dict(schema='independent-tube-course-audit/v1',fitReportSha256=sha(folder/'run.json'),specSha256=sha(fit['specPath']),sourceCodeSha256=sha(__file__),qualifiedCourse=fit['qualifiedCourse'],protectedCellsLiteral=protected,yCoordinatesLiteral=y_literal,sourceProfilesChecked=len(source['profiles']),sourceGaugeFailures=len(gauge_errors),sourceGaugeFailureDetails=gauge_errors,maximumSourcePairDistanceError_m=max_pair,maximumSourceAreaRelativeError=max_area,allSampledSectionsClear=bool(clear),sectionEvidence=sections,maximumFixedPlaneAreaRelativeError=fixed_plane_area_max,maximumFixedPlaneRadiusChange_m=radius_change_max,fixedPlaneAreasAreNotNormalLumenMeasurements=True,introducedZeroAreaFaces=introduced_zero,readyForWholeSurfaceReview=bool(fit['qualifiedCourse']and protected and y_literal and not gauge_errors and clear and not introduced_zero),selfAndWholeBodySurfaceContactsNotChecked=True,normalsAndPackedEncodingNotYetChecked=True,installed=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start)
    write_json(out,report);print(json.dumps({k:val for k,val in report.items()if k not in ['sectionEvidence','sourceGaugeFailureDetails']}))
    if not report['readyForWholeSurfaceReview']:raise SystemExit(1)


if __name__=='__main__':main()
