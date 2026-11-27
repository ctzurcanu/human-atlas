"""Assess an explicit upward basal cage; save only a source-bound input plan."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'


def down(values,full,stop):
    t=np.clip((values-full)/(stop-full),0.,1.)
    return 1-t*t*(3-2*t)


def column_heights(points,triangles):
    # Literal triangle vertical intersections; no closing planes are invented.
    x,z=triangles[:,:,0],triangles[:,:,2];y=triangles[:,:,1]
    d=(z[:,1]-z[:,2])*(x[:,0]-x[:,2])+(x[:,2]-x[:,1])*(z[:,0]-z[:,2])
    keep=abs(d)>1e-20;x,z,y,d=x[keep],z[keep],y[keep],d[keep]
    low=np.full(len(points),np.nan);high=low.copy()
    for start in range(0,len(points),256):
        p=points[start:start+256]
        a=((z[:,1]-z[:,2])*(p[:,0,None]-x[:,2])+(x[:,2]-x[:,1])*(p[:,1,None]-z[:,2]))/d
        b=((z[:,2]-z[:,0])*(p[:,0,None]-x[:,2])+(x[:,0]-x[:,2])*(p[:,1,None]-z[:,2]))/d
        c=1-a-b;inside=(a>=-1e-11)&(b>=-1e-11)&(c>=-1e-11)
        h=a*y[:,0]+b*y[:,1]+c*y[:,2]
        minimum=np.where(inside,h,np.inf).min(1);maximum=np.where(inside,h,-np.inf).max(1)
        low[start:start+len(p)]=np.where(np.isfinite(minimum),minimum,np.nan)
        high[start:start+len(p)]=np.where(np.isfinite(maximum),maximum,np.nan)
    return low,high


def clipped_max_height(triangle,point,radius):
    polygon=[v.copy() for v in triangle]
    for axis,edge,sign in [(0,point[0]-radius,1),(0,point[0]+radius,-1),(2,point[1]-radius,1),(2,point[1]+radius,-1)]:
        if not polygon:return None
        result=[]
        for a,b in zip(polygon,polygon[1:]+polygon[:1]):
            da=(a[axis]-edge)*sign;db=(b[axis]-edge)*sign
            if da>=0:result.append(a)
            if (da>=0)!=(db>=0):result.append(a+(b-a)*(da/(da-db)))
        polygon=result
    return max(float(v[1]) for v in polygon) if polygon else None


def main(output):
    began=time.monotonic();deadline=began+45;output=output.resolve()
    assert output.is_relative_to(WORK) and not output.exists()
    source_path=WORK/'source-v1/source-audit.json';source=json.loads(source_path.read_text())
    catalogue_path=Path(source['sourceCatalogues']['male']['path']);catalogue_sha=digest_file(catalogue_path)
    assert catalogue_sha==source['sourceCatalogues']['male']['sha256']
    atlas=json.loads(catalogue_path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres'
    tongue=reader.read_part(atlas,parts['ZA:Tongue'],frame=frame,catalogue_path=catalogue_path,catalogue_sha=catalogue_sha)
    assert tongue.geometry_sha==source['contextParts'][0]['geometrySha256']
    bone_rows=[];bone_triangles=[]
    for row in source['contextParts']:
        if not row['partId'].startswith('TA98HYOID:'):continue
        assert parts[row['partId']]==row['sourcePart'] and digest_file(row['bankPath'])==row['bankSha256']
        mesh=reader.read_part(atlas,parts[row['partId']],frame=frame,catalogue_path=catalogue_path,catalogue_sha=catalogue_sha)
        assert mesh.geometry_sha==row['geometrySha256'];bone_triangles.append(mesh.triangles)
        bone_rows.append(dict(part_id=row['partId'],bank_path=row['bankPath'],bank_sha256=row['bankSha256'],geometry_sha256=mesh.geometry_sha,topology=topology(mesh)))
    bones=np.concatenate(bone_triangles);radius=.00025;step=.0005
    lower=np.maximum(tongue.vertices[:,[0,2]].min(0),bones[:,:,[0,2]].reshape(-1,2).min(0)-radius)
    upper=np.minimum(tongue.vertices[:,[0,2]].max(0),bones[:,:,[0,2]].reshape(-1,2).max(0)+radius)
    xx,zz=np.meshgrid(np.arange(lower[0],upper[0]+step/2,step),np.arange(lower[1],upper[1]+step/2,step));points=np.c_[xx.ravel(),zz.ravel()]
    basal,dorsal=column_heights(points,tongue.triangles)
    bl=bones[:,:,[0,2]].min(1);bh=bones[:,:,[0,2]].max(1);target=np.full(len(points),np.nan)
    for i,p in enumerate(points):
        if time.monotonic()>=deadline:raise TimeoutError('Basal cage assessment reached its 45 s internal deadline')
        if not np.isfinite(basal[i]):continue
        candidates=np.flatnonzero(np.all(bh>=p-radius,axis=1)&np.all(bl<=p+radius,axis=1))
        values=[v for k in candidates if (v:=clipped_max_height(bones[k],p,radius)) is not None]
        if values:target[i]=max(values)+radius
    valid=np.isfinite(target)&np.isfinite(basal)&np.isfinite(dorsal)
    required=np.where(valid,np.maximum(0.,target-basal),np.nan)
    lift=.010*down(points[:,1],.040,.055)*down(basal,1.610,1.633)
    upper_lift=.010*down(points[:,1],.040,.055)*down(dorsal,1.610,1.633)
    margin=basal+lift-target;crossing=valid&(required>0)
    vertex_dy=.010*down(tongue.vertices[:,2],.040,.055)*down(tongue.vertices[:,1],1.610,1.633)
    protected=(tongue.vertices[:,2]>=.055)|(tongue.vertices[:,1]>=1.633)
    assert np.all(vertex_dy[protected]==0)
    audit_path=WORK/'coverings-native-bone-audit-v1/summary.json';audit=json.loads(audit_path.read_text());face_ids=set()
    for row in audit['pairs']:
        if row['roi']!='native-whole' or not row['defect_pairs']:continue
        assert row['roi_geometry_sha256']==tongue.geometry_sha and digest_file(row['all_witnesses_npz'])==row['all_witnesses_sha256']
        with np.load(row['all_witnesses_npz'],allow_pickle=False) as z:face_ids.update(map(int,z['defect_pairs'][:,0]))
    faces=sorted(face_ids);vertices=np.unique(tongue.faces[faces]);extent=tongue.triangles[faces].reshape(-1,3)
    output.mkdir(parents=True)
    samples=output/'column-assessment.npz';atomic_npz(samples,xz=points,basal_y=basal,dorsal_y=dorsal,conservative_hyoid_clearance_top_y=target,required_upward_lift=required,proposed_analytic_lift=lift,proposed_upper_envelope_lift=upper_lift,proposed_analytic_margin=margin,valid=valid)
    plan=dict(schema='tongue-native-basal-upward-cage-input-plan/v1',source_catalogue_path=str(catalogue_path),source_catalogue_sha256=catalogue_sha,source_audit_path=str(source_path),source_audit_sha256=digest_file(source_path),source_part_id='ZA:Tongue',source_bank_path=tongue.source['bank_path'],source_bank_sha256=tongue.source['bank_sha256'],source_geometry_sha256=tongue.geometry_sha,frame=frame,
        fixed_hyoid=bone_rows,native_crossing_audit_path=str(audit_path),native_crossing_audit_sha256=digest_file(audit_path),native_crossing_pairs=280,affected_source_face_ids=faces,affected_source_vertex_ids=vertices.tolist(),affected_face_bounds_m=[extent.min(0).tolist(),extent.max(0).tolist()],
        field=dict(direction=[0,1,0],maximum_lift_m=.010,z_full_through_m=.040,z_zero_from_m=.055,y_full_through_m=1.610,y_zero_from_m=1.633,transition='C1 cubic smoothstep',x_and_z_literal=True,protected_source_vertex_ids=np.flatnonzero(protected).tolist(),source_vertices_with_nonzero_lift=int(np.count_nonzero(vertex_dy)),analytic_dy_derivative_lower_bound=1-.010*1.5/.023),
        assessment=dict(grid_step_m=step,sampled_columns=len(points),columns_with_native_tongue_and_hyoid=int(valid.sum()),columns_needing_upward_lift=int(crossing.sum()),maximum_sampled_required_lift_m=float(np.nanmax(required)),minimum_proposed_analytic_margin_m=float(np.nanmin(margin[valid])),sampled_columns_failing_proposed_analytic_clearance=int(np.count_nonzero(valid&(margin<0))),maximum_required_lift_over_native_thickness=float(np.max(required[crossing]/(dorsal[crossing]-basal[crossing]))),minimum_new_analytic_thickness_m=float(np.min(dorsal[valid]+upper_lift[valid]-basal[valid]-lift[valid])),maximum_posterior_upper_envelope_lift_m=float(np.max(upper_lift[valid])),upper_envelope_columns_moved=int(np.count_nonzero(valid&(upper_lift>0))),
            clearance_target_method='At each grid column, clip actual hyoid triangles to an X/Z square of radius0.25 mm, then add0.25 mm vertically. This is a conservative cube-clearance upper height at sampled columns; no bone caps or occupied-volume classifier.',sampled_plan_only=True),
        samples_path=str(samples),samples_sha256=digest_file(samples),limitations=['A343-vertex piecewise-linear emitted mesh is not the exact analytic cage; positive analytic Jacobian is not an emitted-triangle/self-clearance proof.','The plan does not certify mouth-floor, pharyngeal or native GG/HG insertion relationships after moving the base.','Actual emitted all-bone/clearance/self checks, source-normal transport, dorsal/front preservation and whole/section renders are required before regenerating layers or accepting anatomy.'],geometry_constructed=False,geometry_modified=False,skeleton_modified=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    assert digest_file(catalogue_path)==catalogue_sha
    atomic_json(output/'input-plan.json',plan)
    print(json.dumps(dict(output=str(output),affected_source_faces=len(faces),affected_source_vertices=len(vertices),field=plan['field'],assessment=plan['assessment'],geometry_constructed=False,elapsed_seconds=plan['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();main(args.output)
