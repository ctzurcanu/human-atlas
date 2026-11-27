"""Read-only native renal arterial stations, profile incidence and tissue witnesses."""
import json,sys,time
from pathlib import Path
import numpy as np
import point_cloud_utils as pcu
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/renal-blood-frontier-v1'
OUT=BASE/'arterial-course-v2'
PILOT=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-renal-source-pilot-v1'

def area(p):
    q=p-p.mean(0);return np.linalg.norm(np.cross(q,np.roll(q,-1,axis=-2)).sum(-2),axis=-1)/2

def profile(a,b):
    choices=[]
    for rev in [False,True]:
        for k in range(len(b)):
            order=np.roll(np.arange(len(b))[::-1]if rev else np.arange(len(b)),k);d=np.linalg.norm(a-b[order],axis=1)
            choices.append((float(d.max()),float(np.sqrt(np.mean(d*d))),rev,k,order))
    gap,rms,rev,k,order=min(choices,key=lambda q:q[:4]);return dict(max_control_gap_m=gap,rms_control_gap_m=rms,literal_coordinates=bool(np.array_equal(a,b[order])),reverse=rev,phase=k,
      relative_area_difference=float(abs(area(a)-area(b))/max(area(a),area(b))),control_order=order.tolist())

def stats(d):return dict(minimum_m=float(d.min()),maximum_m=float(d.max()),negative_controls=int(np.count_nonzero(d<0)),positive_controls=int(np.count_nonzero(d>0)))

def main():
    began=time.monotonic();assert not OUT.exists();ev=json.loads((BASE/'evidence.json').read_text());assert digest_file(ev['catalogue'])==ev['catalogue_sha256'];OUT.mkdir()
    bindings={r['part_id']:r for r in ev['source_bindings']}
    def load(pid):
        b=bindings[pid];assert digest_file(b['mesh'])==b['mesh_sha256'];return Mesh.load(b['mesh'])
    src=load('ZA:Intrarenal arteries of left kidney');z=np.load(PILOT/'rings-04.npz');v=z['physical_vertices'];f=z['physical_faces']
    assert np.array_equal(src.vertices[src.faces],v[f]),'All current arterial source triangles must replay exactly'
    courses=[];centers=[];controls=[];owners=[];station_ids=[];native=[];mouths=[];arrays={}
    for i in range(8):
        ids=z[f'component_{i}_rings'];r=v[ids];c=r.mean(1);aa=area(r);t=np.gradient(c,axis=0);t/=np.linalg.norm(t,axis=1)[:,None]
        faces=z[f'component_{i}_source_faces'];out=OUT/f'course-{i}.npz'
        atomic_npz(out,physical_vertices=v,source_faces=f[faces],actual_source_face_ids=faces,station_control_ids=ids,profiles=r,centers=c,tangents=t,areas_m2=aa)
        centers.append(c);controls.append(r);owners.extend([i]*len(c));station_ids.extend(range(len(c)));native.append((r,c,t))
        for end in [0,-1]:mouths.append(dict(component=i,station=0 if end==0 else len(c)-1,profile=r[end],outward=-t[0]if end==0 else t[-1]))
        courses.append(dict(component=i,stations=len(c),archive=str(out),sha256=digest_file(out),faces=len(faces),source_control_coordinates_literal=True,
          centerline_length_m=float(np.linalg.norm(np.diff(c,axis=0),axis=1).sum()),equivalent_diameter_m=np.sqrt(4*aa/np.pi)[[0,-1]].tolist(),
          ends_m=c[[0,-1]].tolist(),anatomical_flow_direction_verified=False,term_assignment=None))
    cc=np.concatenate(centers);rr=np.concatenate(controls);points=np.r_[cc,rr.reshape(-1,3)];arrays.update(centers=cc,profiles=rr,course_ids=np.array(owners),station_ids=np.array(station_ids))
    pairrows=[]
    for i,a in enumerate(mouths):
        for j,b in enumerate(mouths):
            if j<=i or a['component']==b['component']:continue
            m=profile(a['profile'],b['profile']);pairrows.append(dict(first_mouth=i,second_mouth=j,first_component=a['component'],second_component=b['component'],first_station=a['station'],second_station=b['station'],
                center_distance_m=float(np.linalg.norm(a['profile'].mean(0)-b['profile'].mean(0))),outward_tangent_dot=float(a['outward']@b['outward']),**m,
                source_graph_edge=bool(m['literal_coordinates']),near_profile_does_not_certify_attachment=True))
    atomic_json(OUT/'mouth-profile-incidence.json',dict(mouths=[{k:(value.tolist()if isinstance(value,np.ndarray)else value)for k,value in r.items()}for r in mouths],all_profile_pairs=pairrows,
      complete_source_components=8,literal_source_component_edges=sum(r['source_graph_edge']for r in pairrows),branch_junctions_in_source=0))
    tissues=['RENALDETAIL:male:left:'+n for n in ['capsule','cortex']+['pyramid:'+c for c in 'abcdefghijk']+['papilla:'+c for c in 'abcdefghijk']+['collecting-wall:'+c for c in ['adventitia','muscularis','mucosa']]]
    context=[]
    for n,pid in enumerate(tissues):
        assert time.monotonic()-began<50,'Tissue witness stage nearing deadline'
        m=load(pid);d,fi,bc=pcu.signed_distance_to_mesh(np.ascontiguousarray(points),np.ascontiguousarray(m.vertices),np.ascontiguousarray(m.faces,dtype=np.int32));d=np.asarray(d)
        control_d=d[len(cc):].reshape(-1,12);arrays['signed_distance_'+str(n)]=d;arrays['closest_face_'+str(n)]=fi;arrays['closest_barycentric_'+str(n)]=bc
        per=[]
        for i in range(8):
            mask=arrays['course_ids']==i;q=control_d[mask]
            per.append(dict(component=i,centers=stats(d[:len(cc)][mask]),profile_controls=stats(q),stations_with_any_negative_control=np.flatnonzero(np.any(q<0,axis=1)).tolist(),stations_with_all_negative_controls=np.flatnonzero(np.all(q<0,axis=1)).tolist()))
        context.append(dict(index=n,part_id=pid,source_geometry_sha256=m.geometry_sha,distances_from_actual_finite_tissue=True,classification='sampled profile containment, not complete swept-tube collision clearance',per_course=per))
        print(json.dumps(dict(stage='tissue',index=n,part=pid,seconds=round(time.monotonic()-began,3))),flush=True)
    # Old complete fluid classifier is a read-only exclusion reference; current
    # mucosal free sheet must be replayed against it before using volume signs.
    rp=ROOT/'reports/renal-ureter-fluid-v317.json';fr=json.loads(rp.read_text());npz=ROOT/fr['candidate'];assert digest_file(npz)==fr['candidateSha256'];q=np.load(npz)
    fluid_keys=q.files;fluid_meta=dict(report=str(rp),report_sha256=digest_file(rp),archive=str(npz),archive_sha256=digest_file(npz),arrays={k:list(q[k].shape)for k in fluid_keys},current_boundary_replay_not_yet_done=True)
    atomic_json(OUT/'collecting-fluid-source.json',fluid_meta)
    witnesses=OUT/'station-tissue-witnesses.npz';atomic_npz(witnesses,**arrays)
    refs=[dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC12408821/',pmid='40903539',use='Qualitative human segmental→interlobar column→arcuate pyramid surface→cortical radiate sequence; not a patient-specific geometry template.'),
          dict(url='https://pubmed.ncbi.nlm.nih.gov/2313791/',pmid='2313791',use='Human endocasts: upper calyceal arterial encirclement and variable anterior/posterior pelvic relations; do not infer a single variant from a name.')]
    report=dict(schema='renal-native-arterial-course-witnesses/v2',source_catalogue_sha256=ev['catalogue_sha256'],source_evidence_sha256=digest_file(BASE/'evidence.json'),source_station_archive=str(PILOT/'rings-04.npz'),source_station_sha256=digest_file(PILOT/'rings-04.npz'),
      courses=courses,contexts=context,witnesses=str(witnesses),witnesses_sha256=digest_file(witnesses),closest_profile_pairs=sorted(pairrows,key=lambda x:x['max_control_gap_m'])[:12],references=refs,
      no_geometry_or_labels_modified=True,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',report);print(json.dumps(dict(output=str(OUT/'evidence.json'),seconds=report['elapsed_seconds'],closest_pairs=report['closest_profile_pairs'][:4])))
if __name__=='__main__':main()
