"""Actual native short-course root identity and full three-course plan only."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology,winding_contains
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.branch_junction_envelopes import close_branch
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops,signed_volume
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';OUT=P/'secondary-port-feasibility-v1';FRAME='model:local-ta98:assembled-metres'

def patch(v,f,faceids,rings,lo,hi):
 mask=np.isin(np.arange(len(f)),faceids)&np.all(np.isin(f,rings[lo:hi+1]),axis=1);ids=np.flatnonzero(mask);used=np.unique(f[ids]);lookup=np.full(len(v),-1,int);lookup[used]=np.arange(len(used));outer=v[used];inner=outer.copy()
 for ring in rings[lo:hi+1]:
  q=v[ring];center=q.mean(0);R=np.linalg.norm(q-center,axis=1).mean();thickness=min(.00012*np.sqrt(R/.001),.35*R);inner[lookup[ring]]=center+(q-center)*((R-thickness)/R)
 return outer,inner,lookup[f[ids]],ids,used,lookup

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);z=np.load(P/'source-port-charts.npz');v=z['branch_vertices'];f=z['branch_faces'];groups=recover_station_rings(v,f);g=[]
 for row in groups:
  rings,phase=align_station_correspondence(v,f[row['face_ids']],row['rings']);g.append(dict(**row,aligned_rings=np.array(rings),phase=phase))
 assert sorted(len(r['aligned_rings'])for r in g)==[25,181];main=next(r for r in g if len(r['aligned_rings'])==181);twig=next(r for r in g if len(r['aligned_rings'])==25);mr=main['aligned_rings'];tr=twig['aligned_rings'];assert np.array_equal(mr,z['branch_root_station_rings']);mc=v[mr].mean(1);tc=v[tr].mean(1)
 arr,summary=pair_contacts(Mesh(v,f[twig['face_ids']],FRAME),Mesh(v,f[main['face_ids']],FRAME),0.);assert summary['complete_within_pair_budget'];pairs=arr['pairs'][np.isin(arr['relation'],['crossing','coplanar-overlap'])];assert len(pairs)==58;actual=np.column_stack([twig['face_ids'][pairs[:,0]],main['face_ids'][pairs[:,1]]]);ts=twig['station_of_vertex'][f[actual[:,0]]];ms=main['station_of_vertex'][f[actual[:,1]]];trange=[int(ts.min()),int(ts.max())];mrange=[int(ms.min()),int(ms.max())];assert trange[1]==24
 distances=np.linalg.norm(mc-tc[-1],axis=1);station=int(distances.argmin());root=tc[-1];assert float(distances[station])<2e-7
 m0=max(0,mrange[0]-2);m1=min(len(mr)-1,mrange[1]+2);t0=max(0,trange[0]-2);t1=24;assert m1<176,'Local root domains may not overlap';mo,mi,mf,mids,mused,mlookup=patch(v,f,main['face_ids'],mr,m0,m1);to,ti,tf,tids,tused,tlookup=patch(v,f,twig['face_ids'],tr,t0,t1)
 # No finite parent classifier or source identity is guessed. Explicit root
 # center matching and actual inner fan containment are both prerequisites.
 _,parentfluid=close_branch(mo,mi,mf,{},origin=root)
 if signed_volume(parentfluid['positions'],parentfluid['faces'])<0:parentfluid['faces']=parentfluid['faces'][:,[0,2,1]]
 rootloop=set(tlookup[tr[-1]].tolist());loops=boundary_loops(tf);which=[i for i,loop in enumerate(loops)if set(loop.tolist())==rootloop];assert len(which)==1
 _,twigfluid=close_branch(to,ti,tf,{which[0]:root},origin=root);fan=twigfluid['positions'][np.unique(twigfluid['faces'][twigfluid['junction_faces']])];contains=winding_contains(Mesh(parentfluid['positions'],parentfluid['faces'],FRAME),np.r_[root[None],fan],boundary_tolerance_m=1e-10);clear=all(contains['inside'])and not any(contains['boundary_ambiguous'])
 arrays=dict(branch_vertices=v,branch_faces=f,main_rings=mr,twig_rings=tr,main_secondary_patch_face_ids=mids,twig_root_patch_face_ids=tids,main_source_face_ids=main['face_ids'],twig_source_face_ids=twig['face_ids'],actual_secondary_crossing_pairs=actual,main_root_to_M3_patch_face_ids=z['branch_local_patch_face_ids'],M3_vertices=z['parent_vertices'],M3_faces=z['parent_faces'],M3_rings=z['parent_station_rings'],M3_root_patch_face_ids=z['parent_local_patch_face_ids'],frame=np.array(FRAME));atomic_npz(OUT/'three-course-source-charts.npz',**arrays)
 proof=dict(schema='source-bound-postcentral-secondary-root-feasibility/v1',installed=False,geometry_modified=False,frame=FRAME,source_chart_sha256=digest_file(P/'source-port-charts.npz'),actual_contact_pairs=len(actual),main_contact_station_range=mrange,twig_contact_station_range=trange,main_station_at_twig_root=station,root_center_m=root.tolist(),root_center_match_residual_m=float(distances[station]),inner_root_fan_strict_containment=contains,finite_construction_feasible=bool(clear),explicit_secondary_local_domains=dict(main_station_interval=[m0,m1],main_source_faces=len(mids),twig_station_interval=[t0,t1],twig_source_faces=len(tids)),existing_M3_domains=dict(M3=[33,41],main=[176,180]),full_assembly_plan=dict(courses=['Complete native M3 (157 rings)','Complete native main postcentral source course (181 rings)','Complete native secondary postcentral source course (25 rings)'],junctions=['M3 station 36 to main postcentral root station 180','Main postcentral station '+str(station)+' to secondary root station 24'],inner_profiles='All source ring phases reused with authorized t(R)=min(0.12 mm*sqrt(R/1 mm),0.35 R), an educational assumption, not measured native wall.',topology='Outer union and blood union of three closed full-course classifiers; outer-minus-blood gives one tissue cell with four native annular terminal mouths. No artificial internal patch end rims.',protected='Every original outer face outside the two declared junction domains; original indices/packed normal/UV/material lineage and all brain/skeleton/context remain fixed.',stages=['Separate deterministic full classifiers and unions/difference, <=40 s, checkpoints before native operations.','Separate full intrinsic, four-mouth patency, full-source outer replay and source attribute checks, <=40 s.','Separate bounded recipient contact batches and actual source-normal surface/longitudinal section renders before installation.']),elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'feasibility.json',proof);print(json.dumps({k:proof[k]for k in ['actual_contact_pairs','main_contact_station_range','twig_contact_station_range','main_station_at_twig_root','root_center_match_residual_m','finite_construction_feasible','explicit_secondary_local_domains','elapsed_seconds']}))
 assert clear,'Derived root fluid fan is not wholly inside literal main-course fluid; stop before construction'

if __name__=='__main__':main()
