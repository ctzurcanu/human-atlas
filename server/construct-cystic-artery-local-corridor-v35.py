"""One bounded source-profile-preserving local cystic arterial course correction."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
from scipy.optimize import lsq_linear
from scipy.spatial.distance import cdist
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.algorithms.tube_sections import bend_diagnostics
from geometry_tools.algorithms.tube_prism_constraints import StationPrismConstraints
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
sp=importlib.util.spec_from_file_location('v14',ROOT/'server/construct-admissible-cystic-departure-v14.py');v14=importlib.util.module_from_spec(sp);sp.loader.exec_module(v14);v10=v14.v10
WORK=v14.WORK;OUT=WORK/'artery-local-corridor-v35';CURRENT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345';FRAME=v14.FRAME

def rotation(a,b):
 v=np.cross(a,b);c=float(a@b);assert c>-1+1e-10;x,y,z=v;K=np.array([[0,-z,y],[z,0,-x],[-y,x,0]]);return np.eye(3)+K+K@K/(1+c)

def main():
 began=time.monotonic();assert digest_file(CURRENT)==EXPECTED;assert not OUT.exists();OUT.mkdir();a=json.loads(CURRENT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();old=np.load(WORK/'admissible-departure-v14/sections.npz');old_outer=old['outer_sections'];old_centers=old['centers'];N=len(old_centers);first,last=11,21;variable=np.arange(first,last+1);fixed=np.setdiff1d(np.arange(N),variable);tangent=np.gradient(old_centers,axis=0);tangent/=np.linalg.norm(tangent,axis=1)[:,None];ref=np.array([0.,0.,1.]);x=ref-tangent*(tangent@ref)[:,None];x/=np.linalg.norm(x,axis=1)[:,None];y=np.cross(tangent,x);values=np.arange(-10,11)*.0005;grid=np.array([[u,v]for u in values for v in values]);K=len(grid);zero=int(np.argmin(np.linalg.norm(grid,axis=1)));radius=np.linalg.norm(old_outer-old_centers[:,None],axis=2).max(1);candidates=old_centers[:,None]+grid[None,:,0,None]*x[:,None]+grid[None,:,1,None]*y[:,None];safe=np.ones((N,K),bool)
 blockers=[('TA98HEPATIC:male:biliary',reader.read_part(a,parts['TA98HEPATIC:male:biliary'],frame=FRAME,catalogue_path=CURRENT)),('TA98CYSTICBLOOD:male:vein',reader.read_part(a,parts['TA98CYSTICBLOOD:male:vein'],frame=FRAME,catalogue_path=CURRENT)),('TA98HEPATIC:male:portal',reader.read_part(a,parts['TA98HEPATIC:male:portal'],frame=FRAME,catalogue_path=CURRENT)),('TA98HEPATIC:male:hepatic_vein',reader.read_part(a,parts['TA98HEPATIC:male:hepatic_vein'],frame=FRAME,catalogue_path=CURRENT)),('GB-innermost',Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz'))]
 source_guards=[];query_evidence=[]
 for name,m in blockers:
  if name in parts:source_guards.append(dict(id=name,geometrySha256=m.geometry_sha,bankSha256=m.source['bank_sha256']))
  index=SurfaceIndex(m);hit,ev=index.closest_points(candidates[variable].reshape(-1,3),certify=False);distance=hit['distance_m'].reshape(len(variable),K);safe[variable]&=distance>=radius[variable,None]+.0004;query_evidence.append(dict(id=name,**ev));atomic_npz(OUT/(name.replace(':','_')+'-lattice-distances.npz'),distances_m=distance)
 safe[fixed]=False;safe[fixed,zero]=True;transition=cdist(grid,grid);cost_of_step=3*(transition/.0005)**2;cost_of_step[transition>.00101]=np.inf;cost=np.full(K,np.inf);cost[zero]=0;predecessors=np.full((N,K),-1,int)
 for j in variable:
  choices=cost[:,None]+cost_of_step;predecessors[j]=choices.argmin(0);cost=choices.min(0)+np.sum((grid/.001)**2,axis=1)*.03;cost[~safe[j]]=np.inf
  if not np.isfinite(cost).any():
   proof=dict(schema='local-cystic-artery-corridor/v35',stage=1,maximumStages=3,sourceCatalogueSha256=EXPECTED,firstDisconnectedStation=int(j),safeCounts=safe[variable].sum(1).tolist(),queryEvidence=query_evidence,stoppedBeforeConstruction=True,installed=False,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'gates.json',proof);print(json.dumps(dict(disconnectedStation=int(j),elapsedSeconds=proof['elapsedSeconds'])));return
 cost+=cost_of_step[:,zero];assert np.isfinite(cost).any(),'Connected corridor does not return to fixed ring22';selected=np.zeros(N,int);selected[last]=int(cost.argmin())
 for j in range(last,first,-1):selected[j-1]=predecessors[j,selected[j]]
 target=candidates[np.arange(N),selected]-old_centers;target[fixed]=0;D=np.diff(np.eye(N),n=2,axis=0);A=np.r_[np.eye(N),8*D];delta=np.zeros((N,3));success=[]
 for d in range(3):
  rhs=np.r_[target[:,d]*1000,np.zeros(N-2)];center=target[variable,d]*1000;fit=lsq_linear(A[:,variable],rhs,bounds=(center-.1,center+.1),max_iter=60,tol=1e-9);delta[variable,d]=fit.x/1000;success.append(bool(fit.success))
 centers=old_centers+delta;outer=old_outer.copy();newtangent=np.gradient(centers,axis=0);newtangent/=np.linalg.norm(newtangent,axis=1)[:,None]
 for j in variable:outer[j]=(old_outer[j]-old_centers[j])@rotation(tangent[j],newtangent[j]).T+centers[j]
 fluid=old['fluid_sections'].copy()
 for j in variable:fluid[j]=(old['fluid_sections'][j]-old_centers[j])@rotation(tangent[j],newtangent[j]).T+centers[j]
 assert np.array_equal(outer[fixed],old_outer[fixed])and np.array_equal(fluid[fixed],old['fluid_sections'][fixed]);ov,of,_=v10.swept_profiles(outer);iv,iff,_=v10.swept_profiles(fluid);side=of[:-36];fn=np.cross(ov[side[:,1]]-ov[side[:,0]],ov[side[:,2]]-ov[side[:,0]]);radial=ov[side].mean(1)-np.repeat(centers,36,axis=0)[:len(side)]
 if np.median(np.sum(fn*radial,axis=1))<0:side=side[:,::-1]
 nominal=centers[:,None]+.82*(outer-centers[:,None]);prism=StationPrismConstraints(np.array([outer.reshape(-1,3),nominal.reshape(-1,3)]),side,np.repeat(np.arange(N),18),np.ones(N,bool),allow_source_faults=True);bad=np.flatnonzero(~prism.sourceHealthy[0]);invalid=np.unique(np.repeat(np.arange(N),18)[side[bad]].min(1)).tolist();area=lambda p:np.array([projected_area(x)for x in p]);diam=lambda p:np.linalg.norm(p[:,:,None]-p[:,None,:],axis=3).max((1,2));ae=float(abs(area(outer)/area(old_outer)-1).max());de=float(abs(diam(outer)/diam(old_outer)-1).max());bend=bend_diagnostics(centers,outer[None]);atomic_npz(OUT/'sections.npz',outer_sections=outer,fluid_sections=fluid,centers=centers,source_outer_sections=old_outer,source_fluid_sections=old['fluid_sections'],fixed_stations=fixed,station_displacements_m=delta)
 proof=dict(schema='local-cystic-artery-corridor/v35',stage=1,maximumStages=3,sourceCatalogueSha256=EXPECTED,sourceStage='admissible-departure-v14',movableRings=variable.tolist(),fixedRings=fixed.tolist(),first8RootSectionsLiteral=True,terminalAnnulusAndGBEntryLiteral=True,sourceSectionAreaRelativeError=ae,sourceSectionDiameterRelativeError=de,sourceCalibreErrorBound=.03756634960527+de,sourceCalibreTolerance=.15,continuousPrismInvalidIntervals=invalid,maximumContinuousPrismDeterminant=float(prism.sourceMaximum.max()),maximumCenterDisplacement_m=float(np.linalg.norm(delta,axis=1).max()),maximumStationTurnDegrees=float(np.rad2deg(bend['station_turn_rad']).max()),maximumCurvatureSupport=float(bend['curvature_support'].max()),smoothSolveSuccess=all(success),fixedSourceGuards=source_guards,queryEvidence=query_evidence,installed=False,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-began)
 if invalid or ae>.15 or de>.15:
  proof['stoppedBeforeCSG']=True;atomic_json(OUT/'gates.json',proof);print(json.dumps(dict(stoppedBeforeCSG=True,invalidPrisms=invalid,elapsedSeconds=proof['elapsedSeconds'])));return
 host=Mesh.load(WORK/'fixed-target-20.npz');blood=Mesh.load(WORK/'native-parent-envelopes-v9/artery-blood.npz');lineage=json.loads((WORK/'native-parent-envelopes-v9/artery-inner-lineage.json').read_text());root=centers[0];_,ho,hb,residual=v10.native_outer(host,blood,lineage,root);bo=v10.previous.solid(ov,of,root);bb=v10.previous.solid(iv,iff,root);whole=(ho+bo)-(hb+bb);network=Mesh(*v10.previous.arrays(whole,root),FRAME);network.save(OUT/'artery-joined-wall.npz');branch=Mesh(*v10.previous.arrays((bo-hb)-bb,root),FRAME);branch.save(OUT/'artery-branch-wall.npz');top=topology(network);blood_top=topology(Mesh(*v10.previous.arrays(hb+bb,root),FRAME));escape=[]
 for i,solid in enumerate((bb-bo).decompose()):
  p,f=v10.previous.arrays(solid,root);nearest=int(np.linalg.norm(centers-p.mean(0),axis=1).argmin());volume=float(solid.volume())/1e9
  if nearest<N-2 and volume>1e-23:escape.append(dict(component=i,nearestStation=nearest,volume_m3=volume))
 checks=[]
 for name,context in blockers+[(pid,reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CURRENT))for pid in ['FJ3080','TA98PART:FJ3079:A05.8.01.061']]:
  data,check=pair_contacts(branch,context,0.,max_pairs=500000,max_tests=5000000);assert check['complete_within_pair_budget'];atomic_npz(OUT/(name.replace(':','_')+'-contacts.npz'),**data);checks.append(dict(id=name,geometrySha256=context.geometry_sha,**check))
 overlap=(hb^bb).decompose();passed=not escape and all(c['crossing_face_pairs']==0 and c['coplanar_overlap_face_pairs']==0 for c in checks)and top['closed_oriented_manifold']and len(top['components'])==1 and top['components'][0]['closed_orientable_genus']==1 and len(overlap)==1 and len(blood_top['components'])==1 and blood_top['components'][0]['closed_orientable_genus']==0;proof.update(joinedTopology=top,composedBloodTopology=blood_top,parentBloodOverlapComponents=len(overlap),lateralBloodEscapeComponents=escape,contextChecks=checks,allRequiredGeometryGatesPassed=passed,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'gates.json',proof);assert digest_file(CURRENT)==EXPECTED;print(json.dumps(dict(passed=passed,ductCrossings=next(c['crossing_face_pairs']for c in checks if c['id']=='TA98HEPATIC:male:biliary'),fixedVeinCrossings=next(c['crossing_face_pairs']for c in checks if c['id']=='TA98CYSTICBLOOD:male:vein'),invalidPrisms=invalid,displacementMm=proof['maximumCenterDisplacement_m']*1000,maxTurnDegrees=proof['maximumStationTurnDegrees'],elapsedSeconds=proof['elapsedSeconds'])))
if __name__=='__main__':main()
