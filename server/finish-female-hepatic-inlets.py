"""Join literal native female inlet rims to the uniformly fitted donor liver.

The native foregut/cystic junctions remain the anchors; obsolete intraparenchymal
source branches are replaced, not superimposed. Each result is independently
audited after the actual local Float32 encoding.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.annular_bridge import attach_bridge,trim_named_inlet,find_annular_patch
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank,split_closed_cell
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces
ap=argparse.ArgumentParser();ap.add_argument('--role',choices=['artery','portal','biliary'],required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic()
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/args.stage;out.mkdir(exist_ok=False);role=args.role;frame='construction:female-liver:unscaled-donor-metres'
fp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';mp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(fp.read_text());ma=json.loads(mp.read_text());r=AtlasReader();parts={p['id']:p for p in a['parts']};male={p['id']:p for p in ma['parts']}
fit=json.loads((b/'final-context-v217/fit.json').read_text())['chosen'];old=a['liverTransferIntegration'];scale=old['scale']*fit['additional_scale'];translation=np.array(old['translation_m'])*fit['additional_scale']+(1-fit['additional_scale'])*np.array(fit['center_m'])+fit['translation_m']
def compact(m):
 used=np.unique(m.faces);v,inv=np.unique(m.vertices[used],axis=0,return_inverse=True);ix=np.full(len(m.vertices),-1,int);ix[used]=inv;return Mesh(v,ix[m.faces],frame)
def native(pid):
 m=r.read_part(a,parts[pid],frame='model:local-female-ta98:assembled-metres',catalogue_path=fp);return compact(Mesh((m.vertices-translation)/scale,m.faces,frame))
if role in ['artery','portal']:
 pid='LOCAL:female:proper_hepatic_artery:1718'if role=='artery'else'LOCAL:female:portal_vein:1795';source=native(pid);loops=boundary_loops(source.faces);root=min(loops,key=lambda x:source.vertices[x,1].mean());outer=source.vertices[root];center=outer.mean(0);inner=center+.85*(outer-center);pv=np.r_[outer,inner];n=len(root);pc=[]
 for i in range(n):
  j=(i+1)%n;pc.extend([[i,j,n+j],[i,n+j,n+i]])
 pf=np.asarray(pc,int);tri=pv[pf];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]).sum(0)
 if normal[1]<0:pf=pf[:,[0,2,1]]
 child=compact(r.read_part(ma,male['TA98HEPATIC:male:'+role],frame=frame,catalogue_path=mp))
 if role=='artery':
  q=json.loads((b/'male-parent-walls-v230/artery-result.json').read_text());seed=child.vertices[np.unique(child.faces[q['network_interface_faces']])].mean(0);v,f,cc,trim=trim_named_inlet(child.vertices,child.faces,1.25390625,seed);child=Mesh(v,f,frame)
 else:cc=find_annular_patch(child.vertices,child.faces,[-.0058649,1.242,.032422],[0,1,0]);trim={}
 q=attach_bridge(pv,pf,np.arange(len(pf)),child.vertices,child.faces,cc,stations=20,controls=32);q['faces']=np.r_[q['faces'],pf[:,[0,2,1]]];parent=Mesh(q['vertices'],q['faces'],frame)
 proof=dict(native_root_loop_controls=len(root),native_root_outer_literal=True,old_intrahepatic_source_replaced=True,trim=trim)
else:
 pid='LOCAL:female:bile_duct:1031';source=native(pid);child=Mesh.load(b/'common-hepatic-mouth-v235/child.npz');cc=json.loads((b/'common-hepatic-mouth-v235/junction.json').read_text())['child_interface_faces']
 # Find the common stem immediately above its actual cystic interface. The
 # superior gallbladder attachment is not the cystic/bile confluence.
 normal=np.array([-.25,1.,0.]);row=parts[pid];ids=np.r_[tuple(np.arange(g['start']//3,(g['start']+g['count'])//3)for g in row['tissueDrawGroups']if g.get('neighbor')=='LOCAL:female:cystic_duct:authored')];interface=source.triangles[ids]*scale+translation;offset=float((interface@normal).max()+.001*np.linalg.norm(normal));limit=(offset-normal@translation)/scale
 bank=SharedCutBank(source.vertices);cut=split_closed_cell(bank,source.faces,normal,limit,'native-common-hepatic',allow_holes=True);parent=Mesh(np.asarray(bank.positions),cut.negative_faces,frame);pc=np.flatnonzero(cut.negative_ancestry<0)
 q=attach_bridge(parent.vertices,parent.faces,pc,child.vertices,child.faces,cc,stations=20,controls=40,straight_shear=True);parent=Mesh(q['vertices'],q['faces'],frame);proof=dict(native_cut_normal=normal.tolist(),native_cut_offset_m=offset,native_cystic_and_gallbladder_attachments_literal=True)
origin=np.array([0,1.25,0]);surfaces,_,conform=conform_quantized_surfaces(parent.vertices,{'parent':parent.faces},origin);parent=Mesh((parent.vertices-origin).astype('<f4').astype(float)+origin,surfaces['parent'],frame);parent.save(out/'parent.npz');child.save(out/'child.npz')
check=audit(parent,require_closed=True,max_pairs=1000000,max_tests=12000000);atomic_json(out/'parent-check.json',check)
atomic_json(out/'junction.json',dict(passed=check['passed'],role=role,parent_id=pid,parent_interface_faces=q['interface_faces'].tolist(),child_interface_faces=list(map(int,cc)),curve=q['path'],lumen_disks_added=0,scale=scale,translation_m=translation.tolist(),conform=conform,proof=proof,elapsed_seconds=time.monotonic()-start));print(role,check['passed'],check['failures'],flush=True)
if not check['passed']:raise ValueError('Inlet candidate failed; not installed')
