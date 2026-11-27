"""Reconstruct a proposal-only outer classifier from the literal installed duct wall."""
import numpy as np,json
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.audit import topology
from geometry_tools.contact import self_contact_summary
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

def reconstruct(actual,out):
 w=ROOT/'.local-models/geometry-workbench/liver-completion';source=Mesh.load(w/'packed-hepatic-biliary-v128/TA98HEPATIC_male_biliary.npz');field=json.loads((w/'hepatic-contained-banks-v181/biliary-evidence.json').read_text());assert actual.geometry_sha==field['geometrySha256'];parents=np.load(w/'hepatic-contained-banks-v181/biliary-parents.npz')['parent'];assert len(parents)==len(actual.faces);q=source.triangles.mean(1);distance=[]
 for filename in ['outer.npz','fluid.npz']:
  mesh=Mesh.load(w/'rebuilt-biliary-v107'/filename);hit,evidence=SurfaceIndex(mesh).closest_points(q,certify=False);distance.append(hit['distance_m'])
 d=np.array(distance).T;which=d.argmin(1);assert np.all(d.min(1)<1e-7),'Source boundary lineage requires an unresolved side classification';assert np.all(abs(d[:,0]-d[:,1])>1e-8),'Ambiguous original tissue boundary';faces=actual.faces[which[parents]==0];used,inverse=np.unique(actual.vertices,axis=0,return_inverse=True);faces=inverse[faces];loops=boundary_loops(faces);assert len(loops)==1,'Actual outer duct surface does not have one documented upstream mouth';directed={tuple(e)for e in np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]].tolist()};positions=used.tolist();caps=[]
 for ring in loops:
  ring=list(map(int,ring));center=len(positions);positions.append(used[ring].mean(0).tolist())
  for a,b in zip(ring,ring[1:]+ring[:1]):caps.append([center,b,a]if(a,b)in directed else[center,a,b])
 volume=Mesh(np.array(positions),np.r_[faces,np.array(caps)],actual.frame);top=topology(volume);assert top['closed_oriented_manifold']and len(top['components'])==1 and top['components'][0]['closed_orientable_genus']==0;selfcheck,_=self_contact_summary(volume,max_pairs=500000,max_tests=5000000,stop_on_defect=True);assert selfcheck['complete_within_pair_budget']and selfcheck['crossing_face_pairs']==0 and selfcheck['coplanar_overlap_face_pairs']==0;volume.save(out/'actual-biliary-outer-computational-envelope.npz');proof=dict(schema='actual-installed-biliary-proposal-envelope/v1',actualWallGeometrySha256=actual.geometry_sha,actualWallBankSha256=actual.source['bank_sha256'],sourceWallGeometrySha256=source.geometry_sha,sourceParentsSha256=digest_file(w/'hepatic-contained-banks-v181/biliary-parents.npz'),originalBoundaryWitnessMaximumResidual_m=float(d.min(1).max()),sourceBoundaryClassification='Literal native face ancestry and actual original outer/fluid triangle witnesses; all source faces classified',outerSourceFaces=len(faces),computationalUpstreamMouthCapTriangles=len(caps),capsAreNumericalNotAnatomical=True,proposalsOnlyNotAcceptance=True,topology=top,selfIntersections=selfcheck,geometrySha256=volume.geometry_sha,installed=False);atomic_json(out/'actual-biliary-envelope-proof.json',proof);return volume,proof


def vessel_envelope(actual,source_partition,out,*,label="fixed-vein",bank_sha256=None):
 z=np.load(source_partition);ids=z['branch_source_face_ids'];shared=z['shared_branch_faces'];take=np.ones(len(ids),bool);take[shared]=False;order=[]
 for sheet in sorted(set((ids//1000000)[take])):order.extend(np.flatnonzero(take&(ids//1000000==sheet)))
 order.extend(shared);order=np.array(order,int);assert len(order)==len(actual.faces);mask=np.isin(ids[order]//1000000,[1,6]);faces=actual.faces[mask];positions,inverse=np.unique(actual.vertices,axis=0,return_inverse=True);faces=inverse[faces];loops=boundary_loops(faces);assert len(loops)==2,'Expected two actual cystic-vein blood ports';directed={tuple(e)for e in np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]].tolist()};p=positions.tolist();caps=[]
 for ring in loops:
  ring=list(map(int,ring));center=len(p);p.append(positions[ring].mean(0).tolist())
  for a,b in zip(ring,ring[1:]+ring[:1]):caps.append([center,b,a]if(a,b)in directed else[center,a,b])
 result=Mesh(np.array(p),np.r_[faces,np.array(caps)],actual.frame);top=topology(result);assert top['closed_oriented_manifold']and len(top['components'])==1 and top['components'][0]['closed_orientable_genus']==0;result.save(out/(label+'-actual-outer-computational-envelope.npz'));return result,dict(actualVeinGeometrySha256=actual.geometry_sha,actualVeinBankSha256=bank_sha256 or actual.source['bank_sha256'],actualOuterAndAnnulusFaces=len(faces),numericalBloodPortCapTriangles=len(caps),capsAreOracleOnlyNeverEmitted=True,topology=top,geometrySha256=result.geometry_sha)
