"""Read-only original-boundary/corner audit, with no tolerance enlargement."""
from pathlib import Path
import sys,time,json
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance

W=ROOT/'.local-models/geometry-workbench/liver-completion'
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
OUT=B/'biliary-actual-boundary-recovery-v38'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345'

def query(points,meshes):
    return np.c_[*[SurfaceIndex(m).closest_points(points,certify=False)[0]['distance_m']for m in meshes]]

def corner_class(vertices,faces,dist):
    support=dist<=1e-7
    allside=support[faces].all(1)
    roles=np.full(len(faces),-1,np.int8)
    unique=allside.sum(1)==1
    roles[unique]=allside[unique].argmax(1)
    return roles,allside

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED;assert not OUT.exists();OUT.mkdir()
    a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']}
    actual=AtlasReader().read_part(a,parts['TA98HEPATIC:male:biliary'],frame='model:local-ta98:assembled-metres',catalogue_path=CAT)
    source=Mesh.load(W/'packed-hepatic-biliary-v128/TA98HEPATIC_male_biliary.npz')
    original=[Mesh.load(W/'rebuilt-biliary-v107'/s)for s in ['outer.npz','fluid.npz']]
    oldwall=Mesh.load(W/'rebuilt-biliary-v107/wall.npz')
    p=np.load(W/'hepatic-contained-banks-v181/biliary-parents.npz');parent=p['parent'];sv=p['source_vertices']
    assert len(parent)==len(actual.faces)and len(sv)==len(actual.vertices)
    centroid=query(source.triangles.mean(1),original);distance=query(source.vertices,original)
    role,support=corner_class(source.vertices,source.faces,distance)
    curdistance=query(sv,original);currole,cursupport=corner_class(sv,actual.faces,curdistance)
    inherited=role[parent]
    same=currole==inherited;bothknown=(currole>=0)&(inherited>=0)
    oldbad=centroid.min(1)>=1e-7;oldbad_current=oldbad[parent]
    z=np.load(B/'admissible-departure-v14/sections.npz');outer=z['outer_sections'];centers=z['centers']
    # The whole ±5mm transverse search domain is enclosed by this conservative
    # axis-aligned box, including source profile radius plus 0.25mm clearance.
    radius=np.linalg.norm(outer-centers[:,None],axis=2).max(1)
    margin=np.sqrt(2)*.005+radius[8:22].max()+.00025
    roi=np.array([centers[8:22].min(0)-margin,centers[8:22].max(0)+margin])
    tri=actual.triangles;touch=np.all(tri.max(1)>=roi[0],axis=1)&np.all(tri.min(1)<=roi[1],axis=1)
    # Re-evaluate only the original saved Boolean, tagging exact source roles.
    origin=np.array([0.,1.25,0.]);os,oid=make_solid((original[0].vertices-origin)*1000,original[0].faces)
    fs,fid=make_solid((original[1].vertices-origin)*1000,original[1].faces)
    rebuilt=(os-fs).to_mesh64();origins,faceids=provenance(rebuilt)
    rv=np.asarray(rebuilt.vert_properties)[:,:3]/1000+origin;rf=np.asarray(rebuilt.tri_verts)
    rr=np.where(origins==oid,0,np.where(origins==fid,1,-1)).astype(np.int8);assert np.all(rr>=0)
    t=rv[rf];tree=cKDTree(t.mean(1));cd,nearest=tree.query(oldwall.triangles.mean(1))
    # All corners must literally agree (within 1pm arithmetic replay error),
    # with orientation kept; a nearest centroid is never enough.
    target=t[nearest];original_tri=oldwall.triangles
    errors=np.stack([np.linalg.norm(original_tri-np.roll(target,k,axis=1),axis=2).max(1)for k in range(3)],axis=1).min(1)
    exact=errors<=1e-12
    replayrole=np.full(len(oldwall.faces),-1,np.int8);replayrole[exact]=rr[nearest[exact]]
    atomic_npz(OUT/'ancestry-and-roi.npz',source_vertex_boundary_distances_m=distance,source_centroid_boundary_distances_m=centroid,source_corner_role=role,current_source_corner_boundary_distances_m=curdistance,current_corner_role=currole,current_parent_role=inherited,current_parent=parent,previously_unclassified_source_faces=np.flatnonzero(oldbad),previously_unclassified_current_faces=np.flatnonzero(oldbad_current),current_faces_intersecting_conservative_search_roi=np.flatnonzero(touch),conservative_search_roi_m=roi,original_wall_boolean_role=replayrole,original_wall_boolean_replay_corner_error_m=errors)
    proof=dict(schema='actual-biliary-boundary-recovery/v38',sourceCatalogueSha256=EXPECTED,actualWallGeometrySha256=actual.geometry_sha,actualWallBankSha256=actual.source['bank_sha256'],sourceWallGeometrySha256=source.geometry_sha,boundaryWitnessTolerance_m=1e-7,toleranceNotEnlarged=True,sourceCentroidUnclassified=int(oldbad.sum()),currentCentroidUnclassified=int(oldbad_current.sum()),sourceCornerClassified=int((role>=0).sum()),sourceCornerUnclassified=int((role<0).sum()),currentCornerClassified=int((currole>=0).sum()),currentCornerUnclassified=int((currole<0).sum()),currentCornerVsParentKnownDisagreements=int((~same&bothknown).sum()),previouslyUnclassifiedCurrentFacesInConservativeSearchRoi=int((oldbad_current&touch).sum()),currentCornerUnclassifiedInSearchRoi=int(((currole<0)&touch).sum()),conservativeSearchRoi_m=roi.tolist(),roiIncludesAllOriginalRings8To21PlusEvery5mmTransverseLatticePositionAndSourceProfileRadius=True,booleanReplay=dict(originalWallFaces=len(oldwall.faces),replayedFaces=len(rf),exactCornerMatchedFaces=int(exact.sum()),maximumMatchedCornerError_m=float(errors[exact].max())if exact.any()else None,unmatchedFaces=np.flatnonzero(~exact).tolist(),outerFaces=int((replayrole==0).sum()),innerFaces=int((replayrole==1).sum()),replayIsReadOnlyNotAReplacement=True),geometryModified=False,installed=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',proof);assert digest_file(CAT)==EXPECTED;print(json.dumps(proof))
if __name__=='__main__':main()
