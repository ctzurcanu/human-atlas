"""Recover actual installed duct envelopes from literal Boolean face ancestry."""
from pathlib import Path
import sys,time,json
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import topology
from geometry_tools.contact import self_contact_summary
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,oriented_key,same_oriented_triangles

W=ROOT/'.local-models/geometry-workbench/liver-completion'
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
PRE=B/'biliary-actual-boundary-recovery-v38';OUT=B/'biliary-actual-envelope-v39'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345'

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED;assert not OUT.exists();OUT.mkdir()
    a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};actual=AtlasReader().read_part(a,parts['TA98HEPATIC:male:biliary'],frame='model:local-ta98:assembled-metres',catalogue_path=CAT)
    source=Mesh.load(W/'packed-hepatic-biliary-v128/TA98HEPATIC_male_biliary.npz');old=Mesh.load(W/'rebuilt-biliary-v107/wall.npz');z=np.load(PRE/'ancestry-and-roi.npz');oldrole=z['original_wall_boolean_role'];origin=np.array([0.,1.25,0.])
    quantized=(old.vertices-origin).astype('<f4').astype(float)+origin
    d,vertexmap=cKDTree(quantized).query(source.vertices);assert np.all(d==0),'Packed source contains a vertex outside exact original Float32 controls'
    keys={oriented_key(f):int(r)for f,r in zip(old.faces,oldrole)}
    sourcekeys=[oriented_key(f)for f in vertexmap[source.faces]];missing=[i for i,k in enumerate(sourcekeys)if k not in keys]
    atomic_json(OUT/'source-facet-replay.json',dict(sourceFaces=len(source.faces),originalFaces=len(old.faces),exactQuantizedVertexError_m=float(d.max()),unmatchedSourceFaces=missing,BooleanReplaySha256=digest_file(PRE/'proof.json')))
    assert not missing,'Float32 repair changed original facet topology: requires explicit further ancestry'
    role=np.array([keys[k]for k in sourcekeys],np.int8)
    assert np.array_equal(role,z['source_corner_role']),'Independent Boolean and corner role classifications differ'
    p=np.load(W/'hepatic-contained-banks-v181/biliary-parents.npz');parent=p['parent'];currole=role[parent];evidence=json.loads((W/'hepatic-contained-banks-v181/biliary-evidence.json').read_text());assert actual.geometry_sha==evidence['geometrySha256']
    # A source-side diagonal may change parent facet but must never bridge the
    # outer and inner tissue boundaries. Verify every recorded diagonal pair.
    flips=[]
    for change in evidence['sliverTriangulation']['changes']:
        ids=np.array(change['faceIds']);labels=currole[ids];assert len(set(labels))==1,'Transport diagonal joins opposite original sheets';flips.append(dict(faceIds=ids.tolist(),sheet=int(labels[0])))
    known=z['current_corner_role']>=0;assert np.all(z['current_corner_role'][known]==currole[known])
    positions,inverse=np.unique(actual.vertices,axis=0,return_inverse=True);outer=inverse[actual.faces[currole==0]];inner=inverse[actual.faces[currole==1]]
    ol=boundary_loops(outer);il=boundary_loops(inner);assert len(ol)==len(il)==1
    assert set(ol[0])==set(il[0]),'Outer and inner duct mouth do not share the literal cut boundary'
    ring=ol[0];points=np.r_[positions,positions[ring].mean(0)[None]];center=len(positions)
    caps=np.array([[center,int(b),int(a)]for a,b in zip(ring,np.roll(ring,-1))])
    of=np.r_[outer,caps];iff=np.r_[inner[:,::-1],caps]
    # The closed exterior minus the closed fluid must reconstruct the wall,
    # with the exact same computational cap cancelling algebraically.
    assert same_oriented_triangles(np.r_[of,iff[:,::-1]],np.r_[inverse[actual.faces],caps,caps[:,::-1]])
    envelopes=[];checks=[]
    for name,faces in [('outer',of),('fluid',iff)]:
        mesh=Mesh(points,faces,actual.frame);top=topology(mesh)
        assert top['closed_oriented_manifold']and len(top['components'])==1 and top['components'][0]['closed_orientable_genus']==0,(name,top)
        selfcheck,_=self_contact_summary(mesh,max_pairs=1000000,max_tests=8000000,stop_on_defect=True)
        assert selfcheck['complete_within_pair_budget']and selfcheck['crossing_face_pairs']==selfcheck['coplanar_overlap_face_pairs']==0,(name,selfcheck)
        mesh.save(OUT/('actual-biliary-'+name+'-computational-envelope.npz'));envelopes.append(mesh);checks.append(dict(name=name,geometrySha256=mesh.geometry_sha,topology=top,selfIntersections=selfcheck))
    os,_=make_solid((points-origin)*1000,of);fs,_=make_solid((points-origin)*1000,iff);ws,_=make_solid((actual.vertices-origin)*1000,actual.faces)
    escape=float((fs-os).volume());difference=os-fs;symmetric=float(((difference-ws)+(ws-difference)).volume());assert escape==0 and symmetric<1e-9,(escape,symmetric)
    atomic_npz(OUT/'source-boundary-ownership.npz',source_wall_roles=role,actual_wall_roles=currole,actual_outer_face_ids=np.flatnonzero(currole==0),actual_inner_face_ids=np.flatnonzero(currole==1),computational_mouth_cap_faces=caps,positions=points)
    proof=dict(schema='actual-biliary-envelope-from-exact-Boolean-ancestry/v39',sourceCatalogueSha256=EXPECTED,actualWallGeometrySha256=actual.geometry_sha,actualBankSha256=actual.source['bank_sha256'],originalBooleanReplayAll31794FacesLiteral=True,packedSourceAll31794FacetsExactlyMatchQuantizedOriginalBooleanTriangles=True,independent100nmCornerClassificationAgrees=True,actualFaceParentMapSha256=digest_file(W/'hepatic-contained-banks-v181/biliary-parents.npz'),allRecordedSliverFlipsStayOnSameOriginalBoundary=flips,actualOuterFaces=len(outer),actualInnerFaces=len(inner),computationalMouthCapTriangles=len(caps),outerInnerMouthCoordinatesSharedExactly=True,numericalCapsCancelExactlyFromWall=True,fluidOutsideOuterEnvelope_mm3=escape,envelopeDifferenceVsActualWallSymmetricVolume_mm3=symmetric,envelopes=checks,wholeGlobalEnvelope=True,localSearchRoiNotUsedToOmitAnyFace=True,emittedDuctCaps=0,geometryModified=False,installed=False,passed=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',proof);assert digest_file(CAT)==EXPECTED;print(json.dumps(dict(passed=True,outerFaces=len(outer),innerFaces=len(inner),capFaces=len(caps),symmetricWallDifference_mm3=symmetric,elapsedSeconds=proof['elapsedSeconds'])))
if __name__=='__main__':main()
