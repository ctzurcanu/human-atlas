"""Read-only exact native ACA/ACom witnesses; no invented junction mouths."""
from pathlib import Path
import sys,time,json,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.distance import triangle_pairs
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
_spec=importlib.util.spec_from_file_location('aca_readonly_levels',Path(__file__).with_name('prepare-native-aca-bounded-references-v3.py'))
_mod=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_mod);levels=_mod.levels

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1/native-aca-bounded-references-v3/exact-witnesses-v4'
FRAME='model:local-ta98:assembled-metres'

def provenance(m,p):
    return dict(partId=p['id'],sourceCatalogueSha256=EXPECTED,sourceBankSha256=m.source['bank_sha256'],sourceGeometrySha256=m.geometry_sha,sourceAttributeSha256=m.attributes_sha,license=p.get('license'),frame=m.frame)

def witness(d,j,first,second):
    af,bf=map(int,d['pairs'][j]);x=d['first'][j];y=d['second'][j]
    a=first.triangles[af];b=second.triangles[bf]
    replay=max(np.linalg.norm(d['first_barycentric'][j]@a-x),np.linalg.norm(d['second_barycentric'][j]@b-y),abs(np.linalg.norm(x-y)-d['distance_m'][j]))
    assert replay<1e-12
    return dict(firstSourceFace=af,secondSourceFace=bf,firstSourceVertices=first.faces[af].tolist(),secondSourceVertices=second.faces[bf].tolist(),firstTriangle_m=a.tolist(),secondTriangle_m=b.tolist(),firstPoint_m=x.tolist(),secondPoint_m=y.tolist(),firstBarycentric=d['first_barycentric'][j].tolist(),secondBarycentric=d['second_barycentric'][j].tolist(),distance_m=float(d['distance_m'][j]),relation=str(d['relation'][j]),normalCosine=float(d['normal_cosine'][j]),numericalRoundoff_m=float(d['numerical_roundoff_m'][j]),witnessReplayError_m=float(replay))

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED;assert not (OUT/'proof.json').exists();OUT.mkdir(exist_ok=True)
    atlas=json.loads(CAT.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    p=parts['ZA:Anterior communicating artery'];comm=reader.read_part(atlas,p,frame=FRAME,catalogue_path=CAT)
    cm=provenance(comm,p);cm.update(topology=topology(comm),actualOpenBoundaryLoops=[q.tolist()for q in boundary_loops(comm.faces)])
    assert not cm['actualOpenBoundaryLoops'];rows=[]
    for side in ['l','r']:
        p=parts['ZA:Anterior cerebral artery.'+side];m=reader.read_part(atlas,p,frame=FRAME,catalogue_path=CAT)
        ip=parts['ZA:Internal carotid artery.'+side];ica=reader.read_part(atlas,ip,frame=FRAME,catalogue_path=CAT)
        loops=boundary_loops(m.faces);ends=np.array([m.vertices[q].mean(0)for q in loops]);hit,_=BVH(ica).closest_points(ends);prox=int(np.argmin(hit['distance_m']))
        ring,centers=levels(m,loops[prox]);course=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(centers,axis=0),axis=1))]
        d,summary=pair_contacts(m,comm,.100,max_pairs=250000,max_tests=2000000)
        assert summary['complete_within_pair_budget'] and len(d['pairs']);j=int(np.argmin(d['distance_m']));w=witness(d,j,m,comm)
        af=w['firstSourceFace'];w.update(acaSourceRingRange=[int(ring[m.faces[af]].min()),int(ring[m.faces[af]].max())],acaCourseDistanceRangeFromIcaEnd_m=[float(course[ring[m.faces[af]].min()]),float(course[ring[m.faces[af]].max()])])
        # Every triangle pair closer than 100 mm was exhaustively queried. Since
        # the returned minimum is smaller, no excluded pair can improve it.
        assert w['distance_m']<.100
        looprows=[]
        for n,q in enumerate(loops):
            points=m.vertices[q];edges=np.stack([points,np.roll(points,-1,axis=0),points],axis=1)
            aa=np.repeat(edges,len(comm.faces),axis=0);bb=np.tile(comm.triangles,(len(q),1,1))
            ed=triangle_pairs(aa,bb);z=int(np.argmin(ed['distance_m']));ei=z//len(comm.faces);cf=z%len(comm.faces)
            looprows.append(dict(sourceBoundaryLoop=n,sourceVertexIds=q.tolist(),sourcePoints_m=points.tolist(),sourceRing=int(ring[q[0]]),role='ICA-proximal-open-rim'if n==prox else'distal-native-open-rim',finiteWallAnnulus=False,centroid_m=ends[n].tolist(),centroidToActualIcaSurface_m=float(hit['distance_m'][n]),rimToAComMinimum_m=float(ed['distance_m'][z]),rimEdgeSourceVertexIds=[int(q[ei]),int(q[(ei+1)%len(q)])],acomSourceFace=cf,rimWitness_m=ed['first'][z].tolist(),acomWitness_m=ed['second'][z].tolist(),acomBarycentric=ed['second_barycentric'][z].tolist(),method='All actual source boundary line segments versus all ACom triangles; a degenerate query triangle encodes each line segment, and no geometry is emitted'))
        rows.append(dict(side=side,aca=provenance(m,p),ica=provenance(ica,ip),closestActualSurfaceWitness=w,exactGlobalMinimumWithinFloat64TriangleArithmetic=True,querySummary=summary,sourceOpenRims=looprows,nativeCourseRings=len(centers),nativeCourseLength_m=float(course[-1]),sourceMouthAnnulusPairExists=False,reason='ACA has two open single-surface rims; ACom has no open source boundary and no inner wall, so a literal paired-wall annulus pair does not exist'))
    proof=dict(schema='native-aca-acom-exact-witnesses/v4',sourceCatalogue=str(CAT),sourceCatalogueSha256=EXPECTED,acom=cm,rows=rows,noAnatomicalBoundaryInvented=True,noA1A2ControlsEmitted=True,noGeometryOrLiveMetadataChanged=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',proof);assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(output=str(OUT/'proof.json'),acomOpenRims=0,rows=[dict(side=r['side'],gap_mm=1000*r['closestActualSurfaceWitness']['distance_m'],faces=[r['closestActualSurfaceWitness']['firstSourceFace'],r['closestActualSurfaceWitness']['secondSourceFace']],acaRingRange=r['closestActualSurfaceWitness']['acaSourceRingRange'],rimGaps_mm=[1000*q['rimToAComMinimum_m']for q in r['sourceOpenRims']])for r in rows],elapsedSeconds=proof['elapsedSeconds'])))

if __name__=='__main__':main()
