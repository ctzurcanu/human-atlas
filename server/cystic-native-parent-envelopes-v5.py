"""Recover the actual single-root hepatic wall sheets, without moving points."""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.algorithms.annular_pipe_envelopes import face_adjacency,face_components,boundary_loops,projected_area,close_with_shared_disks,signed_volume

WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
OUT=WORK/'native-parent-envelopes-v9'

def recover_single_root(mesh,parent,source_interface):
    v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces];pairs,count=face_adjacency(f);assert np.all(count==2)
    tri=v[f];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    similar=np.einsum('ij,ij->i',normal[pairs[:,0]],normal[pairs[:,1]])>=np.cos(.04)
    n,labels=face_components(f,pairs[similar]);tol=2*np.linalg.norm(np.spacing(v.astype('<f4')).astype(float),axis=1).max();caps=[]
    for label in np.flatnonzero(np.bincount(labels)>=12):
        ids=np.flatnonzero(labels==label);ff=f[ids]
        try:loops=boundary_loops(ff)
        except ValueError:continue
        if len(loops)!=2:continue
        points=v[np.unique(ff)];_,_,basis=np.linalg.svd(points-points.mean(0),full_matrices=False)
        residual=float(abs((points-points.mean(0))@basis[-1]).max())
        if residual>tol:continue
        loops.sort(key=lambda x:projected_area(v[x]));caps.append(dict(faces=ids,inner=loops[0],outer=loops[1],planeResidual_m=residual))
    literal=set(map(tuple,parent.vertices));scores=[sum(tuple(p)in literal for p in v[c['outer']])for c in caps]
    atomic_json(OUT/('annular-candidates-'+('artery'if'Coeliac'in parent.source.get('record',{}).get('id','')else'vein')+'.json'),dict(candidates=[dict(rootCenter_m=v[c['outer']].mean(0).tolist(),outerControls=len(c['outer']),literalSharedParentControls=score,faces=len(c['faces']))for c,score in zip(caps,scores)]))
    # Saved native parent construction records identify the actual root
    # faces even where the transported annulus is not planar. Verify those
    # literal indexed faces rather than choosing an annular-looking facet.
    lineage=json.loads(source_interface.read_text());capids=np.array(lineage['network_interface_faces'],int)
    loops=boundary_loops(f[capids]);assert len(loops)==2;loops.sort(key=lambda x:projected_area(v[x]))
    cap=dict(faces=capids,inner=loops[0],outer=loops[1],planeResidual_m=None)
    points=v[np.unique(f[capids])];_,_,basis=np.linalg.svd(points-points.mean(0),full_matrices=False);cap['planeResidual_m']=float(abs((points-points.mean(0))@basis[-1]).max())
    saved=reader=AtlasReader();bank=Path(lineage['bank']);record=dict(lineage['part'],chunk=0)
    saved_mesh=reader.read_part(dict(chunks=[dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')))]),record,frame=mesh.frame)
    source_keys={tuple(sorted(map(tuple,tri)))for tri in v[f[capids]]}
    saved_keys={tuple(sorted(map(tuple,tri)))for tri in saved_mesh.triangles[np.array(lineage['parent_interface_faces'],int)]}
    from scipy.spatial import cKDTree
    saved_tri=saved_mesh.triangles[np.array(lineage['parent_interface_faces'],int)]
    actual_tri=v[f[capids]];tree=cKDTree(saved_tri.reshape(-1,3));distance,nearest=tree.query(actual_tri.reshape(-1,3))
    maximum_root_encoding_residual=float(distance.max());assert maximum_root_encoding_residual<5e-9,'Saved root faces differ beyond their explicit Float32 origin-encoding allowance'
    replay=tuple(sorted(tuple(sorted(map(int,row)))for row in nearest.reshape(-1,3)))
    _,saved_inverse=np.unique(saved_tri.reshape(-1,3),axis=0,return_inverse=True)
    # A nearest corner must identify the corresponding source triangle;
    # coincident surface vicinity cannot substitute for face ancestry.
    saved_tree=cKDTree(np.unique(saved_tri.reshape(-1,3),axis=0));actual_ids=saved_tree.query(actual_tri.reshape(-1,3))[1].reshape(-1,3)
    assert {tuple(sorted(t))for t in actual_ids}=={tuple(sorted(t))for t in saved_inverse.reshape(-1,3)}
    sideids=np.setdiff1d(np.arange(len(f)),cap['faces']);nc,ls=face_components(f[sideids]);assert nc==2
    surfaces=[]
    for k in range(nc):
        ids=sideids[ls==k];used=set(f[ids].ravel());kind='inner'if set(cap['inner'])<=used else'outer';assert set(cap[kind])<=used
        positions,closed,disks=close_with_shared_disks(v,f[ids],expected_mouths=1)
        volume=signed_volume(positions,closed)
        if volume<0:closed=closed[:,::-1]
        surfaces.append((kind,Mesh(positions,closed,mesh.frame,source=dict(kind='actual native '+kind+' hepatic sheet plus temporary root classifier disk',parentGeometrySha256=mesh.geometry_sha)),dict(kind=kind,sourceFaceIds=ids.tolist(),signedNativeVolume_m3=volume,temporaryRootDiskFaces=len(disks))))
    assert {x[0]for x in surfaces}=={'inner','outer'}
    assert all((e['signedNativeVolume_m3']<0)==(kind=='inner')for kind,m,e in surfaces),'Native inner/outer orientation does not establish fluid cavity'
    assert all(len(e['sourceFaceIds'])>len(f)*.1 for kind,m,e in surfaces),'A local incidental patch is not the whole native blood sheet'
    return surfaces,dict(annularRootFaceIds=cap['faces'].tolist(),annularRootLoops=dict(inner=v[cap['inner']].tolist(),outer=v[cap['outer']].tolist()),rootPlaneResidual_m=cap['planeResidual_m'],sourceCoordinatesLiteral=True,oneRootAndTwoActualSourceSheets=True,namedUpstreamParentInterfaceEstablished=False,namedParentLiteralMatchControls=max(scores),actualRootSourceLineage=str(source_interface.relative_to(ROOT)),actualRootSourceLineageSha256=digest_file(source_interface),actualCurrentRootFacesSameSavedNativeAncestry=True,maximumSavedInterfaceEncodingResidual_m=maximum_root_encoding_residual)

def main():
    start=time.monotonic();assert not OUT.exists();OUT.mkdir();rows=[]
    cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';atlas=json.loads(cat.read_text());reader=AtlasReader()
    base=ROOT/'.local-models/geometry-workbench/liver-completion'
    for role,index,pid,interface in [('artery',20,'ZA:Coeliac trunk',base/'male-parent-walls-v230/artery-result.json'),('vein',21,'TA98HEPATIC:male:portal-distal',base/'male-parent-walls-v252/portal-result.json')]:
        part=next(p for p in atlas['parts']if p['id']==pid);parent=reader.read_part(atlas,part,frame='model:local-ta98:assembled-metres',catalogue_path=cat)
        mesh=Mesh.load(WORK/('fixed-target-'+str(index)+'.npz'));surfaces,proof=recover_single_root(mesh,parent,interface)
        for kind,surface,lineage in surfaces:
            surface.save(OUT/(role+'-'+('blood'if kind=='inner'else'outer')+'.npz'));atomic_json(OUT/(role+'-'+kind+'-lineage.json'),lineage)
        rows.append(dict(role=role,parentGeometrySha256=mesh.geometry_sha,proof=proof,surfaces=[dict(kind=k,volume_m3=abs(e['signedNativeVolume_m3']))for k,m,e in surfaces]))
    atomic_json(OUT/'proof.json',dict(rows=rows,elapsedSeconds=time.monotonic()-start));print(json.dumps(dict(output=str(OUT),rows=[dict(role=r['role'],surfaces=r['surfaces'])for r in rows],elapsedSeconds=time.monotonic()-start)))

if __name__=='__main__':main()
