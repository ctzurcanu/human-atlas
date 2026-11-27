"""One proper similarity of the actual HRA cystic artery/vein closure.

Correspondences are explicit geometric proposals at named anatomical sites;
neither the resulting proximity nor preserved source topology proves a patent
male recipient junction. No target geometry or catalogue is changed.
"""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.fit import landmark_transform,transform_mesh
from geometry_tools.common import ROOT,digest_file,atomic_json
from geometry_tools.spatial import BVH

WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
TARGET=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='9cf8b30456620f7b7bd6f52311ba818fc94f935beee24e6e07cb953d14106f64'

def centroid(mesh):
    tri=mesh.triangles;area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    return np.average(tri.mean(1),axis=0,weights=area)

def extremity(mesh,axis,upper):
    values=mesh.vertices[:,axis];limit=np.quantile(values,.95 if upper else .05)
    return mesh.vertices[values>=limit if upper else values<=limit].mean(0)

def main():
    began=time.monotonic();assert digest_file(TARGET)==EXPECTED
    source=json.loads((WORK/'source-audit.json').read_text());rows=source['rows'];out=WORK/'first-similarity';assert not out.exists();out.mkdir()
    meshes={r['id']:Mesh.load(ROOT/r['path']) for r in rows};a=json.loads(TARGET.read_text());reader=AtlasReader()
    portal=next(p for p in a['parts'] if p['id']=='TA98HEPATIC:male:region:A12.3.12.002')
    portal_mesh=reader.read_part(a,portal,frame='model:local-ta98:assembled-metres',catalogue_path=TARGET,catalogue_sha=EXPECTED)
    portal_mesh.save(out/'fixed-right-portal-region.npz')
    donor_artery,donor_vein=rows[:2]
    source_gb=meshes['HRA:VH_F_gallbladder'];target_gb=meshes['TA98PART:ZA:Gallbladder:GALLWALL:covering:A05.8.02.004']
    artery_target=meshes['TA98HEPATIC:male:region:A12.2.12.030']
    x=np.array([donor_artery['boundaryProfiles'][0]['center_m'],donor_vein['boundaryProfiles'][0]['center_m'],extremity(source_gb,1,True),extremity(source_gb,1,False),centroid(source_gb)])
    y=np.array([centroid(artery_target),centroid(portal_mesh),extremity(meshes['TA98PART:ZA:Gallbladder:GALLWALL:covering:A05.8.02.005'],1,True),extremity(meshes['TA98PART:ZA:Gallbladder:GALLWALL:covering:A05.8.02.002'],1,False),centroid(target_gb)])
    matrix,metrics=landmark_transform(x,y,weights=[4,4,2,4,2],similarity=True,scale_bounds=(.6,1.6))
    transformed=[]
    for r in rows[:2]:
        original=meshes[r['id']];candidate=transform_mesh(original,matrix,'model:local-ta98:assembled-metres');candidate.save(out/('artery.npz' if 'artery' in r['id'] else 'vein.npz'))
        assert np.array_equal(candidate.faces,original.faces)
        points=np.array(r['boundaryProfiles'][0]['coordinates_m']);mapped=points@matrix[:3,:3].T+matrix[:3,3]
        parent=artery_target if 'artery' in r['id'] else portal_mesh
        near,_=BVH(parent).closest_points(mapped)
        transformed.append(dict(id=r['id'],geometrySha256=candidate.geometry_sha,sourceGeometrySha256=original.geometry_sha,sourceFacesLiteral=True,transformedMouthProfile_m=mapped.tolist(),parentSurfaceGapRange_m=[float(near['distance_m'].min()),float(near['distance_m'].max())],parentSurfaceGapMean_m=float(near['distance_m'].mean())))
    mapped_gb=transform_mesh(source_gb,matrix,'model:local-ta98:assembled-metres');mapped_gb.save(out/'registered-source-gallbladder-context.npz')
    result=dict(schema='cystic-blood-first-similarity/v1',targetCatalogueSha256=EXPECTED,sourceAuditSha256=digest_file(WORK/'source-audit.json'),sourceProximalProofSha256=digest_file(WORK/'source-proximal-attachments.json'),matrix=matrix.tolist(),metrics=metrics,landmarks=dict(names=['cystic artery right-hepatic junction','cystic vein right-portal junction','gallbladder superior neck extremity','gallbladder inferior fundus extremity','gallbladder body area center'],source_m=x.tolist(),target_m=y.tolist(),residual_m=(x@matrix[:3,:3].T+matrix[:3,3]-y).tolist(),correspondenceStatus='Named anatomical locations; sampled centroids/extremities are explicit provisional landmarks requiring native render review'),candidateParts=transformed,sourceParentProof=dict(artery='Every original proximal-ring point lies exactly on original right hepatic artery surface',vein='Original proximal-ring right-portal-branch distance 0 to 0.0288569 mm; source simplification tolerance retained'),closure=['HRA:VH_F_cystic_artery','HRA:VH_F_cystic_vein'],allTargetGeometryFixed=True,skeletonFixed=True,sourceMouthProfilesPreservedBySinglePositiveSimilarity=True,trueMaleParentMouthsEstablished=False,finiteVesselWallsEstablished=False,biliaryLumenExclusionProven=False,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(out/'fit.json',result)
    print(json.dumps(dict(output=str(out),metrics=metrics,candidateParts=[{k:r[k] for k in ['id','parentSurfaceGapRange_m','parentSurfaceGapMean_m']} for r in transformed],eligibleForInstallation=False,elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
