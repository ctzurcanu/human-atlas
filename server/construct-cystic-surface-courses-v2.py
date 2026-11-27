"""One bounded, separately routed cystic artery/vein construction experiment.

Original source sections determine distal calibre. Actual fixed gallbladder
surfaces constrain the new courses; independent CSG tissue/fluid operations
make openings in fixed hepatic parent walls. This never edits a catalogue.
"""
from pathlib import Path
import sys,json,time
from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra,connected_components
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.audit import topology,audit
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.communicating_vascular_network import tube_classifier
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.review import prepare

WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
OUT=WORK/'surface-courses-v4'
FRAME='model:local-ta98:assembled-metres'

def source_sections(mesh,profile):
    v=mesh.vertices;f=mesh.faces;edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0)
    lengths=np.linalg.norm(v[edges[:,1]]-v[edges[:,0]],axis=1)
    graph=coo_matrix((np.r_[lengths,lengths],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(v),len(v))).tocsr()
    lookup={tuple(p):i for i,p in enumerate(v)};root=np.array([lookup[tuple(p)]for p in profile]);dist=dijkstra(graph,directed=False,indices=root,min_only=True)
    levels=np.linspace(float(dist.max())*.012,float(dist.max())*.94,46);centers=[profile.mean(0)];areas=[projected_area(profile)];counts=[]
    diameter=lambda p:float(np.linalg.norm(p[:,None]-p[None,:],axis=2).max())
    diameters=[diameter(profile)];profiles=[profile.tolist()]
    for level in levels:
        cross=edges[(dist[edges[:,0]]<level)!=(dist[edges[:,1]]<level)]
        points={tuple(map(int,e)):v[e[0]]+(v[e[1]]-v[e[0]])*(level-dist[e[0]])/(dist[e[1]]-dist[e[0]])for e in cross}
        adjacency=defaultdict(set)
        for tri in f:
            keys=[tuple(sorted((int(a),int(b))))for a,b in zip(tri,np.roll(tri,-1))if (dist[a]<level)!=(dist[b]<level)]
            if len(keys)==2:adjacency[keys[0]].add(keys[1]);adjacency[keys[1]].add(keys[0])
        assert adjacency and all(len(x)==2 for x in adjacency.values()),'Source geodesic section is not a cycle'
        first=min(adjacency);order=[first];previous=None
        while True:
            nxt=min(adjacency[order[-1]]-({previous}if previous is not None else set()))
            if nxt==first:break
            assert nxt not in order,'Self repeating source section';previous=order[-1];order.append(nxt)
        counts.append(len(order));assert len(order)==len(adjacency),'Source has multiple section cycles: preserve real branches instead of inventing one tube'
        ring=np.array([points[k]for k in order]);centers.append(ring.mean(0));areas.append(projected_area(ring));diameters.append(diameter(ring));profiles.append(ring.tolist())
    return np.array(centers),np.sqrt(np.array(areas)/np.pi),dict(method='literal source edge-geodesic level intersections',singleCycleAtEveryLevel=True,sectionVertexCounts=counts,sectionAreas_m2=areas,sectionDiameters_m=diameters,sectionProfiles_m=profiles,proximalProfileControls=len(root),sampledFractionOfSourceGeodesicLength=[.012,.94],terminalCapExcludedRatherThanForcedIntoStationRecovery=True,maximumGeodesicLength_m=float(dist.max()))

def normal_field(mesh):
    tri=mesh.triangles;normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);values=np.zeros_like(mesh.vertices)
    for k in range(3):np.add.at(values,mesh.faces[:,k],normal)
    length=np.linalg.norm(values,axis=1)
    return np.divide(values,length[:,None],out=np.zeros_like(values),where=length[:,None]>0)

def rigid_profile(profile,center,direction,scale):
    old=profile-profile.mean(0);_,_,basis=np.linalg.svd(old);normal=basis[-1]
    if normal@direction<0:normal=-normal
    cross=np.cross(normal,direction);s=np.linalg.norm(cross);c=normal@direction
    if s<1e-12:rotation=np.eye(3)
    else:
        cross/=s;K=np.array([[0,-cross[2],cross[1]],[cross[2],0,-cross[0]],[-cross[1],cross[0],0]])
        rotation=np.eye(3)+K*s+K@K*(1-c)
    return old@rotation.T*scale+center

def solid(v,f,origin):
    import manifold3d as m
    # Normal-chart seam duplicates share literal coordinates; exact welding
    # changes no anatomical point and avoids treating chart seams as tears.
    used=np.unique(f);points,inverse=np.unique(v[used],axis=0,return_inverse=True)
    mapping=np.full(len(v),-1,np.int64);mapping[used]=inverse
    result=m.Manifold(m.Mesh64((points-origin)*1000,mapping[f].astype(np.uint64)))
    assert result.status()==m.Error.NoError,str(result.status());return result

def arrays(solid,origin):
    m=solid.to_mesh64();return np.asarray(m.vert_properties)[:,:3]/1000+origin,np.asarray(m.tri_verts,dtype=np.int64)

def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir()
    current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(current)=='9cf8b30456620f7b7bd6f52311ba818fc94f935beee24e6e07cb953d14106f64'
    info=json.loads((WORK/'source-audit.json').read_text());fit=json.loads((WORK/'first-similarity/fit.json').read_text());matrix=np.array(fit['matrix']);scale=fit['metrics']['scale']
    exterior=Mesh.load(WORK/'first-similarity/fixed-covering-free-surface.npz');inner=Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz');bvh=BVH(exterior);normals=normal_field(exterior);inner_bvh=BVH(inner)
    results=[];review_meshes=[dict(mesh=str(WORK/'first-similarity/fixed-covering-free-surface.npz'),name='Fixed male gallbladder exterior',color=[.48,.55,.2,.6])]
    for number,role in enumerate(['artery','vein']):
        original=Mesh.load(WORK/('donor-'+str(number)+'.npz'));profile=np.array(info['rows'][number]['boundaryProfiles'][0]['coordinates_m'])
        saved=WORK/'surface-courses-v2'/(role+'-source-sections.json')
        if saved.exists():
            prior=json.loads(saved.read_text());assert prior['originalGeometrySha256']==original.geometry_sha
            source_centers=np.array(prior['centers_m']);radii=np.array(prior['radii_m']);section_evidence=dict(prior['evidence'],resumedFromSavedSha256=digest_file(saved))
        else:source_centers,radii,section_evidence=source_sections(original,profile)
        atomic_json(OUT/(role+'-source-sections.json'),dict(evidence=section_evidence,centers_m=source_centers.tolist(),radii_m=radii.tolist(),originalGeometrySha256=original.geometry_sha))
        path=source_centers@matrix[:3,:3].T+matrix[:3,3];radii*=scale
        host=Mesh.load(WORK/('fixed-target-20.npz'if role=='artery'else'fixed-target-21.npz'))
        host_region=Mesh.load(WORK/'fixed-target-22.npz'if role=='artery'else WORK/'first-similarity/fixed-right-portal-region.npz')
        # The named host region contains both wall faces. Its area centroid is
        # an explicit initial blood-space proposal, checked by the CSG result.
        tri=host_region.triangles;weight=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
        root=np.average(tri.mean(1),axis=0,weights=weight)
        nearest,_=bvh.closest_points(path);distance=nearest['distance_m'];recipient=np.flatnonzero(distance<np.maximum(radii*2,.002))
        assert len(recipient)>5,'Source has no supported recipient section interval'
        start=int(recipient[0]);path[:start]+=np.linspace(1,0,start)[:,None]*(root-path[0]);path[0]=root
        take=np.arange(start,len(path));near,_=bvh.closest_points(path[take]);fn=normals[exterior.faces[near['face_ids']]];n=np.einsum('ij,ijk->ik',near['barycentric'],fn);n/=np.linalg.norm(n,axis=1)[:,None]
        profile_max=np.linalg.norm(profile-profile.mean(0),axis=1).max()/np.sqrt(projected_area(profile)/np.pi)
        offsets=radii[take]*profile_max+.00016
        # Only the final small-calibre recipient ends may enter the outer
        # wall. Their finite radius is retained; no whole-vessel flattening.
        blend=np.clip((np.linspace(0,1,len(take))-.8)/.2,0,1);blend=blend*blend*(3-2*blend)
        offsets-=blend*(radii[take]*profile_max*.4+.00016)
        path[take]=near['points']+offsets[:,None]*n
        path[1:-1]=gaussian_filter1d(path,1.,axis=0,mode='nearest')[1:-1];path[0]=root
        # Smooth local curves use existing section-preserving sweep frames.
        s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))];assert np.all(np.diff(s)>1e-8)
        parameter=np.linspace(0,s[-1],64);curve=CubicSpline(s,path,bc_type='natural')(parameter);radius=np.interp(parameter,s,radii)
        direction=curve[1]-curve[0];direction/=np.linalg.norm(direction)
        first=rigid_profile(profile,curve[0],direction,scale)
        outer_v,outer_f,outer_caps=tube_classifier(curve,radius,ring_size=18,start_ring=first,extend_joints=False)
        ratio=.82 # Explicit illustrative wall construction, not measured HRA histology.
        inner_first=curve[0]+ratio*(first-curve[0]);fluid_curve=curve.copy()
        fluid_curve[0]-=direction*radius[0]*1.6
        distal=curve[-1]-curve[-2];distal/=np.linalg.norm(distal);fluid_curve[-1]+=distal*radius[-1]*1.6
        fluid_first=inner_first-direction*radius[0]*1.6
        inner_v,inner_f,inner_caps=tube_classifier(fluid_curve,radius*ratio,ring_size=18,start_ring=fluid_first,extend_joints=False)
        atomic_npz(OUT/(role+'-section-fields.npz'),centers=curve,radii=radius,proximal_source_profile=profile,proximal_candidate_profile=first,outer_vertices=outer_v,outer_faces=outer_f,fluid_vertices=inner_v,fluid_faces=inner_f,outer_classifier_caps=np.concatenate(outer_caps),fluid_classifier_caps=np.concatenate(inner_caps))
        origin=root.copy();outer_solid=solid(outer_v,outer_f,origin);fluid_solid=solid(inner_v,inner_f,origin)
        native_solid=solid(host.vertices,host.faces,origin)
        joined=(native_solid+outer_solid)-fluid_solid
        joined_v,joined_f=arrays(joined,origin);network=Mesh(joined_v,joined_f,FRAME,source=dict(kind='diagnostic actual fixed-host CSG outer union minus blood-lumen cutter',hostGeometrySha256=host.geometry_sha,sourceDonor=original.geometry_sha))
        network.save(OUT/(role+'-parent-with-cystic-wall.npz'))
        # Save the original source-profile coordinates and new section fields
        # independently of boolean re-triangulation for exact replay.
        # Vessel wall only; classifier disks never become anatomy.
        branch_solid=outer_solid-fluid_solid;bv,bf=arrays(branch_solid,origin);branch=Mesh(bv,bf,FRAME);branch.save(OUT/(role+'-branch-wall.npz'))
        contact,proof=pair_contacts(branch,inner,0.,max_pairs=250000,max_tests=2000000)
        atomic_npz(OUT/(role+'-innermost-witness.npz'),**contact);atomic_json(OUT/(role+'-innermost-contacts.json'),proof)
        top=topology(network);branch_top=topology(branch)
        surface=Mesh(outer_v,outer_f[np.setdiff1d(np.arange(len(outer_f)),np.concatenate(outer_caps))],FRAME)
        surface.save(OUT/(role+'-outer-surface.npz'))
        minimum=inner_bvh.closest_points(np.r_[curve,outer_v],backend='bvh')[0]['distance_m'].min()
        result=dict(role=role,sourceSections=section_evidence,scale=scale,calibreRelativeTolerance=.15,calibreRelativeError=0.,proximal18ControlShapePreservedByPositiveSimilarity=True,sourceProfileSha256=digest_file(OUT/(role+'-section-fields.npz')),recipientStationRange=[start,len(path)-1],hostFixedSourceGeometrySha256=host.geometry_sha,networkGeometrySha256=network.geometry_sha,actualHostWallCSGPerformed=True,joinedTopology=top,branchTopology=branch_top,innermostContacts=proof,minimumSampledInnerSurfaceDistance_m=float(minimum),finiteWallRatio=ratio,wallHistologyProvisional=True,eligibleForInstallation=False)
        atomic_json(OUT/(role+'-result.json'),result);results.append(result)
        review_meshes.append(dict(mesh=str(OUT/(role+'-parent-with-cystic-wall.npz')),name='Actual fixed-host '+role+' with separately rounded cystic branch',color=[.82,.16,.12,1]if role=='artery'else[.15,.25,.78,1]))
        # Stop this finite experiment on a new lumen crossing. Do not tune or
        # retry within the stage and do not start a second global fit.
        if proof['crossing_face_pairs']:
            atomic_json(OUT/'stopped.json',dict(reason='New branch crosses actual innermost gallbladder wall',role=role,crossingFacePairs=proof['crossing_face_pairs'],completedRoles=[r['role']for r in results],elapsedSeconds=time.monotonic()-began));break
    spec=dict(frame=FRAME,meshes=review_meshes,views=[dict(name='surface',direction=[0,0,1]),dict(name='oblique',direction=[1,.2,1]),dict(name='section',direction=[0,0,1],section=dict(normal=[0,0,1],point=exterior.bounds.mean(0).tolist(),clear_outer=True))],resolution=640,wall_seconds=45)
    atomic_json(OUT/'review-spec.json',spec);queued=prepare(None,spec);atomic_json(OUT/'review-queued.json',queued)
    atomic_json(OUT/'construction.json',dict(schema='separate-cystic-surface-courses/v2',roles=results,review=queued,eligibleForInstallation=False,installed=False,allOrgansAndSkeletonFixed=True,elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(output=str(OUT),roles=[dict(role=r['role'],innermostCrossings=r['innermostContacts']['crossing_face_pairs'],joinedComponents=len(r['joinedTopology']['components']),joinedClosed=r['joinedTopology']['closed_oriented_manifold'])for r in results],reviewRequest=queued['request_id'],elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
