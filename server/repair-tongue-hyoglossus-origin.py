"""One source-footprint HG origin repair; fixed hyoid, protected insertion.

Uses native superior body/greater-horn facets and a finite broad sheet merging
into the actual HG section. No external mesh is copied and the old failed
boundary loft is not reused. Bony adjacency is explicit, submicrometre and
unjoined; it does not establish histological enthesis.
"""
from __future__ import annotations
import argparse,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.tissue import layer_arrangement
from geometry_tools.audit import audit,topology
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.distance import point_triangle

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(path));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
patch=module('mucosal_patch','propose-tongue-mucosal-patch.py')
fit=module('mucosal_fit','fit-tongue-mucosal-insertions.py')
encoding=module('origin_encoding','encode-tongue-genioglossus-attachments.py')
WORK=patch.WORK

def interval_at(section,z):
    hits=[]
    for a,b in section:
        d=a[2]-z;e=b[2]-z
        if d==0:hits.append(a[0])
        if d*e<0:hits.append(a[0]+d/(d-e)*(b[0]-a[0]))
    if len(hits)<2:raise ValueError('Actual native HG section has no finite interval at this footprint station')
    return float(min(hits)),float(max(hits))

def broad_sheet(footprint,muscle,stationed_rim=False):
    used,ff=np.unique(footprint.faces,return_inverse=True);ff=ff.reshape(-1,3);base=footprint.vertices[used]
    y_top=1.619;section=patch.segments(muscle,1,y_top)
    z0,z1=base[:,2].min(),base[:,2].max();tz0=section[:,:,2].min()+.002;tz1=section[:,:,2].max()-.002
    if tz1<=tz0 or base[:,1].max()>=y_top-.001:raise ValueError('Actual HG/hyoid extents leave no finite protected origin transition')
    # The native footprint is a height graph: all actual source facet normals
    # point superiorly. Map its ordered Z/X chart into an interior actual HG
    # section, using a positive transverse scale and monotone longitudinal map.
    tri=base[ff];norm=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm/=np.linalg.norm(norm,axis=1)[:,None]
    if norm[:,1].min()<=.4:raise ValueError('Footprint is not the witnessed superior hyoid graph')
    knots=np.linspace(z0+1e-8,z1-1e-8,49);centers=[]
    compact=Mesh(base,ff,muscle.frame)
    for z in knots:
        line=patch.segments(compact,2,z).reshape(-1,3)
        if not len(line):raise ValueError('Disconnected hyoid footprint chart')
        centers.append((line[:,0].min()+line[:,0].max())/2)
    source_trend=np.polyfit(knots,centers,1)
    target=np.empty_like(base);target[:,1]=y_top;target[:,2]=tz0+(base[:,2]-z0)/(z1-z0)*(tz1-tz0)
    target_knots=np.linspace(tz0,tz1,49);target_centers=[sum(interval_at(section,z))/2 for z in target_knots]
    target_trend=np.polyfit(target_knots,target_centers,1)
    # A strictly positive affine Z/X chart avoids reversing a coarse footprint
    # triangle while fitting a curved, thin native cross-section. Its actual
    # points, triangle centres and edge midpoints still must lie inside HG.
    target[:,0]=np.polyval(target_trend,target[:,2])+.08*(base[:,0]-np.polyval(source_trend,base[:,2]))
    samples=np.r_[target,target[ff].mean(1),(target[ff[:,0]]+target[ff[:,1]])/2,(target[ff[:,1]]+target[ff[:,2]])/2]
    winding=np.abs(igl.winding_number(muscle.vertices,muscle.faces,samples))
    if winding.min()<.99:raise ValueError('Broad-sheet top is not entirely inside the actual native HG course')
    # A small overlap into native HG is intentional construction data; the
    # finite union removes every internal overlap face from the result.
    intermediate=base.copy();intermediate[:,1]=1.617
    if base[:,1].max()>=1.616:raise ValueError('Native footprint has no room for the protected intermediate level')
    if stationed_rim:
        # Same positive affine narrowing, with a literal common station at the
        # reviewed source cut. Only rim tessellation changes; physical endpoints
        # and all anatomical dimensions remain unchanged.
        stations=sorted(set(np.r_[np.linspace(1.,0.,14),.9]),reverse=True)
        levels=[]
        for t in stations:
            q=intermediate+t*(target-intermediate)
            if t==.9:q[:,1]=1.6188
            levels.append(q)
        fields=np.array(levels+[base]);names=[f'slab-{i}'for i in range(len(fields)-1)]
        if len(fields)>16 or np.max(np.abs(np.diff([q[0,1]for q in levels])))>.0005:raise ValueError('HG station count/height exceeds the explicit cap')
        projected=base.copy();projected[:,1]=0.
        _,projection=pair_contacts(Mesh(projected,ff,muscle.frame),Mesh(projected,ff,muscle.frame),0.,selfcheck=True,max_pairs=100000,max_tests=1000000)
        if projection['crossing_face_pairs']or projection['coplanar_overlap_face_pairs']or projection['truncated']:raise ValueError('Native footprint X/Z projection overlaps or is not completely tested')
        z_scale=(tz1-tz0)/(z1-z0)
        determinants=[float((1-t+t*.08)*(1-t+t*z_scale))for t in stations]
        if min(determinants)<=0:raise ValueError('Global affine narrowing interpolant reverses')
    else:
        stations=None;projection=None;determinants=None
        fields=np.array([target,intermediate,base]);names=['upper-fan','basal-sheet']
    v,surfaces,cells,prisms=layer_arrangement(fields,ff,names)
    faces=np.vstack([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]
        for uses in cells.values()for u in uses if u['neighbor']is None])
    if not prisms['finite_wall_cells_pass']:raise ValueError('Source-footprint muscular sheet has reversed/collapsed finite prisms')
    return Mesh(v,faces,muscle.frame),dict(sourceFootprintGeometrySha256=footprint.geometry_sha,
        sourceFootprintFaces=len(ff),nativeSuperiorFaceNormalYMinimum=float(norm[:,1].min()),
        footprintAPExtentM=float(z1-z0),targetActualSectionY=y_top,allTargetSamplesInsideNativeHG=True,
        targetChart='positive affine X/Z chart; transverse derivative .08',intermediateFlatLevelY=1.617,
        stationedRim=stationed_rim,stationFractions=stations,fieldCount=len(fields),
        globalAffineStationDeterminants=determinants,completeProjectedFootprintCheck=projection,
        finitePrisms=prisms,notThinWire=True,reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
        referenceScope='HG thin planar origin spans body and greater cornu, then fans into the tongue')

def encode_bone_adjacency(candidate,bone,tree):
    points,alias=np.unique(candidate.vertices,axis=0,return_inverse=True);near,_=tree.closest_points(points,backend='open3d')
    controls=np.flatnonzero(near['distance_m']<1e-6);proposed=points.copy();tri=bone.triangles
    area=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal=area/np.linalg.norm(area,axis=1)[:,None];moves=[]
    for control in controls:
        if near['distance_m'][control]>=.45e-6:continue
        p=points[control];ids=tree.query_bounds(p,p,distance_m=1e-9);d,_,_=point_triangle(np.repeat(p[None],len(ids),axis=0),tri[ids]);ids=ids[d<1e-18]
        if not len(ids):raise ValueError('New HG contact lacks actual cortical facet witnesses')
        signed=np.einsum('ij,ij->i',p-tri[ids,0],normal[ids]);delta=encoding.minimum_offset(normal[ids],signed);proposed[control]+=delta
        moves.append(dict(control=int(control),actualCorticalFaces=ids.tolist(),displacementM=float(np.linalg.norm(delta))))
    origin=np.array([0.,1.625,0.]);stored=(proposed-origin).astype('<f4');actual=stored.astype(float)+origin
    result=Mesh(actual[alias],candidate.faces,candidate.frame,candidate.attributes,
        dict(derived_from_geometry=candidate.geometry_sha,geometricAttachment='adjacent-unjoined',nominalClearanceM=.5e-6,maximumToleranceM=1e-6))
    contacts,_=tree.closest_points(actual[controls],backend='open3d')
    return result,stored[alias],dict(contactControls=len(controls),nominalTissueSideClearanceM=.5e-6,
        actualContactGapRangeM=[float(contacts['distance_m'].min()),float(contacts['distance_m'].max())],moves=moves,
        literalSharedCoordinates=False,histologicalEnthesis=False)

def adjacent_native_footprint(footprint,bone,tree):
    """Move new tissue controls outside all incident literal cortical facets."""
    v=footprint.vertices.copy();tri=bone.triangles;area=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal=area/np.linalg.norm(area,axis=1)[:,None]
    moves=[]
    for control in np.unique(footprint.faces):
        p=v[control];ids=tree.query_bounds(p,p,distance_m=1e-9);d,_,_=point_triangle(np.repeat(p[None],len(ids),axis=0),tri[ids]);ids=ids[d<1e-18]
        if not len(ids):raise ValueError('Native origin footprint control lacks literal cortical witnesses')
        signed=np.einsum('ij,ij->i',p-tri[ids,0],normal[ids]);delta=encoding.minimum_offset(normal[ids],signed);v[control]+=delta
        moves.append(dict(sourceFootprintControl=int(control),actualCorticalFaces=ids.tolist(),displacementM=float(np.linalg.norm(delta))))
    return Mesh(v,footprint.faces,footprint.frame),dict(sourceNativeFootprintGeometrySha256=footprint.geometry_sha,
        newTissueControlsOnly=True,nominalTissueSideClearanceM=.5e-6,maximumAllowedGapM=1e-6,moves=moves)

def main(output,replace_lower_origin=False,repair_rim_topology=False,stationed_rim=False):
    began=time.monotonic();output=output.resolve()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new local tongue stage')
    bones=[Mesh.load(WORK/f'source-v1/meshes/native-{i:02}.npz')for i in [13,14,15,16,17]]
    v=[];f=[];offset=0
    for bone in bones:v.append(bone.vertices);f.append(bone.faces+offset);offset+=len(bone.vertices)
    hyoid=Mesh(np.vstack(v),np.vstack(f),bones[0].frame);hyoid_solid,_=make_solid(hyoid.vertices,hyoid.faces);tree=BVH(hyoid)
    jaw=Mesh.load(WORK/'source-v1/meshes/native-05.npz');mucosa=Mesh.load(WORK/'exposed-mucosa-v1/packed/TA98TONGUE_male_mucosa.npz')
    temporal=[Mesh.load(WORK/f'source-v1/meshes/native-{i:02}.npz')for i in [6,7]]
    output.mkdir(parents=True);rows=[]
    for side in ['left','right']:
        if time.monotonic()-began>48:raise TimeoutError('HG worker stopped before 60-second stage deadline')
        source=WORK/f'mucosal-insertions-v2-certified-rim/HG-{side}.npz';muscle=Mesh.load(source)
        footprint_path=WORK/f'native-attachment-candidate-v1/hyoglossus-{side}-bone-origin-surface.npz';footprint=Mesh.load(footprint_path)
        adjacency_input=None
        if replace_lower_origin:footprint,adjacency_input=adjacent_native_footprint(footprint,hyoid,tree)
        if repair_rim_topology and stationed_rim:raise ValueError('Choose one reviewed producer correction, never competing variants')
        sheet,plan=broad_sheet(footprint,muscle,stationed_rim)
        if stationed_rim:
            standalone=audit(sheet,require_closed=True,max_pairs=100000,max_tests=1000000)
            local=(sheet.vertices-[0.,1.625,0.]).astype('<f4').astype(float)+[0.,1.625,0.]
            standalone_stored=audit(Mesh(local,sheet.faces,sheet.frame),require_closed=True,max_pairs=100000,max_tests=1000000)
            preflight=dict(raw=standalone,stored=standalone_stored,coordinatesFollowSameAffineField=True)
            atomic_json(output/f'{side}-standalone-preflight.json',preflight)
            components=standalone['topology']['components']
            if not standalone['passed']or not standalone_stored['passed']or len(components)!=1 or components[0]['closed_orientable_genus']!=0:
                raise ValueError('Stationed standalone sheet fails global self/component/genus gate before union')
            plan['standalonePreflight']=preflight
        if repair_rim_topology:
            rim_module=module('hg_rim_topology','tongue_hg_rim_triangulation.py')
            sheet,rim_evidence=rim_module.repair(sheet,side,seconds=18.)
            atomic_json(output/f'{side}-standalone-rim-enumeration.json',rim_evidence)
            if not rim_evidence['passed']:raise ValueError('Standalone sheet failed the authorized fixed-boundary rim repair before union')
            plan['fixedBoundaryRimRepair']=rim_evidence
        sheet.save(output/f'{side}-original-broad-origin-sheet.npz')
        if adjacency_input:plan['nativeCorticalAdjacency']=adjacency_input
        native,oid=make_solid(muscle.vertices,muscle.faces,properties=muscle.attributes['normals'][muscle.faces]);addition,_=make_solid(sheet.vertices,sheet.faces)
        if replace_lower_origin:
            retained_upper=native.trim_by_plane([0.,1.,0.],1.6188)
            repaired=retained_upper+addition
            plan.update(lowerNativeOriginReplaced=True,nativeCutPlaneY=1.6188,shortConnectedOverlapM=.0002,
                        originalLowerMuscleNotSuperimposed=True)
        else:repaired=(native+addition)-hyoid_solid
        before,retained,_=fit.decoded_candidate(repaired,muscle,oid,muscle.frame)
        # Use the Float64 Boolean positions until contact controls have moved;
        # do not quantize their shared cortical seam prematurely.
        raw=repaired.to_mesh64();tri=np.asarray(raw.vert_properties,float)[np.asarray(raw.tri_verts,int),:3]
        before=Mesh(tri.reshape(-1,3),before.faces,before.frame,before.attributes,before.source)
        candidate,stored,adjacency=encode_bone_adjacency(before,hyoid,tree)
        dest=output/f'HG-{side}.npz';candidate.save(dest);atomic_npz(output/f'HG-{side}-stored.npz',positions=stored,faces=candidate.faces,
            geometryOrigin=np.array([0.,1.625,0.]),frame=np.array(candidate.frame),attribute_normals=candidate.attributes['normals'])
        checks=audit(candidate,require_closed=True,max_pairs=100000,max_tests=1000000);contacts={}
        for name,other in [('hyoid',hyoid),('mandible',jaw),('temporal-left',temporal[0]),('temporal-right',temporal[1]),('mucosa',mucosa)]:
            arrays,summary=pair_contacts(candidate,other,0.,max_pairs=100000,max_tests=1000000);contacts[name]=summary;atomic_npz(output/f'{side}-{name}-contact.npz',**arrays)
        keys={fit.oriented_triangle(t)for t in candidate.triangles};protected=np.all(muscle.triangles[:,:,1]>=1.621,axis=1);preserved=sum(fit.oriented_triangle(t)in keys for t in muscle.triangles[protected])
        topo=checks['topology'];topology_pass=len(topo['components'])==1 and topo['components'][0]['closed_orientable_genus']==0
        passed=checks['passed']and topology_pass and all(x['crossing_face_pairs']==0 and x['coplanar_overlap_face_pairs']==0 and not x['truncated']for x in contacts.values())and preserved==int(protected.sum())and adjacency['actualContactGapRangeM'][1]<=1e-6
        row=dict(side=side,candidate=str(dest),candidateSha256=digest_file(dest),source=str(source),sourceSha256=digest_file(source),
            nativeFootprint=str(footprint_path),nativeFootprintSha256=digest_file(footprint_path),construction=plan,adjacency=adjacency,
            checks=checks,oneComponentGenusZero=topology_pass,contacts=contacts,protectedSuperiorFaces=int(protected.sum()),actualStoredSuperiorFacesLiteral=preserved,
            originalVolumeMl=float(native.volume()*1e6),repairedVolumeMl=float(repaired.volume()*1e6),passed=passed)
        rows.append(row);atomic_json(output/f'{side}-audit.json',row)
    report=dict(schema='hyoglossus-native-broad-footprint-repair/v2',rows=rows,passed=all(r['passed']for r in rows),
        fixedBonesModified=False,predicatesRelaxed=False,anatomicalAcceptance=False,installed=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(output/'candidate.json',report)
    print(dict(passed=report['passed'],elapsed=report['elapsedSeconds'],rows=[dict(side=r['side'],self=r['checks']['passed'],contacts={k:v['crossing_face_pairs']for k,v in r['contacts'].items()},protected=[r['protectedSuperiorFaces'],r['actualStoredSuperiorFacesLiteral']])for r in rows]),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--replace-lower-origin',action='store_true');p.add_argument('--repair-rim-topology',action='store_true');p.add_argument('--stationed-rim',action='store_true');a=p.parse_args();main(a.output,a.replace_lower_origin,a.repair_rim_topology,a.stationed_rim)
