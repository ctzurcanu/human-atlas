"""Reconstruct current source charts for computation, then scatter ownership.

The suppressed old whole carrier is never read. Original master vertex IDs
and explicit face ancestry gather the CURRENT charts, even across different
packed banks. Distinct master IDs are not welded because their UVs/normals
may intentionally differ. Exact-equal positions are welded only by the
existing source-section recovery's incidence checks.
"""
from dataclasses import dataclass
import numpy as np
from ..model import Mesh
from .source_profile_axis import source_profile_graph


@dataclass
class RegionAssembly:
    whole: Mesh
    meshes: dict
    local_to_whole: dict
    face_rows: dict
    proof: dict


def gather_current_regions(records,meshes,*,family_id,expected_owner_ids):
    owners=sorted(meshes)
    if (not owners or set(records)!=set(meshes)or set(owners)!=set(expected_owner_ids)
            or not isinstance(family_id,str)or not family_id):
        raise ValueError('Explicit complete current owners, records and source family required')
    if len({mesh.frame for mesh in meshes.values()})!=1:
        raise ValueError('Current regional charts must share one explicit world frame')
    schemas={(key,values.dtype.str,values.shape[1:])for key,values in meshes[owners[0]].attributes.items()}
    master={};ancestry={};source_ids=set();contracts=set()
    for pid in owners:
        record=records[pid];derived=record.get('sourceDerivedRegion',{})
        if derived.get('sourcePartId')!=family_id or not derived.get('sourceId'):
            raise ValueError('Current chart source-family provenance disagrees')
        source_ids.add(derived['sourceId']);contracts.add(derived.get('boundaryContractSha256'))
        m=meshes[pid];ids=np.asarray(m.source.get('source_vertex_ids',[]));faces=np.asarray(derived.get('sourceFaceIds',[]))
        if (ids.ndim!=1 or len(ids)!=len(m.vertices)or not np.issubdtype(ids.dtype,np.integer)
                or np.any(ids<0)or len(np.unique(ids))!=len(ids)
                or faces.ndim!=1 or len(faces)!=len(m.faces)or not np.issubdtype(faces.dtype,np.integer)
                or np.any(faces<0)or len(np.unique(faces))!=len(faces)):
            raise ValueError('One original master ID per current control and face ancestry per triangle required')
        if {(key,values.dtype.str,values.shape[1:])for key,values in m.attributes.items()}!=schemas:
            raise ValueError('Current charts disagree on attribute schema or word dtype')
        master[pid]=ids;ancestry[pid]=faces
    if len(source_ids)!=1 or len(contracts)!=1:
        raise ValueError('Mixed source identities or partition contracts are not one original family')
    all_ids=np.unique(np.concatenate(list(master.values())))
    face_ids=np.sort(np.concatenate(list(ancestry.values())))
    if len(np.unique(face_ids))!=len(face_ids):
        raise ValueError('A source face has duplicate current tissue ownership')
    template=meshes[owners[0]];positions=np.empty((len(all_ids),3));occupied=np.zeros(len(all_ids),bool)
    attrs={key:np.empty((len(all_ids),)+value.shape[1:],dtype=value.dtype)for key,value in template.attributes.items()}
    faces=np.empty((len(face_ids),3),int);maps={};face_rows={};members=[]
    for pid in owners:
        m=meshes[pid];dest=np.searchsorted(all_ids,master[pid]);seen=occupied[dest]
        if not np.array_equal(positions[dest[seen]],m.vertices[seen]):
            raise ValueError('Shared original master controls have inconsistent current coordinates')
        for key in attrs:
            if not np.array_equal(attrs[key][dest[seen]],m.attributes[key][seen]):
                raise ValueError('Shared original master controls have inconsistent current '+key)
            attrs[key][dest]=m.attributes[key]
        positions[dest]=m.vertices;occupied[dest]=True;maps[pid]=dest
        target=np.searchsorted(face_ids,ancestry[pid]);faces[target]=dest[m.faces];face_rows[pid]=target
        members.append(dict(partId=pid,geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,
                            controls=len(m.vertices),faces=len(m.faces)))
    assert occupied.all()
    source=dict(kind='current-owned-source-regions',familyId=family_id,sourceId=next(iter(source_ids)),
                source_vertex_ids=all_ids.tolist(),source_face_ids=face_ids.tolist(),currentOwners=owners,
                suppressedOldWholeCoordinatesUsed=False)
    whole=Mesh(positions,faces,template.frame,attrs,source)
    proof=dict(complete=True,familyId=family_id,owners=owners,members=members,currentControls=len(all_ids),
        currentFaces=len(face_ids),sourceVertexIds=all_ids.tolist(),sourceFaceIds=face_ids.tolist(),
        allCurrentFacesAndAttributesRetained=True,sourceFaceOwnershipUnique=True,
        ancestryContiguousFromZero=bool(np.array_equal(face_ids,np.arange(len(face_ids)))),
        sharedMasterCoordinatesNormalsAndSignsExact=True,suppressedOldWholeCoordinatesUsed=False,
        computationOnly=True,geometryModified=False,anatomicalAcceptance=False)
    return RegionAssembly(whole,dict(meshes),maps,face_rows,proof)


def recover_current_sections(assembly,*,ring_size=6):
    """Recover actual annuli in original face ancestry order, not chart rims."""
    ids=np.asarray(assembly.proof['sourceFaceIds'],int);lookup={int(x):i for i,x in enumerate(ids)}
    pairs=np.array([(lookup[int(i)],lookup[int(i)+1])for i in ids if i%2==0 and int(i)+1 in lookup],int).reshape(-1,2)
    if not len(pairs):
        raise ValueError('No explicitly consecutive original exporter triangle pairs survive')
    return source_profile_graph(assembly.whole.vertices,assembly.whole.faces,triangle_pairing=pairs,ring_size=ring_size)


def scatter_region_pose(assembly,posed_whole):
    """Scatter already-reviewed whole coordinates/normal words to current owners."""
    whole=assembly.whole
    if posed_whole.frame!=whole.frame or posed_whole.vertices.shape!=whole.vertices.shape:
        raise ValueError('Posed current whole must retain its frame and all master controls')
    if not np.array_equal(posed_whole.faces,whole.faces)or set(posed_whole.attributes)!=set(whole.attributes):
        raise ValueError('Whole pose must retain every current face and attribute field')
    for key in whole.attributes:
        if posed_whole.attributes[key].dtype!=whole.attributes[key].dtype:
            raise ValueError('Whole attribute word dtype changed')
        if key not in ['normals','normalSigns']and not np.array_equal(posed_whole.attributes[key],whole.attributes[key]):
            raise ValueError('Pose changed a non-normal source attribute')
    result={}
    for pid,m in assembly.meshes.items():
        take=assembly.local_to_whole[pid]
        result[pid]=Mesh(posed_whole.vertices[take],m.faces,m.frame,
                        {key:value[take]for key,value in posed_whole.attributes.items()},dict(m.source))
    return result


def scatter_whole_constraints(assembly,raw_ids,target_points):
    """Scatter point constraints; fragments are never called complete profiles."""
    ids=np.asarray(raw_ids);q=np.asarray(target_points,float)
    if (ids.ndim!=1 or not np.issubdtype(ids.dtype,np.integer)or q.shape!=(len(ids),3)
            or not np.isfinite(q).all()or len(np.unique(ids))!=len(ids)
            or np.any(ids<0)or np.any(ids>=len(assembly.whole.vertices))):
        raise ValueError('Unique whole control IDs and corresponding finite targets required')
    lookup={int(i):j for j,i in enumerate(ids)};result={}
    for pid,take in assembly.local_to_whole.items():
        raw=np.flatnonzero(np.isin(take,ids));targets=np.array([q[lookup[int(take[i])]]for i in raw]).reshape(-1,3)
        result[pid]=dict(rawIds=raw.tolist(),targetPoints_m=targets.tolist())
    return result
