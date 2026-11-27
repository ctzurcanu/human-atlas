"""Source-bound selectable surface annotations, not additional liver solids."""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.open3d_backend import library

start=time.monotonic()
base=ROOT/'.local-models/geometry-workbench/liver-completion'
out=base/'liver-surface-controls-v143';out.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
raw=path.read_bytes();atlas=json.loads(raw);reader=AtlasReader();o3d=library()
parts={p['id']:p for p in atlas['parts']};concepts={c['id']:c for c in atlas['concepts']}
host=parts['TA98PART:ZA:Liver:capsule'];mesh=reader.read_part(atlas,host,frame='model:local-ta98:assembled-metres')
outer=next(g for g in host['tissueDrawGroups']if g['start']==0 and not g.get('neighbor'))
ids=np.arange(outer['count']//3);tri=mesh.triangles[ids];centers=tri.mean(1)
normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
# Adjacency uses exact coordinates, including normal-split source corners.
_,inv=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);faces=inv.reshape(-1,3)
edges=np.vstack([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]);owners=np.tile(np.arange(len(faces)),3)
_,ein,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
order=np.argsort(ein,kind='stable');offset=np.r_[0,np.cumsum(count)[:-1]];pairs=offset[count==2]
adj=np.c_[owners[order[pairs]],owners[order[pairs+1]]]
evidence={};patches={}

def largest(mask,seed=None):
    chosen=np.flatnonzero(mask)
    if not len(chosen):raise ValueError('No actual surface facets for the requested feature')
    keep=mask[adj].all(1);pairs=adj[keep]
    graph=coo_matrix((np.ones(len(pairs)*2),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(mask),len(mask))).tocsr()
    _,labels=connected_components(graph,directed=False)
    label=labels[seed]if seed is not None else max(set(labels[chosen]),key=lambda k:sum(labels[chosen]==k))
    return np.flatnonzero(mask&(labels==label))

def target(name,part_ids):
    selected=[parts[p]for p in part_ids if p in parts and not parts[p].get('suppressed')]
    if not selected:raise ValueError('No current adjoining anatomy: '+name)
    vertices=[];ff=[];n=0;bindings={}
    for p in selected:
        m=reader.read_part(atlas,p,frame=mesh.frame);vertices.append(m.vertices);ff.append(m.faces+n);n+=len(m.vertices)
        bindings[p['id']]=dict(geometry=m.geometry_sha,attributes=m.attributes_sha)
    v=np.concatenate(vertices);f=np.concatenate(ff)
    scene=o3d.t.geometry.RaycastingScene(nthreads=1)
    scene.add_triangles(o3d.core.Tensor(v.astype(np.float32)),o3d.core.Tensor(f.astype(np.uint32)))
    closest=scene.compute_closest_points(o3d.core.Tensor(centers.astype(np.float32)),nthreads=1)['points'].numpy().astype(float)
    distance=np.linalg.norm(closest-centers,axis=1)
    evidence[name]=dict(parts=bindings,closest_distance_m=float(distance.min()))
    return distance,closest

def concept_ids(cid):
    return concepts[cid]['elements']

targets={
 'cardiac': ['ZA:Left ventricle','ZA:Right ventricle'],
 'gastric':concept_ids('ZA:Stomach'),
 'oesophageal':concept_ids('ZA:Oesophagus'),
 'gallbladder':concept_ids('ZA:Gallbladder'),
 'duodenal':concept_ids('ZA:Duodenum'),
 'colic':concept_ids('ZA:Transverse colon')+concept_ids('ZA:Ascending colon'),
 'renal':['RENALDETAIL:male:right:capsule'],
 'suprarenal':['TA98PART:male:Suprarenal:right:capsule'],
 'caval':['TA98HEPATIC:male:hepatic_vein'],
}
distances={name:target(name,pids)for name,pids in targets.items()}
directions=np.array([[0,1,0],[0,0,1],[-1,0,0],[0,0,-1],[0,-1,0]],float)
side=(normals@directions.T).argmax(1)
for number,k in [(3,0),(5,1),(6,2),(7,3)]:patches[number]=largest(side==k)
patches[2]=largest(side!=4);patches[12]=largest(side==4)

def footprint(number,name,allowed=None,width=.004):
    distance,closest=distances[name]
    valid=np.ones(len(centers),bool)if allowed is None else allowed
    # The opposite wall of a concavity must not win by a short back-face ray.
    facing=np.einsum('ij,ij->i',normals,closest-centers)>=-1e-7
    valid &= facing
    seed=int(np.argmin(np.where(valid,distance,np.inf)))
    threshold=float(distance[seed]+width)
    mask=valid&(distance<=threshold)
    patches[number]=largest(mask,seed)
    evidence[name].update(region_facets=len(patches[number]),seed_facet=seed,extent_rule='Connected nearest-approach footprint, provisional anatomical boundary',distance_band_m=width)

footprint(4,'cardiac',side==0,.006)
footprint(13,'gallbladder',side!=0,.003)
for number,name in [(18,'oesophageal'),(19,'gastric'),(20,'duodenal'),(21,'colic'),(22,'renal'),(23,'suprarenal')]:
    footprint(number,name,side!=0,.004)
footprint(9,'caval',side!=1,.004)
# Broad bare-area and the two fissure charts are provisional surface extents
# anchored to the current gross lobe interfaces, not invented solid ligaments.
patches[8]=largest((side==3)&(centers[:,0]<-.015)&(centers[:,1]>1.265))
left=reader.read_part(atlas,parts['TA98PART:ZA:Liver:left'],frame=mesh.frame)
caudate=reader.read_part(atlas,parts['TA98PART:ZA:Liver:caudate'],frame=mesh.frame)
quadrate=reader.read_part(atlas,parts['TA98PART:ZA:Liver:quadrate'],frame=mesh.frame)
from scipy.spatial import cKDTree
dl=cKDTree(left.triangles.mean(1)).query(centers)[0]
dc=cKDTree(caudate.triangles.mean(1)).query(centers)[0]
dq=cKDTree(quadrate.triangles.mean(1)).query(centers)[0]
for number,other,allowed in [(10,dc,side!=1),(14,dq,side!=0)]:
    score=np.maximum(dl,other);seed=int(np.argmin(np.where(allowed,score,np.inf)))
    patches[number]=largest(allowed&(score<=score[seed]+.003),seed)
    evidence[str(number)]=dict(source='Actual left/caudate or left/quadrate interface approach',minimum_max_distance_m=float(score[seed]),boundary_provisional=True)

inventory=json.loads((base/'male-liver-section-install-v142/inventory.json').read_text())
authority={t['authority']['id']:t['authority']for t in inventory['terms']}
groups=[]
for number,faceids in sorted(patches.items()):
    term=f'A05.8.01.{number:03d}';row=authority[term]
    center=centers[faceids].mean(0);face=int(faceids[np.argmin(np.linalg.norm(centers[faceids]-center,axis=1))]);point=centers[face]
    cid='TA98LIVER:male:surface:'+term
    group=dict(term=term,conceptId=cid,partIds=[host['id']],representation='surface-region',
        surface=dict(partId=host['id'],faces=faceids.tolist(),sourceBankSha256=atlas['liverIntegration']['bankSha256']),
        anchor=dict(partId=host['id'],face=face,barycentric=[1/3]*3,point=point.tolist()),
        boundaryProvisional=True)
    groups.append(group)
    parent='TA98LIVER:male:surface:'+row['parent_id']if int(row['parent_id'].split('.')[-1])in patches and row['parent_id'].startswith('A05.8.01.')else'ZA:Liver'
    atlas['concepts'].append(dict(id=cid,name=row['name_en'].capitalize(),elements=[host['id']],ta98Term=term,ta98Kind='exact',partOf=parent,
        description='Selectable region on the actual liver surface. The anatomical boundary remains provisional.'))
plan=dict(sourceCatalogueSha256=hashlib.sha256(raw).hexdigest(),hostGeometrySha256=mesh.geometry_sha,
    hostAttributesSha256=mesh.attributes_sha,adjoining=evidence,definitions='Exact source-facet masks and barycentric anchors; no duplicated anatomical meshes or tissue caps.',
    representation='provisional anatomical surface annotation',anatomicalAcceptance=False)
atomic_json(out/'construction.json',plan)
atlas['liverFeatureIntegration']=dict(schema='liver-surface-features/v1',hostPartId=host['id'],
    bankSha256=atlas['liverIntegration']['bankSha256'],sourcePlanSha256=digest_file(out/'construction.json'),
    termGroups=groups,provisional=True)
atomic_json(out/'candidate.json',atlas)
atomic_json(out/'evidence.json',dict(sourceCatalogueSha256=hashlib.sha256(raw).hexdigest(),
    addedSurfaceAnnotations=len(groups),addedTissueMeshes=0,sourceGeometryLiteral=True,sourceAttributesLiteral=True,
    boundaryProvisional=True,installed=False,elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(annotations=len(groups),facets={g['term']:len(g['surface']['faces'])for g in groups},elapsedSeconds=time.monotonic()-start)))
