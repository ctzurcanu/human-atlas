"""Read-only exact-coordinate donor surface diagnostics, without self crossing tests."""
import argparse, collections, hashlib, itertools, json, re, sys, time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

def main():
 p=argparse.ArgumentParser();p.add_argument('--start',type=int,default=0);p.add_argument('--stop',type=int,default=1140);p.add_argument('--summarize',action='store_true');p.add_argument('--interfaces',action='store_true');a=p.parse_args()
 base=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';src=base/'source';out=base/'source-topology-v1';out.mkdir(exist_ok=True)
 assets=json.loads((src/'assets.json').read_text())['assets'];instances=json.loads((src/'instances.json').read_text())['instances'];bygeometry=collections.defaultdict(list)
 for i in instances:bygeometry[i['geometryId']].append(i)
 if a.interfaces:
  started=time.monotonic();rows={r['geometryId']:r for f in sorted(out.glob('batch-*.json'))for r in json.loads(f.read_text())['assets']};positions=np.memmap(src/'positions.f32',dtype='<f4',mode='r');buckets=collections.defaultdict(list);instance_lookup={i['id']:i for i in instances}
  for asset in assets:
   pr=asset['positions'];v=np.asarray(positions[pr['byteOffset']//4:(pr['byteOffset']+pr['byteLength'])//4]).reshape(-1,3).astype(float);q=np.unique(v,axis=0)
   for instance in bygeometry[asset['geometryId']]:
    matrix=np.asarray(instance['worldMatrix']).reshape((4,4),order='F')
    for number,loop in enumerate(rows[asset['geometryId']]['boundary_loops']):
     points=q[loop['welded_vertex_ids']]@matrix[:3,:3].T+matrix[:3,3];points=points[np.lexsort(points.T[::-1])];signature=(len(points),points.astype('<f8').tobytes());buckets[signature].append(dict(source_id=instance['id'],source_name=instance['name'],geometryId=asset['geometryId'],loop_index=number,vertices=len(points)))
  matches=[owners for owners in buckets.values()if len({r['source_id']for r in owners})>1];major_ids={i['id']for r in rows.values()if r['major_trunk']for i in r['source_instances']};major_matches=[owners for owners in matches if any(r['source_id']in major_ids for r in owners)]
  result=dict(schema='biodigital-exact-boundary-loop-coincidences/v1',world_transform_applied_once=True,coordinate_rounding=False,proper_boundary_cycles=sum(len(v)for v in buckets.values()),exact_loop_pointset_coincidence_groups=len(matches),major_trunk_coincidence_groups=len(major_matches),matches=matches,major_trunk_matches=major_matches,policy='Exact world-coordinate cycle pointset equality only, no nearest-point inference. An equal pointset establishes coincidence, not opposing edge order, shared faces or a patent junction. No equality does not establish separation or intersection.',selfIntersections='UNTESTED',elapsed_seconds=time.monotonic()-started)
  (out/'exact-loop-coincidences.json').write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps({k:v for k,v in result.items()if k not in ['matches','major_trunk_matches']}));return
 if a.summarize:
  rows=[r for f in sorted(out.glob('batch-*.json')) for r in json.loads(f.read_text())['assets']];assert len(rows)==1140 and len({r['geometryId']for r in rows})==1140
  totals={k:sum(r[k]for r in rows)for k in ['triangles','zero_area_faces','welded_vertices','proper_boundary_edges','proper_nonmanifold_edges','proper_winding_conflicts','proper_face_components','simple_boundary_loops','noncycle_boundary_components']}
  summary=dict(uniqueAssets=1140,totals=totals,assetsWithZeroAreaFaces=sum(bool(r['zero_area_faces'])for r in rows),assetsWithOpenProperBoundaries=sum(bool(r['proper_boundary_edges'])for r in rows),assetsWithNoProperBoundaryEdges=sum(not r['proper_boundary_edges']for r in rows),assetsWithProperNonmanifoldEdges=sum(bool(r['proper_nonmanifold_edges'])for r in rows),selfIntersections='UNTESTED',surfaceMeaning='Exact-coordinate surface diagnostics; a closed edge surface alone does not establish a finite wall or patent lumen',topComplexAssets=sorted(rows,key=lambda r:r['proper_boundary_edges'],reverse=True)[:15],majorTrunks=[r for r in rows if r.get('major_trunk')])
  (out/'summary.json').write_text(json.dumps(summary,separators=(',',':'))+'\n');print(json.dumps({k:v for k,v in summary.items()if k not in ['topComplexAssets','majorTrunks']}));return
 started=time.monotonic();rows=[]
 positions=np.memmap(src/'positions.f32',dtype='<f4',mode='r');indices=np.memmap(src/'indices.u32',dtype='<u4',mode='r')
 for asset in assets[a.start:a.stop]:
  pr,ir=asset['positions'],asset['indices'];v=np.asarray(positions[pr['byteOffset']//4:(pr['byteOffset']+pr['byteLength'])//4]).reshape(-1,3).astype(float);f=np.asarray(indices[ir['byteOffset']//4:(ir['byteOffset']+ir['byteLength'])//4]).reshape(-1,3)
  q,inv=np.unique(v,axis=0,return_inverse=True);wf=inv[f];tri=q[wf];zero=np.flatnonzero(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)==0);keep=np.ones(len(f),bool);keep[zero]=False;pf=wf[keep]
  e=np.concatenate([pf[:,[0,1]],pf[:,[1,2]],pf[:,[2,0]]]);ue,ein,counts=np.unique(np.sort(e,axis=1),axis=0,return_inverse=True,return_counts=True);boundary=ue[counts==1];direction=np.where(e[:,0]<e[:,1],1,-1);sums=np.bincount(ein,weights=direction,minlength=len(ue));winding=int(np.sum((counts==2)&(sums!=0)))
  used=np.unique(pf);remap=np.full(len(q),-1,int);remap[used]=np.arange(len(used));graph=coo_matrix((np.ones(len(ue)),(remap[ue[:,0]],remap[ue[:,1]])),shape=(len(used),len(used)));nc=connected_components(graph,directed=False,return_labels=False)if len(used)else 0
  adj=collections.defaultdict(set)
  for x,y in boundary:adj[int(x)].add(int(y));adj[int(y)].add(int(x))
  visited=set();loops=[];noncycles=0
  for start in sorted(adj):
   if start in visited:continue
   stack=[start];vs=[]
   while stack:
    x=stack.pop()
    if x in visited:continue
    visited.add(x);vs.append(x);stack.extend(adj[x]-visited)
   if all(len(adj[x])==2 for x in vs):
    pts=q[vs];loops.append(dict(welded_vertex_ids=sorted(vs),vertex_count=len(vs),centroid_native=pts.mean(0).tolist(),bounds_native=[pts.min(0).tolist(),pts.max(0).tolist()]))
   else:noncycles+=1
  ins=bygeometry[asset['geometryId']];major=any(re.fullmatch(r'(?:(?:Left|Right) )?(?:(?:Ascending|Descending|Abdominal|Thoracic) aorta|Aorta|Aortic arch|(?:Superior|Inferior) vena cava|Pulmonary trunk|(?:Left|Right) pulmonary artery|(?:Common|Internal|External) carotid artery|Femoral artery|Femoral vein)',i['name'],re.I)for i in ins)
  row=dict(geometryId=asset['geometryId'],source_instances=[dict(id=i['id'],name=i['name'])for i in ins],triangles=len(f),welded_vertices=len(q),zero_area_faces=len(zero),zero_area_source_face_indices=zero.tolist(),proper_boundary_edges=len(boundary),proper_nonmanifold_edges=int(np.sum(counts>2)),proper_winding_conflicts=winding,proper_face_components=int(nc),simple_boundary_loops=len(loops),noncycle_boundary_components=noncycles,boundary_loops=loops,major_trunk=major,selfIntersections='UNTESTED',wallLumenStatus='unclassified surface; finite wall and lumen not established by this diagnostic')
  if major:
   # Native units remain undeclared: area/integral keys from the metre toolkit
   # are relabelled explicitly to avoid implying a physical scale.
   qualified=topology(Mesh(v,f,'biodigital-native-undeclared-units'));qualified['surface_area_native_units2']=qualified.pop('surface_area_m2');qualified['oriented_triangle_integral_native_units3']=qualified.pop('oriented_triangle_integral_m3');row['geometry_tools_exact_topology']=qualified
  rows.append(row)
 filename=out/f'batch-{a.start:04d}-{a.stop:04d}.json';assert not filename.exists();filename.write_text(json.dumps(dict(schema='biodigital-exact-source-topology/v1',policy='Exact equal coordinates only; no rounding, proximity weld, source edits, face deletion, self-crossing or patency acceptance. Proper-face edge diagnostics exclude zero-area faces analytically and retain all original face indices.',assets=rows,elapsed_seconds=time.monotonic()-started),separators=(',',':'))+'\n');print(json.dumps(dict(batch=str(filename),assets=len(rows),elapsed_seconds=time.monotonic()-started)))
if __name__=='__main__':main()
