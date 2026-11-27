"""Partition an actual finite arterial tree at witnessed source rings.

Reusable family operation: cut paired source polygons after Boolean junctions,
find the resulting wall cells by connectivity, and store each reciprocal annular
interface once. No blood-lumen disk or overlapping source alias is introduced.
"""
from pathlib import Path
import argparse,json,sys,time
from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim
from geometry_tools.tissue import pack,normal_words

args=argparse.ArgumentParser();args.add_argument('--version',default='v5');args.add_argument('--output',default='owned-bank-v2');args.add_argument('--side',choices=['left','right'],default='left');args=args.parse_args()
T=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/('male-'+args.side+'-lingual-finite-network-'+args.version);O=D/args.output
proof=json.loads((D/'independent-check-v1/evidence.json').read_text());assert proof['numerical_route_gate_pass']
z=np.load(D/'surfaces-and-ancestry.npz');v=z['vertices'];outer=z['outer_faces'];inner=z['inner_faces'];portdata=json.loads((D/'ports.json').read_text());names=portdata['component_names'];root_name=names[0];terms=portdata['component_terms'];frame=str(z['frame'])
face_sets=[outer,inner]+[z[f'annulus_{i}_faces']for i in range(7)];full=np.concatenate(face_sets);edges=defaultdict(list)
for fi,ff in enumerate(full):
    for a,b in zip(ff,np.roll(ff,-1)):edges[tuple(sorted((int(a),int(b))))].append(fi)
assert all(len(ff)==2 for ff in edges.values()),'The actual whole wall must be closed before partitioning'
used=np.unique(full);lookup={tuple(v[i]):int(i)for i in used};assert len(lookup)==len(used)
cuts={};blocked=set()
for component,name in enumerate(names[1:],1):
    field=np.load(D/(name+'-fields.npz'));p=field['fields'];guide=field['guide_points'];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(guide,axis=0),axis=1))];minimum=.008 if component==1 else .003
    for k in np.flatnonzero(arc>=minimum):
        if k>=len(guide)-2:break
        polygons=p[:,k*12:(k+1)*12]
        try:ids=np.array([[lookup[tuple(point)]for point in polygon]for polygon in polygons])
        except KeyError:continue
        ee=[tuple(sorted((int(a),int(b))))for polygon in ids for a,b in zip(polygon,np.roll(polygon,-1))]
        if all(e in edges for e in ee):break
    else:raise AssertionError(('No intact paired source cut after junction',name))
    assert k<len(guide)-2 and all(e in edges for e in ee),(name,k)
    blocked.update(ee);cuts[component]=dict(station=int(k),distance_from_source_start_m=float(arc[k]),outer=ids[0],inner=ids[1])
row=[];col=[]
for edge,(a,b)in edges.items():
    if edge not in blocked:row.extend([a,b]);col.extend([b,a])
number,label=connected_components(coo_matrix((np.ones(len(row)),(row,col)),shape=(len(full),len(full))).tocsr(),directed=False)
if number!=8:
    atomic_json(D/'owner-cut-connectivity-failure.json',dict(components=int(number),sizes=np.bincount(label).tolist(),cut_stations={names[k]:dict(station=q['station'],distance_m=q['distance_from_source_start_m'])for k,q in cuts.items()},outer_components_by_region={str(k):np.unique(z['outer_component'][label[:len(outer)]==k],return_counts=True)[0].tolist()for k in range(number)}));raise AssertionError(('Expected six externally rimmed wall cells plus separate main inner/outer sheets',number))
offset=np.cumsum([0]+[len(x)for x in face_sets]);owner_labels={0:int(label[offset[2]])};assert np.all(label[offset[2]:offset[4]]==owner_labels[0])
for pi,port in enumerate(portdata['ports'][2:],2):
    component=port['component'];ll=np.unique(label[offset[2+pi]:offset[3+pi]]);assert len(ll)==1;owner_labels[component]=int(ll[0])
remaining=set(range(number))-set(owner_labels.values());assert len(remaining)==2
# Main has no external terminal rim: its inner/outer sheets are disconnected
# until the six shared annuli below close its finite wall thickness.
label_owner={ll:ci for ci,ll in owner_labels.items()};label_owner.update({ll:1 for ll in remaining})
face_owner=np.array([label_owner[int(ll)]for ll in label]);surfaces={};cells={name:[]for name in names}
for component,name in enumerate(names):
    for sheet,ff,start in [('outer',outer,0),('inner',inner,len(outer))]:
        key=name+':'+sheet;surfaces[key]=ff[face_owner[start:start+len(ff)]==component];cells[name].append(dict(surface=key,side=1,neighbor=None))
for pi,port in enumerate(portdata['ports']):
    name=names[port['component']];key='external-mouth:'+str(pi);surfaces[key]=z[f'annulus_{pi}_faces'];cells[name].append(dict(surface=key,side=1,neighbor=None))
interfaces=[]
for component,cut in cuts.items():
    name=names[component];loops=boundary_loops(surfaces[name+':outer']);ol=next(r for r in loops if set(r)==set(cut['outer']));il=next(r for r in boundary_loops(surfaces[name+':inner'])if set(r)==set(cut['inner']));order=[int(np.flatnonzero(cut['outer']==i)[0])for i in ol];inner_ids=cut['inner'][order];assert set(inner_ids)==set(il);key='shared-wall-ring:'+name;surfaces[key]=annular_rim(ol,inner_ids)
    # annular_rim closes the directed outer boundary of the distal child.
    cells[name].append(dict(surface=key,side=1,neighbor=root_name if component==1 else'lingual'));parent=root_name if component==1 else'lingual';cells[parent].append(dict(surface=key,side=-1,neighbor=name));interfaces.append(dict(parent=parent,child=name,source_station=cut['station'],distance_from_source_start_m=cut['distance_from_source_start_m'],surface=key,whole_wall_cells_touch_only_at_reciprocal_thickness_annulus=True,blood_lumen_open=True))
main_used=np.unique(np.r_[surfaces['lingual:outer'],surfaces['lingual:inner']]);bank_origin=(v[main_used].min(0)+v[main_used].max(0))/2
encoding_file=D/'position-encoding-v1.json'
if encoding_file.exists():
    selection=json.loads(encoding_file.read_text())['selected'];assert selection is not None and selection['self_contact']['defect_pairs']==0;bank_origin=np.array(selection['origin_m'])
encoded=(v-bank_origin).astype('<f4').astype(float)+bank_origin
words=normal_words(encoded,surfaces);uvs={key:np.zeros((len(ff),3,2),'<f4')for key,ff in surfaces.items()};source_plan=json.loads((B/'male-lingual-source-plan-v1/inputs.json').read_text());parent_record=next(r for r in source_plan['fixed_sources']if r['role']==root_name);native=Mesh.load(ROOT/parent_record['mesh']);mask=face_owner[:len(outer)]==0;anc=z['outer_source_face'][mask];components=z['outer_component'][mask];literal=clipped=0;key=root_name+':outer'
for fi,(face,comp,source_face)in enumerate(zip(surfaces[key],components,anc)):
    if comp!=0:continue
    assert source_face<len(native.faces);tri=native.vertices[native.faces[source_face]];e=tri[1:]-tri[0];gram=e@e.T
    for corner,point in enumerate(v[face]):
        same=np.flatnonzero(np.all(tri==point,axis=1))
        if len(same):weights=np.eye(3)[same[0]];literal+=1
        else:xy=np.linalg.solve(gram,e@(point-tri[0]));weights=np.r_[1-xy.sum(),xy];assert weights.min()>-1e-8 and weights.max()<1+1e-8;clipped+=1
        if 'uvs'in native.attributes:uvs[key][fi,corner]=weights@native.attributes['uvs'][native.faces[source_face]]
        if 'normals'in native.attributes:
            n=weights@native.attributes['normals'][native.faces[source_face]].astype(float);normal=np.rint(n).astype('<i2')if len(same)else np.rint(n/np.linalg.norm(n)*32767).astype('<i2');face_normal=np.cross(v[face[1]]-v[face[0]],v[face[2]]-v[face[0]])
            assert normal@face_normal>0;words[key][fi,corner]=normal
ids={name:'TA98ARTERIAL:male:'+name+':'+args.side+':finite-wall'for name in names};ids[root_name]='ZA:External carotid artery.'+('l'if args.side=='left'else'r');result=pack(v,surfaces,cells,ids,frame,O,origin=bank_origin,corner_words=words,corner_uvs=uvs,max_pairs=250000,max_tests=2000000,checkpoint_directory=D/'owner-validation-checkpoints-v1');assert result['passed'],result['checks']['failures']
evidence=dict(schema='left-lingual-five-term-owned-family/v1',source_geometry_sha256=Mesh.load(D/'finite-wall.npz').geometry_sha,finite_cells=7,represented_unique_terms=sorted(set(x for x in terms if x)),names=names,term_by_component=terms,part_ids=ids,interfaces=interfaces,all_cells_closed_self_crossing_free_disjoint=True,bank_sha256=result['bank_sha256'],bank_bytes=result['total_bytes'],native_parent_literal_UV_normal_corners=literal,clipped_parent_barycentric_corners=clipped,full_patency_and_tissue_gate_passed=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-T)
atomic_json(D/'owned-bank-stage.json',evidence);print(json.dumps(evidence))
