"""Native nerve–membrane crossings define bounded clearance-opening masks.

This reads literal triangle intersections, not an assumed nerve interior.
Damaged source nerves remain untouched. Mask caps are construction data;
they must never become anatomy or close the resulting root openings.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.distance import exact_intersection_segment
from geometry_tools.audit import audit

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--begin',type=int,default=0);ap.add_argument('--count',type=int,default=10);args=ap.parse_args()
start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
sources=base/'female-root-mouths-v2/root-mouths.json';data=json.loads(sources.read_text());groups=defaultdict(list)
for r in data['roots']:groups[(r['level'],'left'if':l_'in r['part_id']else'right')].append(r)
coats=base/'female-common-coats-v1/fields.npz';z=np.load(coats);frame='model:local-female-ta98:assembled-metres'
carriers=[Mesh(z['fields'][i],z['faces'],frame)for i in [0,2]];trees=[BVH(m)for m in carriers];rows=[];clearance=.00002
for level,side in sorted(groups)[args.begin:args.begin+args.count]:
    records=groups[(level,side)];vv=[];ff=[]
    for r in records:
        m=Mesh.load(r['source_mesh']);ff.append(m.faces+sum(len(v)for v in vv));vv.append(m.vertices)
    vertices,inv=np.unique(np.concatenate(vv),axis=0,return_inverse=True);faces=inv[np.concatenate(ff)];q=vertices[faces]
    faces=faces[np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0];nerve=Mesh(vertices,faces,frame);tree=BVH(nerve)
    witnesses=[];ambiguous=0
    for boundary,(carrier,ctree)in zip([0,2],zip(carriers,trees)):
        contact,summary=pair_contacts(nerve,carrier,0.,first_bvh=tree,second_bvh=ctree,max_pairs=50000,max_tests=500000)
        assert not summary['truncated'];ambiguous+=summary['coplanar_overlap_face_pairs']
        for first,second in contact['pairs'][contact['relation']=='crossing']:
            segment=exact_intersection_segment(nerve.triangles[first],carrier.triangles[second])
            if segment['relation']!='crossing':continue
            witnesses.append(dict(nerve_face=int(first),covering_boundary=boundary,covering_face=int(second),points=segment['endpoints'].tolist(),
                maximum_coordinate_roundoff_m=float(segment['coordinate_error_upper_m'].max())))
    masks=[]
    if witnesses:
        centers=np.array([np.mean(w['points'],axis=0)for w in witnesses]);neighbors=cKDTree(centers).query_pairs(.0004);adj=defaultdict(list)
        for a,b in neighbors:adj[a].append(b);adj[b].append(a)
        todo=set(range(len(witnesses)));clusters=[]
        while todo:
            first=min(todo);todo.remove(first);pending=[first];cluster=[]
            while pending:
                i=pending.pop();cluster.append(i)
                for j in adj[i]:
                    if j in todo:todo.remove(j);pending.append(j)
            clusters.append(cluster)
        for k,cluster in enumerate(clusters):
            points=np.concatenate([np.asarray(witnesses[i]['points'])for i in cluster]);center=points.mean(0);_,_,axes=np.linalg.svd(points-center,full_matrices=True)
            local=(points-center)@axes.T;middle=(local.min(0)+local.max(0))/2;center+=middle@axes;local=(points-center)@axes.T
            radius=np.maximum(np.max(np.abs(local[:,:2]),axis=0),.000015);factor=max(1.,float(np.linalg.norm(local[:,:2]/radius,axis=1).max()));radius=radius*factor/np.cos(np.pi/48)+clearance
            length=float(np.max(np.abs(local[:,2]))+.0006);theta=np.arange(48)*2*np.pi/48;ring=np.c_[radius[0]*np.cos(theta),radius[1]*np.sin(theta),np.full(48,-length)]
            v=np.r_[ring,ring+[0,0,2*length],[[0,0,-length],[0,0,length]]];f=[]
            for i in range(48):
                j=(i+1)%48;f.extend([[i,j,48+j],[i,48+j,48+i],[96,j,i],[97,48+i,48+j]])
            if np.linalg.det(axes)<0:axes[2]*=-1
            mask=Mesh(center+v@axes,np.asarray(f),frame);proof=audit(mask,require_closed=True,max_pairs=50000,max_tests=500000);assert proof['passed'],proof['failures']
            path=out/(level+'-'+side+'-'+str(k)+'.npz');mask.save(path)
            enclosed=bool(np.max(np.linalg.norm(local[:,:2]/radius,axis=1))<np.cos(np.pi/48));assert enclosed
            masks.append(dict(path=str(path),geometry_sha256=mask.geometry_sha,source_witness_indices=cluster,center_m=center.tolist(),axes=axes.tolist(),radii_m=radius.tolist(),half_length_m=length,
                clearance_m=clearance,all_crossing_points_enclosed=enclosed,proof=proof,not_an_anatomical_cap=True))
    rows.append(dict(level=level,side=side,source_part_ids=[r['part_id']for r in records],source_geometry_sha256=[r['source_geometry_sha256']for r in records],
        combined_nerve_geometry_sha256=nerve.geometry_sha,witnesses=witnesses,masks=masks,ambiguous_coplanar_pairs=ambiguous,anatomical_interface_certified=False))
    print(json.dumps(dict(level=level,side=side,proper_crossings=len(witnesses),ports=len(masks),ambiguous_coplanar_pairs=ambiguous,elapsed_seconds=time.monotonic()-start)),flush=True)
atomic_json(out/'ports.json',dict(source_report_sha256=digest_file(sources),covering_fields_sha256=digest_file(coats),roots=rows,
    construction_only=True,final_meningeal_ports=False,construction_limits=[
        'These are native contact clusters, including individual rootlets; their count is not a count of anatomical dural openings.',
        'Root sleeves must gather the rootlets into the appropriate dorsal/ventral exits; these masks must not be directly installed as separate dural pores.',
        'No root attachment, source nerve repair or completed sleeve is certified.'],
    anatomy_changed=False,skeleton_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start))
