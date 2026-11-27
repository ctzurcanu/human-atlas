"""Coarse local reference analysis; never exports vendor mesh arrays.

Only named accessors are mapped, in memory. Deliverables contain topology
counts, gross dimensions rounded to 0.5 mm, and at most low-order moments.
They are reference evidence for independently authored geometry, not an
atlas geometry bank or permission to redistribute the supplied models.
"""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import struct
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.spatial import ConvexHull
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

EXT=Path('/Users/christiantzurcanu/Documents/dev/3d/ext')
GLB=EXT/'models/Male_2024.glb'
META=EXT/'meshes/629250ec501a0ed64e00bf2e632c09a4/assets.json'
SOURCE=ROOT/'.local-models/geometry-workbench/tongue-completion/source-v1/meshes'
OUT=ROOT/'.local-models/geometry-workbench/tongue-completion/reference-shape-v1'
WANTED={22838,22825,23844,24094,16478,20150,17818,23706,24587,22122,20872,19491,20769}
START=time.monotonic();DEADLINE=START+45


def check():
    if time.monotonic()>DEADLINE:raise TimeoutError('Reference analysis exceeded 45 second native cap')


def rounded_mm(values, quantum=.5):
    a=np.asarray(values,float);return (np.round(a/quantum)*quantum).tolist()


class SelectedGlb:
    def __init__(self,path):
        self.stream=path.open('rb');head=self.stream.read(20)
        magic,version,total,n,kind=struct.unpack('<4sIIII',head)
        assert magic==b'glTF'and version==2 and kind==0x4e4f534a
        raw=self.stream.read(n);self.json_sha=hashlib.sha256(raw).hexdigest();self.g=json.loads(raw)
        length,kind=struct.unpack('<II',self.stream.read(8));assert kind==0x004e4942
        self.bin_offset=self.stream.tell();self.mapped_bytes=0
        self.parents={child:i for i,node in enumerate(self.g['nodes'])for child in node.get('children',[])}
        self.by_name={n.get('name'):i for i,n in enumerate(self.g['nodes'])};self.world_cache={}
    def accessor(self,i):
        a=self.g['accessors'][i]
        if 'sparse'in a:raise ValueError('Sparse accessors need an explicit reader')
        b=self.g['bufferViews'][a['bufferView']];assert b.get('buffer',0)==0
        dtype=np.dtype({5120:'i1',5121:'u1',5122:'<i2',5123:'<u2',5125:'<u4',5126:'<f4'}[a['componentType']])
        width={'SCALAR':1,'VEC2':2,'VEC3':3,'VEC4':4}[a['type']]
        stride=b.get('byteStride',dtype.itemsize*width);length=(a['count']-1)*stride+dtype.itemsize*width
        self.stream.seek(self.bin_offset+b.get('byteOffset',0)+a.get('byteOffset',0));raw=self.stream.read(length)
        self.mapped_bytes+=len(raw)
        return np.ndarray((a['count'],width),dtype=dtype,buffer=raw,strides=(stride,dtype.itemsize)).copy()
    def world(self,i):
        if i in self.world_cache:return self.world_cache[i]
        node=self.g['nodes'][i]
        if 'matrix'in node:local=np.array(node['matrix']).reshape(4,4,order='F')
        else:
            local=np.eye(4);x,y,z,w=node.get('rotation',[0,0,0,1])
            r=np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
            local[:3,:3]=r@np.diag(node.get('scale',[1,1,1]));local[:3,3]=node.get('translation',[0,0,0])
        result=self.world(self.parents[i])@local if i in self.parents else local
        self.world_cache[i]=result;return result
    def mesh(self,name):
        i=self.by_name[name];node=self.g['nodes'][i];world=self.world(i);vs=[];fs=[];count=0;accessors=[]
        for primitive in self.g['meshes'][node['mesh']]['primitives']:
            mode=primitive.get('mode',4)
            pos=self.accessor(primitive['attributes']['POSITION']).astype(float)
            v=(np.c_[pos,np.ones(len(pos))]@world.T)[:,:3]
            indices=self.accessor(primitive['indices']).reshape(-1).astype(int)
            if np.any(indices>=len(v)):raise ValueError('Out-of-range index or undeclared primitive restart requires explicit handling')
            if mode==4:f=indices.reshape(-1,3)
            elif mode==5:
                f=np.column_stack([indices[:-2],indices[1:-1],indices[2:]])
                f[1::2]=f[1::2][:,[1,0,2]]
            elif mode==6:f=np.column_stack([np.repeat(indices[0],len(indices)-2),indices[1:-1],indices[2:]])
            else:raise ValueError(f'Primitive mode {mode} is not a triangle surface')
            # Strip joins often repeat controls. Remove only literal zero-area
            # triangles, retaining original alternating parity across joins.
            tri=v[f];valid=np.any(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])!=0,axis=1)
            degenerates=int((~valid).sum());f=f[valid]
            if np.linalg.det(world[:3,:3])<0:f=f[:,[0,2,1]]
            vs.append(v);fs.append(f+count);count+=len(v)
            accessors.append(dict(position=primitive['attributes']['POSITION'],indices=primitive['indices'],
                primitive_mode=mode,index_count=len(indices),actual_triangle_count=len(f),literal_zero_area_strip_triangles_removed=degenerates))
        return np.vstack(vs),np.vstack(fs),dict(node_name=name,node_index=i,mesh_index=node['mesh'],
            accessors=accessors,world_transform_applied=True,world_determinant=float(np.linalg.det(world[:3,:3])))


def mesh_summary(v,f):
    t=topology(Mesh(v,f,'local-reference:vendor-world-units'))
    points,inv=np.unique(v,axis=0,return_inverse=True);ff=inv[f]
    edges=np.concatenate([ff[:,[0,1]],ff[:,[1,2]],ff[:,[2,0]]]);e,c=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
    boundary=e[c==1];xyz=points[boundary].reshape(-1,3)if len(boundary)else np.empty((0,3))
    extent=np.ptp(points[np.unique(ff)],axis=0)
    plane_fraction=float(np.mean(abs(xyz[:,0])<max(1e-5,extent[0]*.005)))if len(xyz)else None
    out=dict(triangles=len(f),connected_components=len(t['components']),component_triangle_counts=[r['faces']for r in t['components']],
        open_edges=t['open_edges'],nonmanifold_edges=t['nonmanifold_edges'],inconsistent_winding_edges=t['inconsistent_winding_edges'],
        closed_oriented_manifold=t['closed_oriented_manifold'],euler=[r['euler']for r in t['components']],
        open_boundary_points_within_midline_band_fraction=None if plane_fraction is None else round(plane_fraction,2),
        volumetric_interpretation='Potential closed envelope; no anatomical independence established'if t['closed_oriented_manifold']else'Open surface cannot establish a distinct occupied muscle volume')
    return out


def sections(v,f,scale_mm,axis=2):
    tri=v[f];lower,upper=v[np.unique(f),axis].min(),v[np.unique(f),axis].max();axes=[i for i in range(3)if i!=axis]
    rows=[]
    for fraction in np.linspace(.05,.95,10):
        check();plane=lower+(upper-lower)*fraction;segments=[]
        for triangle in tri:
            h=triangle[:,axis]-plane
            if min(h)>=0 or max(h)<=0:continue
            points=[]
            for i,j in [(0,1),(1,2),(2,0)]:
                if h[i]*h[j]<0:points.append(triangle[i]+(triangle[j]-triangle[i])*h[i]/(h[i]-h[j]))
            if len(points)==2:segments.append(np.array(points)[:,axes]*scale_mm)
        row=dict(fraction=round(float(fraction),2),long_axis=['X','Y','Z'][axis])
        if not segments:rows.append(dict(row,intersects=False));continue
        seg=np.array(segments);a,b=seg[:,0],seg[:,1];length=np.linalg.norm(b-a,axis=1);total=length.sum()
        mean=np.einsum('i,ij->j',length,(a+b)/2)/total
        second=np.zeros((2,2))
        for aa,bb,l in zip(a,b,length):second+=l*((np.outer(aa,aa)+np.outer(bb,bb))/3+(np.outer(aa,bb)+np.outer(bb,aa))/6)
        covariance=second/total-np.outer(mean,mean);ev,vec=np.linalg.eigh(covariance);order=np.argsort(ev)[::-1];ev=ev[order];vec=vec[:,order]
        endpoints=seg.reshape(-1,2);unique,alias=np.unique(np.round(endpoints,7),axis=0,return_inverse=True)
        degree=np.bincount(alias,minlength=len(unique));closed=bool(np.all(degree==2));components=0;adj=defaultdict(set)
        for u,w in alias.reshape(-1,2):adj[int(u)].add(int(w));adj[int(w)].add(int(u))
        seen=set()
        for seed in adj:
            if seed in seen:continue
            components+=1;stack=[seed]
            while stack:
                q=stack.pop()
                if q in seen:continue
                seen.add(q);stack.extend(adj[q]-seen)
        hull=ConvexHull(endpoints)if len(unique)>=3 else None
        row.update(intersects=True,curve_components=components,section_curves_closed=closed,
            section_open_endpoints=int(np.sum(degree==1)),extent_mm=rounded_mm(np.ptp(endpoints,axis=0)),
            boundary_rms_ellipse_diameters_mm=rounded_mm(2*np.sqrt(np.maximum(2*ev,0))),
            major_axis_angle_degrees=float(np.round(np.degrees(np.arctan2(vec[1,0],vec[0,0]))/5)*5),
            hull_area_upper_bound_mm2=None if hull is None else float(np.round(hull.volume/5)*5),
            actual_occupied_area_mm2=None,
            area_note='Hull/ellipse describe a curve footprint, never establish muscle occupancy or close an open section.')
        rows.append(row)
    return rows


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'assessment.json').exists():raise FileExistsError('Keep the completed evidence stage immutable')
    language=json.loads((EXT/'meshes/en.json').read_text());names={r['cid']:r['name']for r in language['items']if r.get('cid')in WANTED}
    meta=json.loads(META.read_text());reader=SelectedGlb(GLB);meshes={};items=[]
    for index,structure in enumerate(meta['structures']):
        if structure.get('cid')not in WANTED:continue
        for source_node in structure.get('sni',[]):
            if source_node['id']not in reader.by_name:continue
            check();v,f,identity=reader.mesh(source_node['id']);key=f"{structure['cid']}:{structure.get('attr')}:{index}"
            meshes[key]=(v,f);items.append(dict(key=key,name=names[structure['cid']],cid=structure['cid'],attr=structure.get('attr'),
                structure_index=index,linked_structure_indices=structure.get('linked',[]),twin_structure_index=structure.get('twin'),
                identity=identity,topology=mesh_summary(v,f)))
    reference_hyoid=next(meshes[r['key']][0]for r in items if r['cid']==20872)
    reference_mandible=next(meshes[r['key']][0]for r in items if r['cid']==22122)
    native_hyoid=np.vstack([Mesh.load(SOURCE/f'native-{i:02d}.npz').vertices for i in [13,14,15,16,17]])
    native_mandible=Mesh.load(SOURCE/'native-05.npz').vertices
    ratios=[]
    for label,ref,native in [('hyoid',reference_hyoid,native_hyoid),('mandible',reference_mandible,native_mandible)]:
        scales=1000*np.ptp(native,axis=0)/np.ptp(ref,axis=0)
        ratios.append(dict(bone=label,native_extent_mm=rounded_mm(np.ptp(native,axis=0)*1000),
            implied_mm_per_vendor_unit_xyz=np.round(scales,1).tolist()))
    # Gross transverse skeleton dimensions set one explicitly provisional
    # uniform conversion. Axial/AP ratios remain reported as an uncertainty.
    scale_mm=float(np.mean([1000*np.ptp(native_hyoid[:,0])/np.ptp(reference_hyoid[:,0]),
        1000*np.ptp(native_mandible[:,0])/np.ptp(reference_mandible[:,0])]))
    all_sections={}
    for row in items:
        v,f=meshes[row['key']];row['gross_extent_mm']=rounded_mm(np.ptp(v[np.unique(f)],axis=0)*scale_mm)
        row['gross_axis_names']=['left-right X','superior-inferior Y','anterior-posterior Z']
        if row['cid']not in{22122,20872,20769}:
            axis=2 if row['cid']in{22838,22825,23844,24094,19491}else 1
            all_sections[row['key']]=sections(v,f,scale_mm,axis)
    report=dict(schema='local-tongue-reference-coarse-shape-analysis/v1',
        source=dict(glb=str(GLB),json_sha256=reader.json_sha,metadata=str(META),
            metadata_sha256=hashlib.sha256(META.read_bytes()).hexdigest(),asset=reader.g.get('asset')),
        rights=dict(open_source_license_established=False,local_reference_only=True,
            vendor_geometry_arrays_exported=False,vendor_normals_indices_or_textures_exported=False),
        calibration=dict(provisional_mm_per_vendor_unit=round(scale_mm,1),method='Mean scale from actual native and reference hyoid/mandible transverse dimensions; all three dimensional ratios are reported.',
            bone_cross_checks=ratios,not_a_rigid_anatomical_registration=True,metric_rounding_mm=.5),
        primary_research=dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
            constraints=['Tongue muscle fascicles interweave; complete muscles are not mutually exclusive slabs.',
                'Published atlas traces transverse and vertical as one combined region where intertwined.',
                'Distinguish an external muscle envelope from its distributed intralingual fascicles.']),
        items=items,section_metrics_file=str(OUT/'cross-sections.json'),
        limitations=['Topology does not establish anatomical correctness or a redistribution license.',
            'Open section footprints and convex hulls must not be interpreted as occupied solid tissue.',
            'These rounded gross measurements constrain independent authorship; no precise donor profile is saved.'],
        selected_accessor_bytes_read=reader.mapped_bytes,elapsed_seconds=time.monotonic()-START)
    (OUT/'assessment.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'cross-sections.json').write_text(json.dumps(all_sections,indent=2)+'\n')
    print(json.dumps(dict(report=str(OUT/'assessment.json'),sections=str(OUT/'cross-sections.json'),
        scale=report['calibration'],items=[dict(name=r['name'],side=r['attr'],extent=r['gross_extent_mm'],topology=r['topology'])for r in items],
        elapsed_seconds=report['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
