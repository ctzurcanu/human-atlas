"""Canonical shared tissue charts from one source-provenance Boolean graph.

Independent exports can triangulate the same interface differently. Overlay
their actual fragments in their common original face plane, once, then use
the same triangles with opposite signs. This does not infer tissue identity.
"""
from collections import defaultdict
import json
from types import SimpleNamespace
import numpy as np
import manifold3d as md
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,CommonPositionBank
from geometry_tools.algorithms.source_face_complement import conform_surface_boundary_edges


class PlanarOverlayError(ValueError):
    def __init__(self, details):
        super().__init__('Planar overlay triangulation does not replay its occupied area: '+
                         str({k:v for k,v in details.items()if k!='polygons'}))
        self.details=details


class SourcePlaneError(ValueError):
    def __init__(self,details):
        super().__init__('Boolean face left its declared source plane: '+str(details['maximum_plane_error_m']))
        self.details=details


def triangulate_overlay(polygons, expected_area, factor, tolerance_m):
    """Respect source uncertainty while rejecting flipped/overlapping ears."""
    flat=np.concatenate(polygons)
    perimeter=sum(float(np.linalg.norm(p-np.roll(p,1,axis=0),axis=1).sum())for p in polygons)/factor
    area_resolution=tolerance_m*perimeter
    if abs(expected_area)<=area_resolution:
        return None
    # epsilon=0 lets near-coincident contour points produce a large reversed
    # ear. The actual uncertainty is the same declared source intersection
    # resolution used by the 3D common bank, expressed in this planar chart.
    tri=np.asarray(md.triangulate(polygons,epsilon=tolerance_m*factor,allow_convex=False),int)
    tq=flat[tri]/factor
    actual_area=float(abs((tq[:,1,0]-tq[:,0,0])*(tq[:,2,1]-tq[:,0,1])-
                          (tq[:,1,1]-tq[:,0,1])*(tq[:,2,0]-tq[:,0,0])).sum()/2)
    floating_bound=256*np.finfo(float).eps*float(np.ptp(flat/factor,axis=0).max())**2*len(tri)
    if expected_area<0 or abs(actual_area-expected_area)>max(area_resolution,floating_bound):
        raise PlanarOverlayError(dict(expected_area_m2=expected_area,actual_area_m2=actual_area,
            area_resolution_m2=area_resolution,floating_bound_m2=floating_bound,
            factor=factor,polygons=[p.tolist()for p in polygons]))
    return flat,tri


class BooleanCells:
    def __init__(self,*,tolerance_m=2e-12):
        self.primitives={};self.tolerance_m=tolerance_m;self.aliases={}

    def source(self,vertices,faces,name):
        v=np.asarray(vertices,float);f=np.asarray(faces,int);solid,oid=make_solid(v,f)
        self.primitives[oid]=(v,f,name);return solid

    def literal_cell(self,vertices,faces,name):
        """Retain an unchanged source cell without a kernel export round-trip.

        Closure/occupied-volume construction is still checked by source().
        Only cells that undergo no Boolean operation may use this exporter.
        Every emitted triangle and source face ID remains literal.
        """
        solid=self.source(vertices,faces,name);oid=max(self.primitives);v,f,_=self.primitives[oid]
        used,inv=np.unique(f,return_inverse=True)
        mesh=SimpleNamespace(vert_properties=v[used],tri_verts=inv.reshape(-1,3).astype(np.uint64),
            run_index=np.array([0,f.size],np.uint64),run_original_id=np.array([oid],np.uint32),
            face_id=np.arange(len(f),dtype=np.uint64))
        return SimpleNamespace(to_mesh64=lambda:mesh,volume=solid.volume,is_empty=solid.is_empty)

    def derived_source(self,solid,name,*,retain_parents=False):
        m=solid.to_mesh64();result=self.source(np.asarray(m.vert_properties)[:,:3],np.asarray(m.tri_verts,int),name)
        if retain_parents:
            oid=max(self.primitives);origins,faces=provenance(m)
            for i,(parent,face)in enumerate(zip(origins,faces)):self.aliases[(oid,i)]=(int(parent),int(face))
        return result

    def save(self,solids,path):
        from .common import atomic_npz
        arrays={};shared={};meta=dict(tolerance_m=self.tolerance_m,primitives=[],cells=[],aliases=[[list(k),list(v)]for k,v in self.aliases.items()])
        def store(key,value):
            value=np.asarray(value);identity=(id(value),value.dtype.str,value.shape)
            if identity in shared:return shared[identity]
            arrays[key]=value;shared[identity]=key;return key
        for i,(oid,(v,f,name))in enumerate(self.primitives.items()):
            prefix='primitive_'+str(i)
            meta['primitives'].append(dict(oid=oid,name=name,prefix=prefix,
                v_key=store(prefix+'_v',v),f_key=store(prefix+'_f',f)))
        for i,(name,solid)in enumerate(solids.items()):
            m=solid.to_mesh64();prefix='cell_'+str(i)
            keys={k:store(prefix+'_'+k,getattr(m,k))for k in ['vert_properties','tri_verts','run_index','run_original_id','face_id']}
            meta['cells'].append(dict(name=name,prefix=prefix,array_keys=keys,volume_m3=float(solid.volume())))
        arrays['manifest_json']=np.array(json.dumps(meta));atomic_npz(path,**arrays)

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False)as z:
            meta=json.loads(str(z['manifest_json']));b=cls(tolerance_m=meta['tolerance_m'])
            cache={}
            def read(key):
                if key not in cache:cache[key]=z[key].copy()
                return cache[key]
            b.primitives={q['oid']:(read(q.get('v_key',q['prefix']+'_v')),read(q.get('f_key',q['prefix']+'_f')),q['name'])for q in meta['primitives']}
            b.aliases={tuple(k):tuple(v)for k,v in meta['aliases']};solids={}
            for q in meta['cells']:
                m=SimpleNamespace(**{k:read(q.get('array_keys',{}).get(k,q['prefix']+'_'+k))for k in ['vert_properties','tri_verts','run_index','run_original_id','face_id']})
                solids[q['name']]=SimpleNamespace(to_mesh64=lambda mesh=m:mesh,volume=lambda volume=q['volume_m3']:volume,is_empty=lambda:False)
        return b,solids

    def assemble(self,solids,*,progress=None,conform_collinear_cycles=False):
        source=np.unique(np.concatenate([v[np.unique(f)]for v,f,_ in self.primitives.values()]),axis=0)
        equivalent={};canonical={}
        for oid,(v,f,_)in self.primitives.items():
            for i,tri in enumerate(f):
                key=tuple(sorted(tuple(p)for p in v[tri]));canonical.setdefault(key,(oid,i));equivalent[(oid,i)]=canonical[key]
        for child,parent in self.aliases.items():
            while parent in self.aliases:parent=self.aliases[parent]
            equivalent[child]=equivalent[parent]
        bank=CommonPositionBank(source,self.tolerance_m);groups=defaultdict(lambda:defaultdict(list));volumes={}
        for name,solid in solids.items():
            if solid.is_empty():raise ValueError('Empty proposed tissue: '+name)
            volumes[name]=float(solid.volume());mesh=solid.to_mesh64();v=np.asarray(mesh.vert_properties)[:,:3]
            f=np.asarray(mesh.tri_verts,int);oids,parent=provenance(mesh);mapped=np.array([bank.insert(p)for p in v],int)
            for tri,p,oid,face in zip(mapped[f],v[f],oids,parent):
                if oid not in self.primitives:raise ValueError('A derived field was not rebound as a common primitive')
                oid,face=equivalent[(int(oid),int(face))]
                pv,pf,_=self.primitives[oid];q=pv[pf[int(face)]];n=np.cross(q[1]-q[0],q[2]-q[0])
                plane_error=float(np.max(abs((p-q[0])@n))/np.linalg.norm(n))
                if plane_error>self.tolerance_m*2:
                    raise SourcePlaneError(dict(cell=name,primitive=self.primitives[oid][2],source_face=int(face),
                        source_triangle_m=q.tolist(),fragment_triangle_m=p.tolist(),maximum_plane_error_m=plane_error,
                        tolerance_m=self.tolerance_m,source_double_area_m2=float(np.linalg.norm(n))))
                sign=1 if np.dot(np.cross(p[1]-p[0],p[2]-p[0]),n)>0 else-1
                groups[(int(oid),int(face))][(name,sign)].append(tri)
        pairs=defaultdict(list);overlays=0;zero_area_fragments=0;thin_overlay_fragments=0;maximum_removed_width=0.;empty_overlay_regions=0
        def record(f,a,b,*,derived_overlay=False):
            nonlocal zero_area_fragments,thin_overlay_fragments,maximum_removed_width
            if not len(f):return
            if b is not None and a>b:a,b=b,a;f=f[:,[0,2,1]]
            # Canonicalizing a newly intersected polygon can make a numerical
            # overlay ear exactly collinear or repeat one endpoint. It has no
            # occupied area. Original nonzero facets are retained. New overlay
            # ears below the explicitly declared intersection resolution may
            # also collapse; conformity/closed-cell checks still reject holes.
            xyz=np.asarray([bank.rows[int(i)]for i in np.asarray(f).ravel()]).reshape(-1,3,3)
            area2=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)
            nonzero=area2>0
            zero_area_fragments+=int((~nonzero).sum())
            if derived_overlay:
                longest=np.linalg.norm(xyz-np.roll(xyz,1,axis=1),axis=2).max(1)
                width=np.divide(area2,longest,out=np.zeros_like(area2),where=longest>0)
                thin=nonzero&(width<=self.tolerance_m)
                thin_overlay_fragments+=int(thin.sum())
                maximum_removed_width=max(maximum_removed_width,float(width[thin].max(initial=0)))
                nonzero&=~thin
            f=np.asarray(f)[nonzero]
            pairs[(a,b)].extend(f)
        for count,((oid,face),owners)in enumerate(groups.items()):
            if progress and count%5000==0:progress(count,len(groups),overlays)
            owners={k:np.asarray(v,int)for k,v in owners.items()}
            if len(owners)==1:
                (name,sign),f=next(iter(owners.items()));record(f,name,None);continue
            if len(owners)==2:
                (a,sa),(b,sb)=owners
                if sa==-sb:
                    fa,fb=owners[(a,sa)],owners[(b,sb)]
                    if {tuple(q)for q in np.sort(fa,axis=1)}=={tuple(q)for q in np.sort(fb,axis=1)}:
                        record(fa,a,b);continue
            overlays+=1;pv,pf,_=self.primitives[oid];q=pv[pf[face]]
            # An arbitrarily numbered first edge can be microscopic beside
            # the other two. Scaling the entire planar overlay by its inverse
            # overflows Clipper's range and loses narrow source fragments.
            # A longest-edge chart has the same physical plane, with every
            # source-triangle coordinate bounded by the declared chart scale.
            edge=int(np.argmax(np.linalg.norm(np.roll(q,-1,axis=0)-q,axis=1)))
            origin=q[edge];u=q[(edge+1)%3]-origin;scale=np.linalg.norm(u);u/=scale
            normal=np.cross(q[1]-q[0],q[2]-q[0]);normal/=np.linalg.norm(normal);w=np.cross(normal,u);basis=np.c_[u,w]
            factor=1e4/scale;regions={}
            for (name,sign),f in owners.items():
                actual=np.array([bank.rows[int(i)]for i in f.ravel()]).reshape(-1,3,3)
                triangles=(actual-origin)@basis*factor
                if sign<0:triangles=triangles[:,[0,2,1]]
                regions[(name,sign)]=md.CrossSection(list(triangles),md.FillRule.NonZero)
            positives={name:r for(name,s),r in regions.items()if s==1};negatives={name:r for(name,s),r in regions.items()if s==-1}
            union_plus=md.CrossSection();union_minus=md.CrossSection()
            for r in positives.values():union_plus=union_plus+r
            for r in negatives.values():union_minus=union_minus+r
            pieces=[]
            for a,ra in positives.items():
                for b,rb in negatives.items():pieces.append((a,b,ra^rb))
                pieces.append((a,None,ra-union_minus))
            for b,rb in negatives.items():pieces.append((None,b,rb-union_plus))
            for a,b,region in pieces:
                if region.is_empty():continue
                polygons=region.to_polygons()
                expected_area=float(region.area())/factor**2
                # Opposite nearly coincident contours can leave a signed
                # roundoff region. Triangulating that nominally empty region
                # may fill a large triangle across its coincident hole.
                # Discard only derived regions below the declared geometric
                # resolution, then require the real triangles to replay area.
                try:
                    triangulated=triangulate_overlay(polygons,expected_area,factor,self.tolerance_m)
                except PlanarOverlayError as e:
                    e.details.update(source_primitive=oid,source_face=face,owners=[a,b])
                    raise
                if triangulated is None:
                    empty_overlay_regions+=1;continue
                flat,tri=triangulated
                xyz=origin+(flat/factor)@basis.T;ids=np.array([bank.insert(p)for p in xyz],int);f=ids[tri]
                if a is None:a,b=b,None;f=f[:,[0,2,1]]
                record(f,a,b,derived_overlay=True)
        # Source cells can already share a literal chart while having distinct
        # primitive IDs (a papilla and its pyramid, for example). Reconcile
        # those actual identical triangles globally as well.
        ownership=defaultdict(dict)
        for (a,b),faces in pairs.items():
            for tri in faces:
                key=tuple(sorted(map(int,tri)));q=tuple(map(int,tri));sign=1 if q in [key,(key[1],key[2],key[0]),(key[2],key[0],key[1])]else-1
                for owner,s in [(a,sign),(b,-sign)]:
                    if owner is not None:ownership[key][owner]=ownership[key].get(owner,0)+s
        pairs=defaultdict(list)
        for tri,owners in ownership.items():
            owners={n:s for n,s in owners.items()if s}
            if not owners:continue
            if any(abs(s)!=1 for s in owners.values())or len(owners)>2:
                p=np.array([bank.rows[i]for i in tri]);area=float(np.linalg.norm(np.cross(p[1]-p[0],p[2]-p[0]))/2)
                raise ValueError('Repeated or multiply occupied tissue face: '+str(dict(owners=owners,area_m2=area,coordinates=p.tolist())))
            names=sorted(owners);a=names[0];b0=names[1]if len(names)>1 else None
            if b0 and owners[a]!=-owners[b0]:
                points=np.asarray([bank.rows[i]for i in tri])
                area=float(np.linalg.norm(np.cross(points[1]-points[0],points[2]-points[0]))/2)
                raise ValueError('Nonreciprocal existing source chart: '+str(dict(owners=owners,
                    area_m2=area,coordinates=points.tolist())))
            pairs[(a,b0)].append(tri if owners[a]>0 else(tri[0],tri[2],tri[1]))
        surfaces={};cells={n:[]for n in solids}
        for i,((a,b),faces)in enumerate(sorted(pairs.items(),key=lambda q:str(q[0]))):
            sid='native-interface:'+str(i);surfaces[sid]=np.asarray(faces,int)
            cells[a].append(dict(surface=sid,side=1,neighbor=b))
            if b:cells[b].append(dict(surface=sid,side=-1,neighbor=a))
        surfaces,conformity=conform_surface_boundary_edges(bank.positions,surfaces,cells,tolerance=self.tolerance_m,
            allow_collinear_cycles=conform_collinear_cycles)
        return dict(positions=bank.positions,surfaces=surfaces,cells=cells,evidence=dict(source_controls=len(source),derived_controls=len(bank.positions)-len(source),
                    maximum_intersection_snap_m=bank.maximum_snap,source_face_groups=len(groups),planar_overlays=overlays,
                    exact_zero_area_overlay_fragments_removed=zero_area_fragments,
                    derived_overlay_fragments_below_construction_resolution=thin_overlay_fragments,
                    maximum_removed_overlay_width_m=maximum_removed_width,
                    derived_overlay_regions_below_construction_resolution=empty_overlay_regions,
                    source_cell_volumes_m3=volumes,boundary_conformity=conformity))
