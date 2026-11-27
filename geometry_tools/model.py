"""The single atlas reader: indexed controls, origins and offsets exactly once."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from collections import OrderedDict
from functools import cached_property
import gzip, json
import numpy as np
from .common import ROOT, digest_file, identity, atomic_npz


@dataclass
class Mesh:
    vertices: np.ndarray
    faces: np.ndarray
    frame: str
    attributes: dict[str, np.ndarray] = field(default_factory=dict)
    source: dict = field(default_factory=dict)

    def __post_init__(self):
        self.vertices = np.array(self.vertices, dtype='<f8',copy=True)
        raw = np.asarray(self.faces)
        if raw.ndim != 2 or raw.shape[1] != 3 or not np.issubdtype(raw.dtype, np.integer):
            raise ValueError('Integer triangle indices with shape [faces,3] required')
        self.faces = raw.astype('<i8', copy=True)
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 3 or not np.isfinite(self.vertices).all():
            raise ValueError('Finite metre coordinates with shape [vertices,3] required')
        if not isinstance(self.frame, str) or not self.frame.strip():
            raise ValueError('An explicit coordinate frame is required')
        if len(raw) and (raw.min() < 0 or raw.max() >= len(self.vertices)):
            raise ValueError('Triangle index outside its actual vertex bank')
        for name, values in self.attributes.items():
            if len(values) != len(self.vertices):
                raise ValueError(f'{name}: attribute does not cover the compact vertex bank')
        self.vertices.flags.writeable=False;self.faces.flags.writeable=False
        self.attributes={name:np.array(values,copy=True)for name,values in self.attributes.items()}
        for values in self.attributes.values():values.flags.writeable=False

    @cached_property
    def triangles(self):
        return self.vertices[self.faces]

    @cached_property
    def bounds(self):
        if not len(self.faces):
            raise ValueError('An empty mesh has no triangle bounds')
        points = self.vertices[np.unique(self.faces)]
        return np.array([points.min(0), points.max(0)])

    @cached_property
    def geometry_sha(self):
        # Exact-coordinate welding only. Keep original face and corner order:
        # BVH witnesses must still identify the correct source face/corner.
        ids=np.unique(self.faces);used,inv=np.unique(self.vertices[ids],axis=0,return_inverse=True)
        mapping=np.full(len(self.vertices),-1,dtype='<i8');mapping[ids]=inv
        return identity('indexed-world-mesh/v1', used.astype('<f8'), mapping[self.faces].astype('<i8'))

    @cached_property
    def attributes_sha(self):
        return identity('mesh-attributes/v1', *[(k, identity(np.asarray(v))) for k,v in sorted(self.attributes.items())])

    @cached_property
    def surface_sha(self):
        ids=np.unique(self.faces);v,inv=np.unique(self.vertices[ids],axis=0,return_inverse=True);mapping=np.full(len(self.vertices),-1,int);mapping[ids]=inv;f=mapping[self.faces]
        if not len(f):
            return identity(v, f)
        rotations = np.stack([np.roll(f, k, axis=1) for k in range(3)], axis=1)
        first = np.argmin(rotations[:,:,0], axis=1); f = rotations[np.arange(len(f)), first]
        return identity('oriented-surface/v1', v, f[np.lexsort(f.T[::-1])])

    def save(self, path):
        arrays = dict(vertices=self.vertices, faces=self.faces, frame=np.array(self.frame),
                      source_json=np.array(json.dumps(self.source, sort_keys=True)))
        arrays.update({'attribute_'+k:np.asarray(v) for k,v in self.attributes.items()})
        return atomic_npz(path, **arrays)

    @classmethod
    def load(cls, path, *, frame=None):
        path=Path(path).resolve()
        with np.load(path, allow_pickle=False) as z:
            v=z['vertices'] if 'vertices' in z else z['positions']
            origin=np.asarray(z['geometryOrigin'],float) if 'geometryOrigin'in z else np.zeros(3)
            actual_frame=frame or (str(z['frame']) if 'frame'in z else None)
            attrs={k[len('attribute_'):]:z[k].copy()for k in z.files if k.startswith('attribute_')}
            source=json.loads(str(z['source_json'])) if 'source_json'in z else {}
            source.update(input_path=str(path),input_sha256=digest_file(path))
            return cls(v.astype(float)+origin,z['faces'].copy(),actual_frame,attrs,source)


class AtlasReader:
    """Read literal file-backed views with bounded mmap handles and hash reuse."""
    def __init__(self, root=ROOT, max_open_banks=16):
        self.root=Path(root).resolve();self.max_open_banks=max_open_banks
        self.maps=OrderedDict();self.hashes={}

    def asset(self, url):
        if url.startswith('/local-models/'):
            path=self.root/'.local-models'/url[len('/local-models/'):]
        elif url.startswith('/') and not url.startswith('//'):
            path=self.root/'public'/url[1:]
        else:
            raise ValueError('Only existing local atlas assets are read; no remote mesh download')
        path=path.resolve()
        if not path.is_relative_to(self.root):
            raise ValueError('Asset escapes the declared project root')
        return path

    def file_sha(self, path):
        p=Path(path);s=p.stat();key=(str(p),s.st_ino,s.st_size,s.st_mtime_ns)
        if key not in self.hashes:self.hashes[key]=digest_file(p)
        return self.hashes[key]

    def _bank(self,path):
        key=str(path)
        if key in self.maps:
            self.maps.move_to_end(key);return self.maps[key]
        value=np.memmap(path,dtype='u1',mode='r');self.maps[key]=value
        while len(self.maps)>self.max_open_banks:self.maps.popitem(last=False)
        return value

    @staticmethod
    def view(bank, offset, count, dtype, width=1):
        dtype=np.dtype(dtype)
        if not isinstance(offset,int)or offset<0 or offset%dtype.itemsize or count<0 or offset+count*width*dtype.itemsize>len(bank):
            raise ValueError('Packed array is misaligned or outside the source bank')
        return np.ndarray((count,width) if width>1 else(count,),dtype=dtype,buffer=bank,offset=offset)

    def read_part(self, atlas, part, *, frame, catalogue_path=None, catalogue_sha=None):
        n=part['vertexCount'];ic=part['indexCount']
        if ic%3:raise ValueError('Triangle index count must be divisible by three')
        path=self.asset(atlas['chunks'][part['chunk']]['url']);bank=self._bank(path)
        f=self.view(bank,part['indices'],ic,'<u4').reshape(-1,3)
        if len(f)and f.max()>=n:raise ValueError('Source triangle index exceeds its declared vertex count')
        used, inv=np.unique(f,return_inverse=True)
        raw=self.view(bank,part['positions'],n,'<f4',3)[used].astype('<f8')
        scale=float(part.get('geometryScale',1.))
        if not np.isfinite(scale)or scale<=0:raise ValueError('Geometry scale must be finite and positive')
        origin=np.asarray(part.get('geometryOrigin',[0.,0.,0.]),float)
        offset=np.asarray(part.get('sourceOffset',[0.,0.,0.]),float)
        if origin.shape!=(3,)or offset.shape!=(3,)or not np.isfinite(np.r_[origin,offset]).all():
            raise ValueError('Both source origin and source offset must be finite metre vectors')
        attrs={}
        for name,dtype,width in [('normals','<i2',3),('normalSigns','i1',1),('uvs','<f4',2),('colors','u1',3)]:
            if name in part:attrs[name]=self.view(bank,part[name],n,dtype,width)[used].copy()
        if 'normalSigns'in attrs and not np.isin(attrs['normalSigns'],[-1,1]).all():
            raise ValueError('Packed normal signs must be exactly -1 or +1')
        from .geometry_instance import instance_pose
        pose=instance_pose(part)
        if pose is not None:
            rotation,instance_scale=pose
            vertices=((raw*scale+offset)@rotation.T)*instance_scale+origin
            if 'normals'in attrs:attrs['normals']=attrs['normals'].astype(float)@rotation.T
        else:vertices=raw*scale+origin+offset
        source=dict(kind='atlas-part',project_root=str(self.root),catalogue_path=str(Path(catalogue_path).resolve())if catalogue_path else None,
                    catalogue_sha256=catalogue_sha,bank_path=str(path),bank_sha256=self.file_sha(path),record=part,
                    geometry_origin=origin.tolist(),geometry_scale=scale,source_offset=offset.tolist(),source_vertex_ids=used.tolist(),
                    source_face_ids=list(range(len(f))),coordinates='indexed source Float32 * positive geometryScale + geometryOrigin + sourceOffset, exactly once; metres',
                    anatomical_registration_verified=False)
        if pose is not None:source.update(geometry_instance=part['geometryInstance'],coordinates='literal indexed source Float32 in a proper positive similarity; local source offset, rotated/scaled once, then geometryOrigin; metre world frame',normal_basis='local integer normal words rotated into the world frame without integer requantization')
        return Mesh(vertices,inv.reshape(-1,3),frame,attrs,source)

    def reload(self, source, frame):
        path=Path(source['bank_path'])
        if self.file_sha(path)!=source['bank_sha256']:
            raise ValueError('Source bank changed; reimport the catalogue before reusing cached geometry')
        # A saved slice is independent of mutable current catalogue metadata.
        record=source['record'];atlas={'chunks':[{'url':'/local-models/'+str(path.relative_to(self.root/'.local-models'))}]} if path.is_relative_to(self.root/'.local-models') else {'chunks':[{'url':'/'+str(path.relative_to(self.root/'public'))}]}
        record=dict(record,chunk=0)
        return self.read_part(atlas,record,frame=frame,catalogue_path=source.get('catalogue_path'),catalogue_sha=source.get('catalogue_sha256'))
