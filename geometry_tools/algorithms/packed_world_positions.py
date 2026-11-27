"""Inverse of the atlas position transform; apply origins and scales once."""
import numpy as np

from geometry_tools.geometry_instance import instance_pose


def encode_world_positions(points,part):
    world=np.asarray(points,float)
    scale=float(part.get('geometryScale',1.))
    origin=np.asarray(part.get('geometryOrigin',[0.,0.,0.]),float)
    offset=np.asarray(part.get('sourceOffset',[0.,0.,0.]),float)
    if world.ndim!=2 or world.shape[1]!=3 or origin.shape!=(3,)or offset.shape!=(3,):
        raise ValueError('Explicit world points and three-dimensional source origin/offset required')
    if not np.isfinite(np.r_[world.ravel(),origin,offset,scale]).all()or scale<=0:
        raise ValueError('Finite coordinates and positive geometry scale required')
    pose=instance_pose(part)
    local=world-origin
    if pose is not None:
        rotation,instance_scale=pose
        local=(local/instance_scale)@rotation
    return ((local-offset)/scale).astype('<f4')


def decode_source_positions(points,part):
    raw=np.asarray(points,float);scale=float(part.get('geometryScale',1.))
    origin=np.asarray(part.get('geometryOrigin',[0.,0.,0.]),float)
    offset=np.asarray(part.get('sourceOffset',[0.,0.,0.]),float)
    pose=instance_pose(part)
    if pose is None:return raw*scale+offset+origin
    rotation,instance_scale=pose
    return ((raw*scale+offset)@rotation.T)*instance_scale+origin
