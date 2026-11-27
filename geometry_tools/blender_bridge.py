"""Run once in an EXISTING Blender console. Local data-only render queue.

No process launch, main-file save, arbitrary queued code, or original edit.
Only owned diagnostic scenes/objects are created. Reload stops our own timer.
"""
from pathlib import Path
import hashlib,json,os,time,traceback
import bpy,bmesh,numpy as np
from mathutils import Vector, Matrix

_ga_root=Path('/Users/christiantzurcanu/Documents/dev/human-atlas').resolve()
_ga_queue=_ga_root/'.local-models/geometry-workbench/review'


def _ga_sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


_ga_implementation_sha=_ga_sha(_ga_root/'geometry_tools/blender_bridge.py')
_ga_request_id=globals().pop('atlas_geometry_review_request_id',None)
if _ga_request_id is not None and (not isinstance(_ga_request_id,str)or len(_ga_request_id)!=64 or any(c not in '0123456789abcdef'for c in _ga_request_id)):
    raise ValueError('An explicit review request SHA identity is required')
_ga_bridge_deadline=time.monotonic()+120.


def _ga_atomic(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)


def _ga_path(value,root=_ga_queue):
    p=Path(value).resolve()
    if not p.is_relative_to(root):raise ValueError('Review I/O is confined to its workbench queue')
    return p


def _ga_original_state():
    result={}
    for obj in bpy.data.objects:
        if obj.get('atlas_geometry_review'):continue
        h=hashlib.sha256()
        if obj.type=='MESH':
            v=np.empty(len(obj.data.vertices)*3,'<f4');obj.data.vertices.foreach_get('co',v);h.update(v.tobytes())
            loops=np.empty(len(obj.data.loops),'<i4');obj.data.loops.foreach_get('vertex_index',loops);h.update(loops.tobytes())
        result[str(obj.as_pointer())]=(obj.name,tuple(float(v)for row in obj.matrix_world for v in row),h.hexdigest(),tuple(m.name if m else None for m in obj.data.materials)if obj.type=='MESH'else None)
    return result


def _ga_render(task):
    started=time.monotonic();deadline=float(task.get('wall_seconds',60))
    if not 1<=deadline<=120:raise ValueError('Bounded native review deadline required')
    def check_deadline():
        if time.monotonic()-started>deadline:raise TimeoutError('Native review exceeded '+str(deadline)+' seconds; original Blender remains open')
    before=_ga_original_state();active=bpy.context.window.scene if bpy.context.window else None
    scene=bpy.data.scenes.new('AtlasGeometryReview-'+task['id'][:10]);scene['atlas_geometry_review']=True
    scene.render.engine='BLENDER_WORKBENCH';res=int(task.get('resolution',640))
    if not 128<=res<=1600:raise ValueError('Bounded render resolution required')
    scene.render.resolution_x=res;scene.render.resolution_y=res;scene.render.resolution_percentage=100;scene.render.image_settings.file_format='PNG'
    scene.display.shading.light='STUDIO';scene.display.shading.color_type='OBJECT';scene.display.shading.show_shadows=True;scene.display.shading.show_cavity=True
    scene.display.shading.cavity_type='BOTH';scene.display.shading.background_type='WORLD';scene.world=bpy.data.worlds.new(scene.name+'-world');scene.world.color=(.04,.05,.06)
    objects=[];source_rows=[];allpoints=[]
    for entry in task['meshes']:
        check_deadline()
        path=_ga_path(entry['path'])
        if _ga_sha(path)!=entry['sha256']:raise ValueError('Bound review mesh changed')
        with np.load(path,allow_pickle=False)as z:
            vertices=z['vertices'].astype(float);faces=z['faces'].copy();frame=str(z['frame'])
            words=z['attribute_normals'].copy()if 'attribute_normals'in z else None
            signs=z['attribute_normalSigns'].copy()if 'attribute_normalSigns'in z else None
            corner_words=(z['corner_normals'].copy()if 'corner_normals'in z else
                z['attribute_corner_normals'].copy()if 'attribute_corner_normals'in z else None)
        if frame!=task['frame']or vertices.ndim!=2 or vertices.shape[1]!=3 or not np.isfinite(vertices).all():raise ValueError('Invalid or differently framed review mesh')
        data=bpy.data.meshes.new(entry['name']);data.from_pydata(vertices.tolist(),[],faces.tolist());data.update();obj=bpy.data.objects.new(entry['name'],data);obj['atlas_geometry_review']=True;scene.collection.objects.link(obj);obj.color=entry['color'];objects.append(obj);allpoints.extend(vertices.tolist())
        if words is not None:
            normals=words.astype(float);normals*=signs[:,None]if signs is not None else 1;length=np.linalg.norm(normals,axis=1);normals=np.divide(normals,length[:,None],out=np.zeros_like(normals),where=length[:,None]>0)
            if np.all(length>0):
                for p in data.polygons:p.use_smooth=True
                data.normals_split_custom_set_from_vertices(normals.tolist())
                attribute=data.attributes.new(name='atlas_review_source_normal',type='FLOAT_VECTOR',domain='POINT')
                attribute.data.foreach_set('vector',normals.ravel().tolist())
        if corner_words is not None:
            if corner_words.shape!=(len(faces),3,3):raise ValueError('One literal normal word required for every face corner')
            corner_normals=corner_words.astype(float).reshape(-1,3);length=np.linalg.norm(corner_normals,axis=1)
            if not np.isfinite(corner_normals).all()or np.any(length<=0):raise ValueError('Finite nonzero corner normals required')
            corner_normals/=length[:,None]
            for p in data.polygons:p.use_smooth=True
            data.normals_split_custom_set(corner_normals.tolist())
            attribute=data.attributes.new(name='atlas_review_corner_normal',type='FLOAT_VECTOR',domain='CORNER')
            attribute.data.foreach_set('vector',corner_normals.ravel().tolist())
        actual=np.array([v.co[:]for v in data.vertices]);source_rows.append(dict(path=str(path),sha256=entry['sha256'],geometry_sha256=entry['geometry_sha256'],maximum_blender_coordinate_error_m=float(np.linalg.norm(actual-vertices,axis=1).max())))
    points=np.asarray(allpoints,float);center=(points.min(0)+points.max(0))/2;span=float(np.linalg.norm(np.ptp(points,axis=0)))
    camera_data=bpy.data.cameras.new(scene.name+'-camera');camera=bpy.data.objects.new(scene.name+'-camera',camera_data);camera['atlas_geometry_review']=True;scene.collection.objects.link(camera);scene.camera=camera;camera_data.type='ORTHO';camera_data.ortho_scale=span*1.15;camera_data.clip_start=max(span*1e-5,1e-7);camera_data.clip_end=max(10,span*10)
    images=[];folder=_ga_path(task['output'])
    try:
        for index,view in enumerate(task['views']):
            check_deadline()
            if index>=12:raise ValueError('At most 12 views per review request')
            direction=Vector(view['direction']).normalized();camera.location=Vector(center)+direction*max(span*2,1)
            # Atlas anatomy uses Y for superior, whereas Blender's generic
            # track quaternion assumes world Z is up. Specify the actual basis
            # so side/oblique diagnostic images do not roll the body sideways.
            up=Vector(view.get('up',[0,1,0])).normalized()
            if abs(up.dot(direction))>.999999:
                up=Vector((0,0,-1))if abs(direction.z)<.99 else Vector((1,0,0))
            right=up.cross(direction).normalized();up=direction.cross(right).normalized()
            camera.rotation_euler=Matrix((right,up,direction)).transposed().to_quaternion().to_euler()
            originals=[]
            if view.get('section'):
                plane=view['section'];normal=Vector(plane['normal']);point=Vector(plane['point'])
                if not normal.length:raise ValueError('Nonzero diagnostic section normal required')
                for obj in objects:
                    originals.append((obj,obj.data));mesh=obj.data.copy();obj.data=mesh;bm=bmesh.new();bm.from_mesh(mesh)
                    bmesh.ops.bisect_plane(bm,geom=list(bm.verts)+list(bm.edges)+list(bm.faces),dist=0.,plane_co=point,plane_no=normal,clear_outer=bool(plane.get('clear_outer',True)),clear_inner=not bool(plane.get('clear_outer',True)))
                    bm.to_mesh(mesh);bm.free();mesh.update()
                    # BMesh does not retain the custom split-normal layer.
                    # A point-domain vector attribute is interpolated at cut
                    # edges; restore it so sections review the emitted shading
                    # rather than a newly generated faceted substitute.
                    corner_attribute=obj.data.attributes.get('atlas_review_corner_normal')
                    original_attribute=obj.data.attributes.get('atlas_review_source_normal')
                    if corner_attribute is not None:
                        values=np.array([v.vector[:]for v in corner_attribute.data],float).reshape(-1,3)
                        lengths=np.linalg.norm(values,axis=1)
                        if len(values)!=len(mesh.loops)or not np.isfinite(values).all()or np.any(lengths<=0):
                            raise ValueError('Diagnostic clipping lost corner-normal interpolation')
                        for p in mesh.polygons:p.use_smooth=True
                        mesh.normals_split_custom_set((values/lengths[:,None]).tolist())
                    elif original_attribute is not None:
                        values=np.array([v.vector[:]for v in original_attribute.data],float).reshape(-1,3)
                        lengths=np.linalg.norm(values,axis=1)
                        if len(values)!=len(mesh.vertices)or not np.isfinite(values).all()or np.any(lengths<=0):
                            raise ValueError('Diagnostic clipping lost source normal interpolation')
                        if len(values):
                            for p in mesh.polygons:p.use_smooth=True
                            mesh.normals_split_custom_set_from_vertices((values/lengths[:,None]).tolist())
            path=folder/('view-'+str(index)+'.png');scene.render.filepath=str(path);bpy.ops.render.render(write_still=True,scene=scene.name)
            check_deadline()
            images.append(dict(path=str(path),sha256=_ga_sha(path),name=view.get('name',str(index)),section=view.get('section'),diagnostic_section_uncapped=bool(view.get('section'))))
            for obj,old in originals:
                temp=obj.data;obj.data=old;bpy.data.meshes.remove(temp)
    finally:
        if active is not None and bpy.context.window:bpy.context.window.scene=active
    if _ga_original_state()!=before:raise ValueError('Original Blender geometry or object transforms changed')
    return dict(request_id=task['id'],state='rendered',scene=scene.name,images=images,sources=source_rows,original_objects_preserved=True,main_file_saved=False,bridge_implementation_sha256=_ga_implementation_sha,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-started,wall_seconds=deadline)


def _ga_step():
    try:
        if _ga_request_id is None or time.monotonic()>=_ga_bridge_deadline:
            _ga_atomic(_ga_queue/'bridge.json',dict(state='stopped',reason='no explicit request or finite session deadline',pid=os.getpid(),existing_Blender=True))
            return None
        # Never replay stale pending work just because Blender was restarted.
        # The bootstrap selects exactly one reviewed data request.
        pending=list((_ga_queue/'pending').glob(_ga_request_id+'.json'))
        if pending:
            p=pending[0];envelope=json.loads(p.read_text());request=_ga_path(envelope['request']);task=json.loads(request.read_text());result_path=_ga_path(task['output'])/'result.json'
            try:
                if _ga_sha(request)!=envelope['sha256']:raise ValueError('Review request changed after enqueue')
                if task.get('action')!='render':raise ValueError('Only data-bound render requests are accepted')
                if not result_path.exists():_ga_atomic(result_path,_ga_render(task))
            except Exception as e:_ga_atomic(result_path,dict(request_id=task.get('id'),state='failed',error=repr(e),traceback=traceback.format_exc(),anatomical_acceptance=False))
            done=_ga_queue/'done'/p.name;done.parent.mkdir(parents=True,exist_ok=True);os.replace(p,done)
            _ga_atomic(_ga_queue/'bridge.json',dict(state='stopped',reason='selected request completed',request_id=_ga_request_id,pid=os.getpid(),existing_Blender=True))
            return None
        heartbeat=_ga_queue/'bridge.json'
        if not heartbeat.exists()or time.time()-heartbeat.stat().st_mtime>5:
            _ga_atomic(heartbeat,dict(state='ready',updated_at=time.time(),pid=os.getpid(),existing_Blender=True,implementation_sha256=_ga_implementation_sha))
    except Exception as e:print('Atlas geometry review bridge:',repr(e),flush=True)
    return .5


_ga_previous=bpy.app.driver_namespace.get('atlas_geometry_review_timer')
if _ga_previous and bpy.app.timers.is_registered(_ga_previous):bpy.app.timers.unregister(_ga_previous)
bpy.app.driver_namespace['atlas_geometry_review_timer']=_ga_step
bpy.app.timers.register(_ga_step,first_interval=.2)
print('Atlas geometry review ready for one explicit request;120s session limit; original scene preserved.',flush=True)
