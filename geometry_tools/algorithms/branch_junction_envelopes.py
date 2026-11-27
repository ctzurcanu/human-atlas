"""Close explicit tube branches for an outer-minus-fluid network Boolean.

Internal mouths continue into a caller-supplied point strictly inside their
parent. The conical continuation is a junction sleeve, not an end disk. A
caller must check containment, communication, normals, topology and anatomy.
External fluid disks are displaced past the native outer end solely by its
measured nonplanarity/storage guard; they must never survive in the wall.
"""
import numpy as np
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

def close_branch(outer,inner,faces,internal_anchors,origin=(0.,0.,0.),*,fluid_through_external=True):
    outer=np.asarray(outer,float);inner=np.asarray(inner,float);faces=np.asarray(faces,int);origin=np.asarray(origin,float);loops=boundary_loops(faces);result=[];evidence=[]
    for which,source in [('outer',outer),('fluid',inner)]:
        points=list(source);out=list(faces);disk_faces=[];junction_faces=[];rows=[]
        for index,loop in enumerate(loops):
            ring=source[loop];center=ring.mean(0)
            if index in internal_anchors:
                anchor=np.asarray(internal_anchors[index],float);tip=len(points);points.append(anchor);ids=[]
                for a,b in zip(loop,np.roll(loop,-1)):ids.append(len(out));out.append([int(b),int(a),tip])
                junction_faces.extend(ids);rows.append(dict(mouth=index,kind='parent-lumen-junction',anchor=anchor.tolist(),newFaceIds=ids));continue
            _,_,axes=np.linalg.svd(outer[loop]-outer[loop].mean(0),full_matrices=False);direction=axes[-1];touching=faces[np.any(np.isin(faces,loop),axis=1)];adjacent=np.setdiff1d(np.unique(touching),loop)
            if not len(adjacent):raise ValueError('A source mouth has no adjacent tube section')
            if direction@(outer[adjacent].mean(0)-center)>0:direction=-direction
            if which=='outer' or not fluid_through_external:
                tip=len(points);points.append(center);ids=[]
                for a,b in zip(loop,np.roll(loop,-1)):ids.append(len(out));out.append([int(b),int(a),tip])
                disk_faces.extend(ids);rows.append(dict(mouth=index,kind='external-'+which+'-classifier',newFaceIds=ids))
            else:
                q=outer[loop];guard=max(2e-9,4*float(np.linalg.norm(np.spacing((q-origin).astype('<f4')).astype(float),axis=1).max()));extension=max(0.,float(((q-center)@direction).max()-((ring-center)@direction).min()))+guard
                if not 0<extension<2e-6:raise ValueError('External fluid continuation exceeds a numerical end-plane guard')
                new=np.arange(len(points),len(points)+len(loop));points.extend(ring+extension*direction);tip=len(points);points.append(center+extension*direction);ids=[]
                for j,k in zip(range(len(loop)),np.roll(np.arange(len(loop)),-1)):
                    a,b,aa,bb=int(loop[j]),int(loop[k]),int(new[j]),int(new[k]);out.extend([[b,a,aa],[b,aa,bb]]);ids.append(len(out));out.append([bb,aa,tip])
                disk_faces.extend(ids);rows.append(dict(mouth=index,kind='external-fluid-classifier',extensionMeters=extension,storageGuardMeters=guard,newFaceIds=ids))
        result.append(dict(positions=np.asarray(points),faces=np.asarray(out,int),disk_faces=np.array(disk_faces,int),junction_faces=np.array(junction_faces,int),mouths=rows))
    return result[0],result[1]

def annular_mouths(paired_loops):
    """Bridge explicit corresponding outer/inner rims, never their fluid disks."""
    faces=[]
    for outer,inner in paired_loops:
        outer=np.asarray(outer,int);inner=np.asarray(inner,int)
        if outer.ndim!=1 or outer.shape!=inner.shape or len(outer)<3 or len(set(outer))!=len(outer)or len(set(inner))!=len(inner)or set(outer)&set(inner):raise ValueError('Distinct equal-count corresponding native rims required')
        for j,k in zip(range(len(outer)),np.roll(np.arange(len(outer)),-1)):
            a,b,aa,bb=int(outer[j]),int(outer[k]),int(inner[j]),int(inner[k]);faces.extend([[b,a,aa],[b,aa,bb]])
    return np.asarray(faces,int).reshape(-1,3)
