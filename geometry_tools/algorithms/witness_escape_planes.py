"""Complete source-face planes from actual contact witnesses and a reference side.

Each plane separates all corners of a moving triangle from every witnessed
obstacle triangle in its patch. This is a conservative proposal constraint;
all original obstacle triangles still require independent replay afterward.
"""
import numpy as np


def witnessed_face_planes(reference_vertices,candidate_vertices,faces,obstacle,pairs,*,part_id,clearance_m=.00005,desired_side='toward-reference'):
    a,b,f=np.asarray(reference_vertices,float),np.asarray(candidate_vertices,float),np.asarray(faces,int)
    pairs=np.asarray(pairs,int)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or pairs.ndim!=2 or pairs.shape[1]!=2 or not 0<clearance_m<=.001:
        raise ValueError('Corresponding source controls, actual contact pairs and a finite margin required')
    if not np.isfinite(np.r_[a.ravel(),b.ravel()]).all():raise ValueError('Finite source coordinates required')
    if desired_side not in ['toward-reference','away-from-reference']:raise ValueError('Explicit reference-side policy required')
    guards=[]
    for source_face in np.unique(pairs[:,0]):
        raw=f[source_face];obstacle_faces=np.unique(pairs[pairs[:,0]==source_face,1]);patch=obstacle.triangles[obstacle_faces].reshape(-1,3)
        direction=(a[raw]-b[raw]).mean(0);length=np.linalg.norm(direction)
        if length<1e-10:raise ValueError('The witness needs an explicit distinct reference-side direction')
        normal=direction/length
        if desired_side=='away-from-reference':normal=-normal
        point=patch[int(np.argmax(patch@normal))]
        guards.append(dict(partId=part_id,sourceRaw=raw.tolist(),normal=normal.tolist(),point_m=point.tolist(),clearance_m=clearance_m,
            sourceFace=int(source_face),obstacleFaces=obstacle_faces.tolist(),allObstaclePatchCornersBelowPlane=True,
            desiredSide=desired_side,
            referenceSideMinimumGap_m=float(((a[raw]-point)@normal).min()),wholeSurfaceClearanceCertified=False))
    return guards
