"""Area-preserving course shear shared with the attached gastric cardia."""
import numpy as np

def transform(points,displacement,*,sex):
    p=np.asarray(points,float);knots=[1.242,1.272,1.303,1.355]if sex=='male'else[1.202,1.231,1.252,1.290]
    t=np.clip((p[:,1]-knots[0])/(knots[1]-knots[0]),0,1);u=np.clip((p[:,1]-knots[2])/(knots[3]-knots[2]),0,1)
    rise=t*t*(3-2*t);fall=1-u*u*(3-2*u);weight=rise*fall
    derivative=6*t*(1-t)/(knots[1]-knots[0])*fall-rise*6*u*(1-u)/(knots[3]-knots[2])
    direction=np.asarray(displacement,float);assert direction.shape==(3,)and direction[1]==0
    q=p+weight[:,None]*direction;J=np.broadcast_to(np.eye(3),(len(p),3,3)).copy();J[:,:,1]+=derivative[:,None]*direction
    return q,J
