"""Smooth size correction that leaves the cardia and pyloric mouths fixed."""
import numpy as np

def expand_body(points,strength,*,anchor=.030,transition=.010):
    p=np.asarray(points,float);u=p[:,0]-anchor;t=np.clip(u/transition,0,1)
    h=np.where(u<=0,0,np.where(u<transition,transition*(t**3-.5*t**4),u-.5*transition))
    derivative=t*t*(3-2*t)
    q=p.copy();q[:,0]+=strength*h
    J=np.broadcast_to(np.eye(3),(len(p),3,3)).copy();J[:,0,0]+=strength*derivative
    return q,J
