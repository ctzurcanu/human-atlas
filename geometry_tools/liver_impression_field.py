"""Positive smooth hepatic packing field; explicit source-coordinate support."""
import numpy as np

def smooth(x,lo,hi):
    t=np.clip((x-lo)/(hi-lo),0,1)
    return t*t*(3-2*t),6*t*(1-t)/(hi-lo)

def transform(points,amplitude,*,lower=1.258,upper=1.278):
    p=np.asarray(points,float);wx,dx=smooth(p[:,0],-.05,.03);wy,dy=smooth(p[:,1],lower,upper);sz,dz=smooth(p[:,2],.025,.065);wz=1-sz
    q=p.copy();q[:,0]-=amplitude*wx*wy*wz
    J=np.broadcast_to(np.eye(3),(len(p),3,3)).copy();J[:,0,0]-=amplitude*dx*wy*wz;J[:,0,1]=-amplitude*wx*dy*wz;J[:,0,2]=amplitude*wx*wy*dz
    return q,J

def posterior_impression(points,amplitude,*,down=0.):
    """Indent the posterior left hepatic surface without shifting its upper dome.

    The native hilar annuli below y=1.255 remain exactly fixed. Coordinates are
    in the original male liver frame. Its broad differentiable support avoids
    moving individual surface vertices independently of the owned tissue.
    """
    p=np.asarray(points,float)
    xl,dxl=smooth(p[:,0],-.055,-.025)
    xr,dxr=smooth(p[:,0],.025,.070)
    yl,dyl=smooth(p[:,1],1.257,1.268)
    yr,dyr=smooth(p[:,1],1.320,1.360)
    zz,dzz=smooth(p[:,2],.020,.070)
    wx=xl*(1-xr); dx=dxl*(1-xr)-xl*dxr
    wy=yl*(1-yr); dy=dyl*(1-yr)-yl*dyr
    wz=1-zz
    weight=wx*wy*wz
    grad=np.c_[dx*wy*wz,wx*dy*wz,-wx*wy*dzz]
    direction=np.array([0.,-down,amplitude])
    q=p+weight[:,None]*direction
    J=np.broadcast_to(np.eye(3),(len(p),3,3)).copy()
    J+=direction[None,:,None]*grad[:,None,:]
    return q,J

def slab_impression(points,amplitude,*,down=0.):
    """Piecewise affine common field with exact, conformable y-plane support."""
    p=np.asarray(points,float);y=p[:,1]
    knots=np.array([1.257,1.277,1.320,1.350]);values=np.array([0.,1.,1.,0.])
    w=np.interp(y,knots,values,left=0.,right=0.)
    dw=np.where((y>knots[0])&(y<knots[1]),1/(knots[1]-knots[0]),np.where((y>knots[2])&(y<knots[3]),-1/(knots[3]-knots[2]),0.))
    direction=np.array([0.,-down,amplitude]);q=p+w[:,None]*direction
    J=np.broadcast_to(np.eye(3),(len(p),3,3)).copy();J[:,:,1]+=dw[:,None]*direction
    return q,J
