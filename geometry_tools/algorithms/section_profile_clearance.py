"""Rigid anterior bypass proposals against literal arterial section contours.

Use only for an explicitly identified anterior branch. Sampled profile-edge
clearance is a proposal check, not continuous3D surface or lumen acceptance.
"""
import numpy as np
import shapely
from .literal_plane_sections import literal_plane_sections

class AnteriorSectionProfileClearance:
 def __init__(self,vertices,faces,*,margin_m=.00005,front_padding_m=.000075,sample_spacing_m=.0001,maximum_translation_m=.0015):
  self.vertices=np.asarray(vertices,float);self.faces=np.asarray(faces,int);self.margin=float(margin_m);self.padding=float(front_padding_m);self.spacing=float(sample_spacing_m);self.maximum=float(maximum_translation_m);self.cache={}
  if not 0<self.margin<=self.padding or self.spacing<=0 or self.maximum<=0 or not np.isfinite([self.margin,self.padding,self.spacing,self.maximum]).all():raise ValueError('Finite positive profile-clearance limits required')
 def section(self,y):
  key=round(float(y),12)
  if key not in self.cache:
   cut=literal_plane_sections(self.vertices,self.faces,[0,y,0],[0,1,0])
   if len(cut['loops'])!=1 or cut['open_components']or cut['unresolved_source_faces']:raise ValueError('Arterial section must be a single actual resolved closed contour')
   loop=cut['loops'][0][:,[0,2]];poly=shapely.polygons(loop)
   if not shapely.is_valid(poly):raise ValueError('Invalid actual arterial contour; no implicit repair permitted')
   self.cache[key]=(loop,poly)
  return self.cache[key]
 def propose(self,profile):
  profile=np.asarray(profile,float)
  if profile.ndim!=2 or profile.shape[1]!=3 or len(profile)<3 or not np.isfinite(profile).all():raise ValueError('A complete finite posed source profile is required')
  samples=[]
  for a,b in zip(profile,np.roll(profile,-1,axis=0)):
   t=np.linspace(0,1,max(2,int(np.ceil(np.linalg.norm(b-a)/self.spacing))),endpoint=False);samples.extend(a+t[:,None]*(b-a))
  need=0.;old=[]
  for point in samples:
   loop,poly=self.section(point[1]);probe=shapely.points(point[[0,2]]);distance=float(shapely.distance(probe,shapely.boundary(poly)));old.append(-distance if shapely.contains(poly,probe)else distance);xmin,_,xmax,_=shapely.bounds(poly)
   if point[0]<xmin-self.padding or point[0]>xmax+self.padding:continue
   x=np.clip(point[0],xmin,xmax);a=loop;b=np.roll(loop,-1,axis=0);dx=b[:,0]-a[:,0];valid=abs(dx)>1e-14;t=(x-a[valid,0])/dx[valid];hit=(t>=-1e-9)&(t<=1+1e-9);z=a[valid,1]+t*(b[valid,1]-a[valid,1])
   if not hit.any():raise ValueError('Actual anterior contour ray is unresolved')
   need=max(need,float(z[hit].max()+self.padding-point[2]))
  if not 0<=need<=self.maximum:raise ValueError('Local branch correction exceeds the declared rigid translation bound')
  gaps=[]
  for point in samples:
   _,poly=self.section(point[1]);probe=shapely.points(point[[0,2]]+[0,need]);distance=float(shapely.distance(probe,shapely.boundary(poly)))
   if shapely.contains(poly,probe)or distance<self.margin:raise ValueError('Proposed profile still fails sampled arterial clearance')
   gaps.append(distance)
  return np.array([0.,0.,need]),dict(edgeSamples=len(samples),maximumSampleSpacing_m=self.spacing,oldMinimumSignedGap_m=min(old),proposedMinimumGap_m=min(gaps),properTranslationOnly=True,continuousSurfaceAcceptance=False)
