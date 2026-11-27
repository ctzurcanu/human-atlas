"""Read-only full-source AP corridor witnesses for the single ITA v8 trial."""
from pathlib import Path
import json, sys, time
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file
from geometry_tools.model import AtlasReader

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1'
OUT=BASE/'ita-candidate-ap-corridor-independent-v1'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME='model:local-ta98:assembled-metres'

def runs(ids):
    result=[]
    for i in ids:
        if not result or i>result[-1][-1]+1: result.append([i])
        else: result[-1].append(i)
    return [[r[0],r[-1]] for r in result]

class AP:
    def __init__(self, mesh):
        t=mesh.triangles
        self.tri=t; self.lo=t[:,:,:2].min(1); self.hi=t[:,:,:2].max(1)
        self.a=t[:,0]; self.e=t[:,1]-t[:,0]; self.f=t[:,2]-t[:,0]
        self.det=self.e[:,0]*self.f[:,1]-self.e[:,1]*self.f[:,0]
    def hits(self, p):
        ids=np.flatnonzero(np.all(self.lo<=p[:2]+1e-12,axis=1)&np.all(self.hi>=p[:2]-1e-12,axis=1)&(np.abs(self.det)>1e-16))
        if not len(ids): return []
        q=p[:2]-self.a[ids,:2];e=self.e[ids];f=self.f[ids];det=self.det[ids]
        u=(q[:,0]*f[:,1]-q[:,1]*f[:,0])/det
        v=(e[:,0]*q[:,1]-e[:,1]*q[:,0])/det
        take=(u>=-1e-10)&(v>=-1e-10)&(u+v<=1+1e-10)
        ids=ids[take];u=u[take];v=v[take];z=self.a[ids,2]+u*self.e[ids,2]+v*self.f[ids,2]
        return [dict(sourceFace=int(i),point_m=[float(p[0]),float(p[1]),float(zz)],barycentric=[float(1-uu-vv),float(uu),float(vv)]) for i,zz,uu,vv in zip(ids,z,u,v)]

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False)
    src=json.loads((BASE/'ita-distal-closure-context-v6/evidence.json').read_bytes())
    sha=src['sourceCatalogueSha256']; assert digest_file(CAT)==sha
    atlas=json.loads(CAT.read_bytes()); parts={p['id']:p for p in atlas['parts']}; reader=AtlasReader()
    context={r['partId']:r for r in src['rows']}; bindings={}; ap={}
    ids=['ZA:Pleura','ZA:Anterior segmental artery of left lung','ZA:Anterior segmental artery of right lung']
    for s in 'lr':
        ids += ['ZA:Transversus thoracis muscle.'+s,'ZA:Internal intercostal muscles.'+s,'ZA:Innermost intercostal muscles.'+s]
        ids += ['ZA:Costal cartilage of '+rib+' rib.'+s for rib in ['first','second','third','fourth','fifth','sixth','seventh']]
    for pid in ids:
        m=reader.read_part(atlas,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha)
        ap[pid]=AP(m);clipped=context[pid]
        bindings[pid]=dict(wholeSourceGeometrySha256=m.geometry_sha,wholeSourceTriangles=len(m.faces),v6SnapshotTriangles=clipped['retainedTriangles'],v6SnapshotWasClipped=clipped['retainedTriangles']!=len(m.faces),fullOriginalSourceUsed=True)
    result=[]
    for side in 'lr':
        file=BASE/'one-parasternal-family-candidate-v8'/f'{side}-ita-field.npz'
        z=np.load(file); centers=z['target_centers'];rings=z['encoded_sections']; stationRows=[]
        anterior=['ZA:Costal cartilage of '+rib+' rib.'+side for rib in ['first','second','third','fourth','fifth','sixth','seventh']]
        back=['ZA:Transversus thoracis muscle.'+side,'ZA:Pleura']
        others=['ZA:Internal intercostal muscles.'+side,'ZA:Innermost intercostal muscles.'+side,'ZA:Anterior segmental artery of left lung'if side=='l'else'ZA:Anterior segmental artery of right lung']
        for i,(c,ring) in enumerate(zip(centers,rings)):
            queries=[]
            for qi,p in enumerate(np.vstack([c,ring])):
                rows=[]
                for pid in anterior+back+others:
                    hits=ap[pid].hits(p)
                    if hits: rows.append(dict(partId=pid,apHits=hits,minZ_m=min(x['point_m'][2]for x in hits),maxZ_m=max(x['point_m'][2]for x in hits)))
                cartilage=[(r['minZ_m'],r)for r in rows if r['partId']in anterior]
                posterior=[(r['maxZ_m'],r)for r in rows if r['partId']in back]
                upper=min(cartilage,key=lambda t:t[0])if cartilage else None
                lower=max(posterior,key=lambda t:t[0])if posterior else None
                gap=upper[0]-lower[0]if upper and lower else None
                queries.append(dict(query='center'if qi==0 else'actualRingVertex',ringVertex=qi-1 if qi else None,point_m=p.tolist(),sources=rows,
                    cartilagePosteriorZ_m=upper[0]if upper else None,posteriorTissueAnteriorZ_m=lower[0]if lower else None,
                    cartilageSource=upper[1]['partId']if upper else None,posteriorTissueSource=lower[1]['partId']if lower else None,
                    twoSidedAPGap_m=gap,lessThan2_6mm=gap is not None and gap<.0026,
                    candidatePosteriorToPosteriorTissueAnteriorSurface=bool(lower and p[2]<lower[0]),
                    candidateAnteriorToCartilagePosteriorSurface=bool(upper and p[2]>upper[0])))
            stationRows.append(dict(station=i,center_m=c.tolist(),queries=queries))
        result.append(dict(side=side,candidateField=str(file),candidateFieldSha256=digest_file(file),stationRows=stationRows,
            summary=dict(centerTwoSidedStations=sum(r['queries'][0]['twoSidedAPGap_m']is not None for r in stationRows),
                centerGapBelow2_6mmIntervals=runs([r['station']for r in stationRows if r['queries'][0]['lessThan2_6mm']]),
                anyRingRayGapBelow2_6mmIntervals=runs([r['station']for r in stationRows if any(x['lessThan2_6mm']for x in r['queries'][1:])]),
                centerPosteriorToPosteriorSurfaceIntervals=runs([r['station']for r in stationRows if r['queries'][0]['candidatePosteriorToPosteriorTissueAnteriorSurface']]))))
    atomic_json(OUT/'evidence.json',dict(schema='ita-candidate-full-source-ap-corridor-independent/v1',sourceCatalogueSha256=sha,
        noGeometryChanged=True,fullSourceBindings=bindings,sides=result,anatomicalAcceptance=False,
        limitation='These are AP literal triangle witnesses, not normal thickness, a closed solid or a full 3D corridor proof. Source tissue envelopes may overlap. Missing ray hits do not establish room. Ring-vertex and center rays do not certify faces between samples.',
        elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),summary=[dict(side=r['side'],**r['summary'])for r in result],elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__': main()
