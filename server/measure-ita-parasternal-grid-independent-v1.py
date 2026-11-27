"""Read-only sparse AP corridor grid at literal parasternal source slices."""
from pathlib import Path
import json, sys, time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file
from geometry_tools.model import AtlasReader
from importlib.util import spec_from_file_location, module_from_spec
spec=spec_from_file_location('ita_ap',Path(__file__).with_name('measure-ita-candidate-ap-corridor-independent-v1.py'))
mod=module_from_spec(spec);spec.loader.exec_module(mod)
BASE=mod.BASE;CAT=mod.CAT;FRAME=mod.FRAME

def edge_at_y(tri,y,sign):
    points=[]
    for aa,bb in [(0,1),(1,2),(2,0)]:
        a,b=tri[:,aa],tri[:,bb];dy=b[:,1]-a[:,1]
        ids=np.flatnonzero(((a[:,1]-y)*(b[:,1]-y)<=0)&(np.abs(dy)>1e-15))
        t=(y-a[ids,1])/dy[ids];q=a[ids]+t[:,None]*(b[ids]-a[ids])
        points += [dict(sourceFace=int(i),edge=[aa,bb],edgeParameter=float(tt),point_m=pt.tolist())for i,tt,pt in zip(ids,t,q)]
    return max(points,key=lambda q:sign*q['point_m'][0])if points else None

def main():
    began=time.monotonic();side=sys.argv[1];assert side in 'lr';sign=1 if side=='l'else-1
    out=BASE/('ita-parasternal-grid-independent-v1-'+side);out.mkdir(exist_ok=False)
    ctx=json.loads((BASE/'ita-distal-closure-context-v6/evidence.json').read_bytes()); sha=ctx['sourceCatalogueSha256'];assert digest_file(CAT)==sha
    atlas=json.loads(CAT.read_bytes());parts={p['id']:p for p in atlas['parts']}; reader=AtlasReader();cr={r['partId']:r for r in ctx['rows']}
    anterior=['ZA:Costal cartilage of '+s+' rib.'+side for s in ['first','second','third','fourth','fifth','sixth','seventh']]+['ZA:Internal intercostal muscles.'+side]
    posterior=['ZA:Pleura']+[p for p in cr if p.startswith('ZA:Transversus thoracis muscle.')and p.rsplit('.',1)[1].endswith(side)]
    bones=['ZA:Body of sternum:Bone-7','ZA:Manubrium of sternum:Bone-6']
    meshes={}; ap={};bindings={}
    for pid in anterior+posterior+bones:
        m=reader.read_part(atlas,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha);meshes[pid]=m
        ap[pid]=mod.AP(m);bindings[pid]=dict(geometrySha256=m.geometry_sha,wholeTriangles=len(m.faces),v6Triangles=cr[pid]['retainedTriangles'],fullSourceUsed=True)
    srcfile=BASE/'one-parasternal-family-candidate-v8'/f'{side}-ita-field.npz';z=np.load(srcfile);centers=z['target_centers'];rows=[]
    for i,c in enumerate(centers):
        edge=[]
        for pid in bones:
            hit=edge_at_y(meshes[pid].triangles,c[1],sign)
            if hit:edge.append(dict(partId=pid,**hit))
        if not edge:rows.append(dict(station=i,noActualSternalSlice=True));continue
        b=max(edge,key=lambda h:sign*h['point_m'][0]);station=dict(station=i,sourceY_m=float(c[1]),sternalLateralWitness=b,grid=[])
        for offset in [.005,.010,.015,.020]:
            x=b['point_m'][0]+sign*offset;queries=[]
            for dx,dy in [(0,0),(-.00134,-.00134),(-.00134,.00134),(.00134,-.00134),(.00134,.00134)]:
                point=np.array([x+dx,c[1]+dy,c[2]]);a=[];p=[]
                for pid in anterior+posterior:
                    hits=ap[pid].hits(point)
                    if not hits:continue
                    witness=min(hits,key=lambda h:h['point_m'][2])if pid in anterior else max(hits,key=lambda h:h['point_m'][2])
                    q=dict(partId=pid,**witness)
                    (a if pid in anterior else p).append(q)
                upper=min(a,key=lambda h:h['point_m'][2])if a else None
                lower=max(p,key=lambda h:h['point_m'][2])if p else None
                gap=upper['point_m'][2]-lower['point_m'][2]if upper and lower else None
                queries.append(dict(queryXY_m=point[:2].tolist(),upperAnteriorTissuePosteriorWitness=upper,lowerPosteriorTissueAnteriorWitness=lower,gap_m=gap,
                    allAnteriorSourceWitnesses=a,allPosteriorSourceWitnesses=p,unknown=gap is None))
            known=[q for q in queries if q['gap_m']is not None];center=queries[0]
            lo=max(q['lowerPosteriorTissueAnteriorWitness']['point_m'][2]for q in known)if len(known)==5 else None
            hi=min(q['upperAnteriorTissuePosteriorWitness']['point_m'][2]for q in known)if len(known)==5 else None
            station['grid'].append(dict(lateralFromActualSternalEdge_mm=offset*1000,centerXY_m=[float(x),float(c[1])],queries=queries,
                allFiveRaysKnown=len(known)==5,centerGap_m=center['gap_m'],centerFits2_46mm=center['gap_m']is not None and center['gap_m']>=.00246,
                centerFits2_68mm=center['gap_m']is not None and center['gap_m']>=.00268,
                fiveRayEnvelopeGap_m=hi-lo if hi is not None else None,
                fiveRayConservative2_68mmCenterInterval_m=[lo+.00134,hi-.00134]if hi is not None and hi-lo>=.00268 else None))
        rows.append(station)
    summary=[]
    for ix,offset in enumerate([5,10,15,20]):
        valid=[r for r in rows if 'grid'in r];known=[r for r in valid if r['grid'][ix]['centerGap_m']is not None]
        summary.append(dict(lateralOffset_mm=offset,knownCenterStations=len(known),centerGapAtLeast2_68mmIntervals=mod.runs([r['station']for r in known if r['grid'][ix]['centerFits2_68mm']]),
            fixedSourceEnvelopeOverlapIntervals=mod.runs([r['station']for r in known if r['grid'][ix]['centerGap_m']<0]),
            conservativeFiveRaySphereBoxFeasibleIntervals=mod.runs([r['station']for r in valid if r['grid'][ix]['fiveRayConservative2_68mmCenterInterval_m']is not None])))
    atomic_json(out/'evidence.json',dict(schema='ita-sparse-lateral-grid-ap-corridor/v1',side=side,sourceCatalogueSha256=sha,candidateFieldSha256=digest_file(srcfile),
        fullSourceBindings=bindings,anteriorBoundaryParts=anterior,posteriorBoundaryParts=posterior,rows=rows,summary=summary,
        noGeometryChanged=True,anatomicalAcceptance=False,
        limitation='Literal AP ray samples only. Missing hits remain unknown. Five rays and a 1.34mm sphere-box guard do not cover the entire swept footprint or prove a coherent continuous 3D curve. Negative gaps establish intersecting AP tissue envelopes at that ray, not clinical tissue volumes.',elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(evidence=str(out/'evidence.json'),summary=summary,elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
