"""Transport complete source profiles with proper frames and harmonic cores.

Usage: python -m geometry_tools.commands.transport_tube_profiles FIT_DIRECTORY
       --output NEW_DIRECTORY
This constructs a review candidate. Old fixed-plane fitting certificates are
invalid after rotation; fresh surface/section/normal/packing audits are required.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from geometry_tools.course_io import sha,load_arrays,load_course_spec,write_json
from geometry_tools.algorithms.profile_chain_frames import branch_anchored_frames
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.model import Mesh
from geometry_tools.audit import topology


def profile_area(p):
    q=p-p.mean(0);return float(np.linalg.norm(np.cross(q,np.roll(q,-1,axis=0)).sum(0))*.5)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('fit_directory');parser.add_argument('--output',required=True);parser.add_argument('--core-method',choices=['harmonic','arap'],default='harmonic');args=parser.parse_args();start=time.monotonic();folder=Path(args.fit_directory).resolve();out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new immutable candidate directory')
    fit=json.loads((folder/'run.json').read_text());spec,newg,_,plan,course,active,fixed,_,_=load_course_spec(fit['specPath'])
    if fit.get('displacementsSha256')!=sha(folder/fit['displacements']):raise ValueError('Stale guide field')
    original_graph=Path(spec['source_curvature_arrays']);source_path=original_graph.with_name('source-graph.json');source=json.loads(source_path.read_text());g=load_arrays(original_graph)
    if source['graphSha256']!=sha(original_graph)or not np.array_equal(g['centres'],newg['centres'][:len(g['centres'])]):raise ValueError('Source profile identities or frame changed')
    delta=load_arrays(folder/fit['displacements'])['delta'][:len(g['centres'])];protected={k for k in fixed if k<len(g['centres'])};target=g['centres']+delta;profile=set(map(int,g['profile_cells']));adj={key:set()for key in profile};degree=np.zeros(len(g['centres']),int);definitions={r['partId']:[]for r in source['records']}
    for a,b in g['edges']:
        degree[a]+=1;degree[b]+=1
        if int(a)in profile and int(b)in profile:adj[int(a)].add(int(b));adj[int(b)].add(int(a))
    for p in source['profiles']:definitions[p['partId']].append(p)
    ends={k for k in profile if len(adj[k])!=2 or degree[k]!=2 or k in protected};seen=set();rotations=np.repeat(np.eye(3)[None],len(g['centres']),axis=0);chains=[]
    for first in sorted(ends):
        for nxt in sorted(adj[first]):
            if tuple(sorted((first,nxt)))in seen:continue
            chain=[first];previous,current=first,nxt
            while True:
                seen.add(tuple(sorted((previous,current))));chain.append(current)
                if current in ends:break
                following=next(x for x in adj[current]if x!=previous);previous,current=current,following
            if len(chain)<3:continue
            rotations[chain]=branch_anchored_frames(g['centres'][chain],target[chain],np.eye(3),np.eye(3));chains.append(chain)
    unresolved_edges=sorted({tuple(sorted((a,b)))for a in adj for b in adj[a]}-seen);out.mkdir();rows=[];maximum_pair=0.;maximum_area=0.
    for record in source['records']:
        pid=record['partId']
        if pid not in spec['moving_parts']:continue
        if sha(record['path'])!=record['sha256']:raise ValueError('Source mesh changed')
        base=load_arrays(record['path']);v=base['vertices'];f=base['faces'];off=record['offset'];cell=g['cells'][off:off+len(v)];_,representative,physical=np.unique(base['source_vertices'],axis=0,return_index=True,return_inverse=True);defs=definitions[pid];raw_profiles=[np.array(p['rawIds'])for p in defs];keys=[p['cell']for p in defs];R=rotations[keys].copy();centres=[]
        fixed_raw=np.flatnonzero(np.isin(cell,list(protected)))
        for i,(raw,key)in enumerate(zip(raw_profiles,keys)):
            centre=v[representative[physical[raw]]].mean(0)
            if key in protected:R[i]=np.eye(3);centres.append(centre)
            else:centres.append(g['centres'][key]+delta[key]+R[i]@(centre-g['centres'][key]))
        loops=boundary_loops(physical[f]);boundary=np.unique(np.concatenate([representative[loop]for loop in loops]))if loops else np.array([],int);covered=set(map(int,physical[np.concatenate(raw_profiles)]))if raw_profiles else set();prescribed=np.array([raw for raw in boundary if int(physical[raw])not in covered],int);prescribed_positions=[]
        for raw in prescribed:
            key=cell[raw]
            if key in protected:prescribed_positions.append(v[raw])
            elif np.array_equal(rotations[key],np.eye(3)):prescribed_positions.append(v[raw]+delta[key])
            else:prescribed_positions.append(g['centres'][key]+delta[key]+rotations[key]@(v[raw]-g['centres'][key]))
        result=transport_profile_cores(v,f,raw_profiles,np.array(centres),R,fixed_vertices=fixed_raw,physical_ids=physical,prescribed_vertices=prescribed,prescribed_positions=np.array(prescribed_positions).reshape(-1,3),core_method=args.core_method);q=result['vertices'];assert np.array_equal(q[fixed_raw],v[fixed_raw]);checks=[]
        for p,raw in zip(defs,raw_profiles):
            a,b=base['source_vertices'][raw],q[raw];pair=float(np.abs(np.linalg.norm(a[:,None]-a[None,:],axis=2)-np.linalg.norm(b[:,None]-b[None,:],axis=2)).max());rel=abs(profile_area(b)/profile_area(a)-1);maximum_pair=max(maximum_pair,pair);maximum_area=max(maximum_area,rel);assert pair<5e-8 and rel<.005;checks.append(dict(kind=p['kind'],cell=p['cell'],pairDistanceError_m=pair,areaRelativeError=rel))
        t0=topology(Mesh(v,f,'atlas-world-metres'));t1=topology(Mesh(q,f,'atlas-world-metres'));keys_topo=['open_edges','nonmanifold_edges','inconsistent_winding_edges','duplicate_faces','zero_area_faces','invalid_vertex_links'];same=all(t0[k]==t1[k]for k in keys_topo)and sorted(c['euler']for c in t0['components'])==sorted(c['euler']for c in t1['components']);assert same,'Transport changed source topology'
        path=out/f'part-{len(rows):03}.npz';np.savez_compressed(path,vertices=q,faces=f,source_vertices=v,source_original_vertices=base['source_vertices'],physical_ids=physical,profile_rotations=result['profile_rotations'],source_normals=base['normals'],normalSigns=base['normalSigns'],**({'uvs':base['uvs']}if'uvs'in base else{}));rows.append(dict(partId=pid,geometryFile=path.name,sha256=sha(path),sourceArraySha256=sha(record['path']),sourceOffset=off,profiles=checks,rotatedProfiles=int(np.any(abs(R-np.eye(3))>1e-12,axis=(1,2)).sum()),fixedRawPoints=len(fixed_raw),prescribedSourcePortPoints=len(prescribed),harmonicResidual_m=result['harmonic_equation_maximum_residual_m'],harmonicResidualAppliesToInitializationOnly=result['harmonic_residual_applies_to_initialization_only'],coreMethod=result['core_method'],arapEvidence=result['arap_evidence'],topologyPreserved=same))
        if time.monotonic()-start>20:raise TimeoutError('20s native proper profile/core transport deadline')
    report=dict(schema='proper-source-profile-transport/v1',parentFitReportSha256=sha(folder/'run.json'),specPath=fit['specPath'],sourceGraphSha256=sha(source_path),originalSourceGraph=str(source_path),rows=rows,frameChains=len(chains),unresolvedClosedProfileGraphEdges=unresolved_edges,maximumSourcePairDistanceError_m=maximum_pair,maximumSourceAreaRelativeError=maximum_area,allRotationsProper=True,sourceFacesAndUVsRetained=True,oldLinearSectionCertificateInvalidated=True,posedNormalsAndPackedEncodingPending=True,wholeBodyAndSelfSurfaceAuditRequired=True,sourceCodeSha256=sha(__file__),installed=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start);write_json(out/'report.json',report);print(json.dumps({k:v for k,v in report.items()if k not in ['rows','unresolvedClosedProfileGraphEdges']}))


if __name__=='__main__':main()
