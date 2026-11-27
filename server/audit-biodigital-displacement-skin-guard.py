"""Bounded moved-vertex skin guard; unproved native distances never prune.

Distance lemma: if corresponding triangle vertices differ by at most epsilon,
their surfaces have Hausdorff distance <= epsilon (same barycentric weights).
Distance to a set is 1-Lipschitz in both its query and Hausdorff perturbations.
Thus d_true >= max(0, d_native_computed - epsilon - query_error - eta), ONLY
when eta certifies native distance's numerical/global-minimum error. Open3D's
installed Float32 API supplies no eta; this guard uses lower=0 and all rays.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.oriented_rays import OrientedRayIndex
from geometry_tools.open3d_backend import library

def conditional_distance_lower_bound(native_distance,vertex_error,query_error,certified_native_error=None):
    d,q=np.broadcast_arrays(np.asarray(native_distance,float),np.asarray(query_error,float));e=float(vertex_error)
    if not np.isfinite(np.r_[d.ravel(),q.ravel(),e]).all()or np.any(d<0)or np.any(q<0)or e<0:raise ValueError('Nonnegative finite distance/error inputs required')
    if certified_native_error is None:return np.zeros_like(d),False
    eta=np.asarray(certified_native_error,float)
    if not np.isfinite(eta).all()or np.any(eta<0):raise ValueError('Certified native error must be finite and nonnegative')
    return np.maximum(0,d-e-q-eta),True

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb')as stream:
        while block:=stream.read(1048576):h.update(block)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--candidate-stage',required=True);ap.add_argument('--baseline-stage',required=True)
    args=ap.parse_args();assert all(Path(v).name==v and v not in ['.','..']for v in vars(args).values())
    began=time.monotonic();work=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
    candidate=work/args.candidate_stage/'candidate-v1';baseline=work/args.baseline_stage/'candidate-v1'
    out=work/args.candidate_stage/'localized-skin-displacement-guard-v1';assert not out.exists();out.mkdir()
    report=dict(schema='batched-moved-vertex-displacement-skin-guard/v1',status='preparing',complete=False,
        candidateStage=args.candidate_stage,baselineStage=args.baseline_stage,geometryModified=False,installed=False,anatomicalAcceptance=False)
    def save():
        report['elapsedSeconds']=time.monotonic()-began;(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    try:
        terminal=json.loads((candidate/'checkpoint.json').read_text());assert terminal['status']=='candidate-numerical-gates-passed-anatomy-unverified'
        opaths=[p/'instance-ownership.json'for p in [baseline,candidate]];ownership=[json.loads(p.read_text())for p in opaths];assert ownership[0]==ownership[1]
        records={r['instanceIndex']:r for r in ownership[0]['instances']};inputs={str(p):sha(p)for p in opaths}
        arrays={k:[]for k in ['before','after','owner','vertex','mesh','packed']};parts=[];frames=set()
        for mid,name in enumerate(['arterial','venous','arterial-remnant-companion']):
            paths=[p/(name+'.npz')for p in [baseline,candidate]];aps=[p/(name+'-ancestry.npz')for p in [baseline,candidate]]
            for p in paths+aps:inputs[str(p)]=sha(p)
            with np.load(paths[0],allow_pickle=False)as b,np.load(paths[1],allow_pickle=False)as a,np.load(aps[0],allow_pickle=False)as ba,np.load(aps[1],allow_pickle=False)as aa:
                bv,av=b['vertices'],a['vertices'];assert bv.shape==av.shape and np.array_equal(b['faces'],a['faces'])
                for key in ba.files:assert key in aa.files and np.array_equal(ba[key],aa[key]),'Source ancestry differs: '+key
                frames.update([str(b['frame']),str(a['frame'])]);ids=np.flatnonzero(np.any(av!=bv,axis=1))
                arrays['before'].append(bv[ids]);arrays['after'].append(av[ids]);arrays['owner'].append(ba['candidate_vertex_source_instance'][ids]);arrays['vertex'].append(ba['candidate_vertex_source_vertex'][ids]);arrays['mesh'].append(np.full(len(ids),mid,np.int8));arrays['packed'].append(ids)
                parts.append(dict(mesh=name,literalChangedVertices=len(ids),threshold_m=0.))
        assert len(frames)==1
        data={k:np.concatenate(v)for k,v in arrays.items()};n=len(data['before']);assert n>0
        np.savez_compressed(out/'moved-source-controls.npz',**data)
        benchmarkpath=work/'body-envelope-fit-v2/current-envelope-ray-benchmark-v1/report.json';benchmark=json.loads(benchmarkpath.read_text())
        skinpath=Path(benchmark['skin_path']);bp=Path(benchmark['boundary_evidence']['path']);br=bp.parent/'report.json';boundaryreport=json.loads(br.read_text())
        assert sha(skinpath)==benchmark['skin_sha256'];assert sha(bp)==benchmark['boundary_evidence']['sha256'];assert sha(br)==benchmark['boundary_evidence']['report_sha256']
        for p in [benchmarkpath,skinpath,bp,br]:inputs[str(p)]=sha(p)
        boundary=np.load(bp,allow_pickle=False);assert boundary['primary_segments_m'].shape==(760,2,3)
        index=OrientedRayIndex(Mesh.load(skinpath),nthreads=2,uncertainty_segments=boundary['primary_segments_m'],
            uncertainty_component_ids=boundary['primary_boundary_component_index'],boundary_tolerance_m=boundaryreport['geometricTolerance_m'],boundary_evidence=benchmark['boundary_evidence'])
        directions=np.asarray(benchmark['directions']);assert directions.shape==(6,3)
        report.update(status='classifying-all-moved-vertices',candidateTerminal=terminal,inputBindings=inputs,
            exactSourceAncestryEqual=True,movedVertices=n,parts=parts,sourceFrame=next(iter(frames)),skinFrame=index.mesh.frame,
            frameContract=benchmark['frame_contract'],sourceIdentityEncoding='source instance index + original source vertex id; exact lookup in bound instance-ownership.json',
            primaryBoundarySegments=760,maxClassifierQueriesPerBatch=8192,
            broadphase=dict(formula='max(0,dNative-epsilon-queryRounding-etaNative)',mathematicalLemma='Triangle barycentric correspondence gives surface Hausdorff bound epsilon; set-distance Lipschitz inequality gives the conditional formula.',
                epsilonSkinVertexRounding_m=index.index.vertex_rounding_m,certifiedNativeErrorAvailable=False,
                nativeAPI='Open3D '+library().__version__+' Float32 compute_distance/compute_closest_points; installed documentation provides no numerical/global-completeness error certificate.',
                wrapperEvidence='geometry_tools.open3d_backend.SurfaceIndex.closest_points(certify=False) returns literal-triangle witness upper bounds only.',
                certifiedLowerBoundUsed_m=0.,prunedVertices=0,unprovedPointsPruned=False),
            processedMovedVertices=0,remainingMovedVertices=n)
        save();labels={-1:'uncertain',0:'outside',1:'inside'};transitions={labels[b]+'_to_'+labels[a]:0 for b in [-1,0,1]for a in [-1,0,1]}
        beforecounts={v:0 for v in labels.values()};aftercounts=beforecounts.copy();perowner={};batches=[];examples=[];pending=0;processed=0;firstproof=None
        report.update(before=beforecounts,after=aftercounts,transitions=transitions)
        for begin in range(0,n,4096):
            # Self-limit leaves time to retain an honest incomplete receipt before the outer 40 s supervisor.
            if time.monotonic()-began>32:break
            end=min(n,begin+4096);count=end-begin;b=data['before'][begin:end];a=data['after'][begin:end]
            local=((b-index.index.origin)/index.index.scale).astype(np.float32)
            queryround=np.linalg.norm(local.astype(float)*index.index.scale+index.index.origin-b,axis=1)
            native=index.index.scene.compute_distance(library().core.Tensor(local),nthreads=2).numpy().astype(float)*index.index.scale
            lower,certified=conditional_distance_lower_bound(native,index.index.vertex_rounding_m,queryround,None)
            assert not certified and not np.any(lower>0)
            query=np.r_[b,a];classified,proof=index.classify(query,directions,near_surface_m=benchmark['near_surface_m'],min_reliable_rays=benchmark['min_reliable_rays'])
            cb,ca=classified['classification'][:count],classified['classification'][count:]
            if firstproof is None:firstproof=proof
            for x in [-1,0,1]:
                beforecounts[labels[x]]+=int((cb==x).sum());aftercounts[labels[x]]+=int((ca==x).sum())
                for y in [-1,0,1]:transitions[labels[x]+'_to_'+labels[y]]+=int(((cb==x)&(ca==y)).sum())
            pending+=int(((cb==-1)|(ca==-1)).sum());owner=data['owner'][begin:end]
            for o in np.unique(owner):
                mask=owner==o;r=perowner.setdefault(int(o),dict(sourceId=records[int(o)]['sourceInstanceId'],name=records[int(o)]['name'],movedVertices=0,
                    before={v:0 for v in labels.values()},after={v:0 for v in labels.values()},newInsideToOutside=0,pendingAnyUncertainty=0))
                r['movedVertices']+=int(mask.sum());r['newInsideToOutside']+=int((mask&(cb==1)&(ca==0)).sum());r['pendingAnyUncertainty']+=int((mask&((cb==-1)|(ca==-1))).sum())
                for label,value in labels.items():r['before'][value]+=int((mask&(cb==label)).sum());r['after'][value]+=int((mask&(ca==label)).sum())
            for j in np.flatnonzero((cb==1)&(ca==0))[:max(0,24-len(examples))]:
                row=records[int(owner[j])];examples.append(dict(movedIndex=begin+int(j),sourceId=row['sourceInstanceId'],name=row['name'],originalSourceVertex=int(data['vertex'][begin+j]),
                    candidateMeshIndex=int(data['mesh'][begin+j]),candidatePackedVertex=int(data['packed'][begin+j]),before_m=b[j].tolist(),after_m=a[j].tolist(),
                    afterNearestSkinDistanceUpperBound_m=float(classified['surface_distance_upper_bound_m'][count+j])))
            batchpath=out/f'batch-{len(batches):04d}.npz';np.savez_compressed(batchpath,movedIndex=np.arange(begin,end),nativeFloat32Distance_m=native,
                queryRounding_m=queryround,certifiedDistanceLowerBound_m=lower,**classified)
            batches.append(dict(path=str(batchpath),begin=begin,end=end,vertices=count,classifierQueries=2*count,
                nativeEvents=proof['native_proposed_events'],elapsedSeconds=proof['elapsed_seconds']))
            processed=end;report.update(processedMovedVertices=processed,remainingMovedVertices=n-processed,transitions=transitions,before=beforecounts,after=aftercounts,pendingAnyUncertainty=pending,batches=batches)
            save()
        complete=processed==n;newexits=transitions['inside_to_outside'];report.update(status='complete-no-new-clean-skin-exits'if complete and not newexits else'complete-skin-exits-found'if complete else'incomplete-wall-budget-stop',
            complete=complete,processedMovedVertices=processed,remainingMovedVertices=n-processed,newCleanInsideToOutside=newexits,
            pendingAnyUncertainty=pending,perVessel=sorted(perowner.values(),key=lambda r:(-r['newInsideToOutside'],-r['pendingAnyUncertainty'],r['name'])),newExitExamples=examples,
            rayPolicyProof=firstproof,guardPassed=complete and newexits==0 and pending==0,
            qualification='All literal changed blood vertices, original source identities retained. Native lower-bound pruning is unavailable and disabled. Six oriented rays replay native events on literal Float64 triangles with bound primary skin uncertainty. No nearest-normal sign test. This vertex gate does not certify moved triangle containment, lumen or anatomical packing.',
            scope='Blood candidate only; independently stored 21-part heart/cusp banks are not part of these source ownership arrays.',scriptSha256=sha(__file__))
        save();print(json.dumps({k:report[k]for k in ['status','complete','movedVertices','processedMovedVertices','remainingMovedVertices','before','after','transitions','newCleanInsideToOutside','pendingAnyUncertainty','guardPassed','elapsedSeconds']}))
    except Exception as exc:
        report.update(status='failed-retained-for-diagnosis',error=repr(exc),traceback=traceback.format_exc());save();raise

if __name__=='__main__':main()
