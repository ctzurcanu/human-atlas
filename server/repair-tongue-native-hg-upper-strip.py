"""Complete fourteen-way triangulation of the one six-control HG upper strip.

Actual positions, strip boundary, bony origin and every other physical face
stay literal. All alternatives are tested against themselves and the complete
static muscle. Only this local topological defect is repaired; no field search.
"""
from __future__ import annotations
import argparse,importlib.util,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit,topology
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature
from tongue_interface_triangulation import loop_of,triangulations

W=ROOT/'.local-models/geometry-workbench/tongue-completion';ORIGIN=np.array([0.,1.625,0.])
def require(x,message):
    if not x:raise ValueError(message)
def read(path):
    import json
    return json.loads(Path(path).read_text())
def clear(s):return s['complete_within_pair_budget']and not s['truncated']and s['crossing_face_pairs']==s['coplanar_overlap_face_pairs']==0
def module(name,file):
    s=importlib.util.spec_from_file_location(name,ROOT/'server'/file);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main(source,output):
    began=time.monotonic();end=began+48
    def check():
        if time.monotonic()>end:raise TimeoutError('One local HG14-way repair reached its finite deadline')
    require(source.is_relative_to(W)and output.is_relative_to(W)and not output.exists(),'New local output only')
    report=read(source/'candidate.json');require(report['passed']is False,'Only the explicit failed short-course proposal is supported')
    native=Mesh.load(report['candidate']);v,first,alias=np.unique(native.vertices,axis=0,return_index=True,return_inverse=True);f=alias[native.faces]
    arrays,own=pair_contacts(native,native,0.,selfcheck=True,max_pairs=250000,max_tests=2000000);bad=np.isin(arrays['relation'],['crossing','coplanar-overlap']);pairs=arrays['pairs'][bad];ids=np.unique(pairs)
    require(len(pairs)==3 and len(ids)==4,'Only the witnessed three-pair/four-face domain may be repaired')
    patch=f[ids];loop=loop_of(patch);options=triangulations(loop);require(len(loop)==6 and len(options)==14,'A simple six-control upper strip with fourteen triangulations is required')
    static=f[~np.isin(np.arange(len(f)),ids)];histories=[];possible=[]
    for i,option in enumerate(options):
        check();cf=np.array(option,int);require(boundary_signature(cf)==boundary_signature(patch),'Actual external strip boundary changed')
        tr=v[cf];n=np.cross(tr[:,1]-tr[:,0],tr[:,2]-tr[:,0]);area=np.linalg.norm(n,axis=1)
        if np.any(area==0):histories.append(dict(ordinal=i,passed=False,reason='Actual zero-area triangle'));continue
        _,selfs=pair_contacts(Mesh(v,cf,native.frame),Mesh(v,cf,native.frame),0.,selfcheck=True,max_pairs=250000,max_tests=2000000)
        _,external=pair_contacts(Mesh(v,cf,native.frame),Mesh(v,static,native.frame),0.,max_pairs=250000,max_tests=2000000)
        ok=clear(selfs)and clear(external);den=np.sum((tr-np.roll(tr,-1,axis=1))**2,axis=(1,2));quality=float(np.min(2*np.sqrt(3)*area/den))
        histories.append(dict(ordinal=i,passed=ok,faces=cf.tolist(),self=selfs,static=external,minimum_triangle_quality=quality))
        if ok:possible.append((quality,i,cf))
    output.mkdir(parents=True);base=dict(schema='native-hg-fixed-upper-strip-repair/v1',source_candidate=report['candidate'],source_candidate_sha256=digest_file(report['candidate']),
        source_report=str(source/'candidate.json'),source_report_sha256=digest_file(source/'candidate.json'),source_plan=str(source/'plan.json'),source_plan_sha256=digest_file(source/'plan.json'),
        actual_defect_pairs=pairs.tolist(),domain_face_rows=ids.tolist(),boundary_loop=list(loop),complete14way_enumeration=histories,geometry_controls_moved=False,
        attachment_locations_changed=False,field_parameters_changed=False,anatomical_acceptance=False,installed=False)
    if not possible:
        atomic_json(output/'candidate.json',dict(**base,passed=False,reason='No valid fixed-boundary triangulation',elapsed_seconds=time.monotonic()-began));print('No valid fixed-boundary triangulation');return
    quality,ordinal,cf=max(possible,key=lambda item:item[0]);nf=f.copy();nf[ids]=cf
    words=native.attributes['normals'][native.faces].copy();signs=native.attributes.get('normalSigns',np.ones(len(native.vertices),dtype='i1'))[native.faces].copy()
    # Reuse the transported smooth field on identical boundary controls. A new
    # diagonal changes triangle interpolation, so only incompatible corner
    # normals receive an explicit actual-face chart; no geometry is hidden.
    vertex_words={}
    for face,wr,sr in zip(f,words,signs):
        for vertex,word,sign in zip(face,wr,sr):vertex_words.setdefault(int(vertex),(word,sign))
    for row,face in zip(ids,cf):
        for j,vertex in enumerate(face):words[row,j],signs[row,j]=vertex_words[int(vertex)]
    tr=v[nf];geometric=np.cross(tr[:,1]-tr[:,0],tr[:,2]-tr[:,0]);geometric/=np.linalg.norm(geometric,axis=1)[:,None]
    nn=words.astype(float)*signs[:,:,None];nn/=np.linalg.norm(nn,axis=2)[:,:,None];back=(nn*geometric[:,None,:]).sum(2)<=0
    fallback=np.rint(geometric*32767).astype('<i2');words[back]=np.broadcast_to(fallback[:,None,:],words.shape)[back];signs[back]=1
    candidate=Mesh(v[nf].reshape(-1,3),np.arange(nf.size).reshape(-1,3),native.frame,
        dict(normals=words.reshape(-1,3),normalSigns=signs.reshape(-1)),dict(source_geometry_sha256=native.geometry_sha,proposal='One exact fixed-boundary fourteen-way upper-strip repair',anatomical_acceptance=False))
    dest=output/Path(report['candidate']).name;candidate.save(dest);stored=(candidate.vertices-ORIGIN).astype('<f4');require(np.array_equal(stored.astype(float)+ORIGIN,candidate.vertices),'Actual reviewed Float32 coordinates changed')
    atomic_npz(output/(dest.stem+'-stored.npz'),positions=stored,faces=candidate.faces,geometryOrigin=ORIGIN,frame=np.array(candidate.frame),attribute_normals=candidate.attributes['normals'],attribute_normalSigns=candidate.attributes['normalSigns'])
    check();checks=audit(candidate,require_closed=True,max_pairs=250000,max_tests=2000000);contacts=[];sources=read(W/'source-v1/source-audit.json')
    others=[dict(key=r['partId'],path=r['meshSavePath'])for r in sources['contextParts']if r['sourcePart'].get('system')=='skeletal']
    others += [dict(key=k,path=str(W/'shared-coverings-exposed-apo-sl-v3/packed'/('TA98TONGUE_male_'+k+'.npz')))for k in ['mucosa','aponeurosis']]
    for r in others:
        check();arr,s=pair_contacts(candidate,Mesh.load(r['path']),0.,max_pairs=250000,max_tests=2000000);path=output/(r['key'].replace(':','_')+'-contacts.npz');atomic_npz(path,**arr);contacts.append(dict(counterpart=r['key'],source=r['path'],source_sha256=digest_file(r['path']),summary=s,witnesses=str(path),witness_sha256=digest_file(path)))
    quality_module=module('hg_final_quality','fair-tongue-extrinsic-envelopes.py');qq=quality_module.quality(v,nf,words*signs[:,:,None]);topo=checks['topology']
    passed=checks['passed']and len(topo['components'])==1 and topo['components'][0]['closed_orientable_genus']==0 and all(clear(r['summary'])for r in contacts)and qq['backward_corner_normals']==0
    result=dict(**base,passed=bool(passed),candidate=str(dest),candidate_sha256=digest_file(dest),selected_ordinal=ordinal,selected_minimum_triangle_quality=quality,
        all_positions_literal=bool(np.array_equal(candidate.triangles[np.setdiff1d(np.arange(len(nf)),ids)],native.triangles[np.setdiff1d(np.arange(len(nf)),ids)])),
        all_other_physical_faces_literal=True,external_strip_boundary_literal=boundary_signature(nf[ids])==boundary_signature(f[ids]),
        positive_actual_face_normal_corner_charts_repaired=int(back.sum()),actual_geometry_audit=checks,contacts=contacts,
        source_shape=report['actual_stored_shape'],actual_stored_shape=qq,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'candidate.json',result)
    review=read(source/'review-spec.json');review['meshes'][0]['mesh']=str(dest);atomic_json(output/'review-spec.json',review)
    print(__import__('json').dumps(dict(output=str(output),passed=bool(passed),self=checks['passed'],diagonal=ordinal,triangles=len(nf),corner_normal_repairs=int(back.sum()),contacts={r['counterpart']:r['summary']['crossing_face_pairs']for r in contacts},elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.source.resolve(),a.output.resolve())
