"""Create a bounded finger-course phase proposal; never authorize installation."""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,complete_source_sections
from geometry_tools.algorithms.source_control_clearance import source_control_basis
from geometry_tools.algorithms.section_axis_phase import phase_rigid_groups


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['candidate','profiles','binding','guides','part','output']:p.add_argument('--'+key,required=True)
    p.add_argument('--component',type=int,required=True);p.add_argument('--digit',type=int,required=True,choices=range(1,6))
    p.add_argument('--angle-degrees',type=float,required=True)
    args=p.parse_args();start=time.monotonic()
    paths={k:Path(getattr(args,k)).resolve()for k in ['candidate','profiles','binding','guides']}
    data={k:json.loads(v.read_text())for k,v in paths.items()};out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new bounded phase stage')
    row=next(r for r in data['candidate']['rows']if r['partId']==args.part)
    if sha(row['path'])!=row['sha256']or sha(row['sourcePath'])!=row['sourceSha256']:raise ValueError('Source ancestry changed')
    m,original=Mesh.load(row['path']),Mesh.load(row['sourcePath'])
    if not np.array_equal(m.faces,original.faces):raise ValueError('Literal original source faces required')
    owned=next(r for r in data['binding']['rows']if r['partId']==args.part and r['physicalComponent']==args.component)
    if not owned['anatomicallyReviewed']or args.digit not in owned['coveredDigits']:raise ValueError('Source digit ownership is unreviewed')
    fids=np.array(owned['sourceFaceIds']);raw=np.unique(m.faces[fids]);remap=np.full(len(m.vertices),-1,int);remap[raw]=np.arange(len(raw))
    bank=next(r for r in data['profiles']['rows']if r['partId']==args.part)
    if sha(bank['path'])!=bank['sha256']:raise ValueError('Source section bank changed')
    with np.load(bank['path'])as z:
        if not np.array_equal(z['source_faces'],original.faces):raise ValueError('Source section face indices are stale')
        pairing=z['retained_source_triangle_pairs'].copy()
    graph=source_profile_graph(original.vertices,original.faces,triangle_pairing=pairing,extend_annuli=True)
    physical=np.unique(original.vertices,axis=0,return_inverse=True)[1]
    sections=[remap[s]for s in complete_source_sections(graph,physical)if np.all(remap[s]>=0)]
    W,_,_,_=source_control_basis(m.vertices[raw],remap[m.faces[fids]],sections,physical[raw],micro_coherence_m=.00002,core_patch_diameter_m=.001)
    port=[]
    for pair in data['candidate']['pairs']:
        if pair['firstPart']==args.part:port.extend(pair['firstRaw'])
        if pair['secondPart']==args.part:port.extend(pair['secondRaw'])
    port=np.intersect1d(raw,port)
    bones=next(r['namedBoneCorrespondences']for r in data['guides']['rows']if r['side']==row['side']and r['digit']==args.digit)
    points=np.array([r['targetCentre_m']for r in bones]);curve=PchipInterpolator(np.arange(len(points)),points)
    t=np.linspace(-.5,len(points)-.75,257);axis=curve(t);tangent=curve.derivative()(t);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    q,R,proof=phase_rigid_groups(m.vertices[raw],W.argmax(1),axis,tangent,args.angle_degrees,remap[port])
    gauge=max(float(abs(np.linalg.norm(q[s,None]-q[None,s],axis=2)-np.linalg.norm(m.vertices[raw][s,None]-m.vertices[raw][None,s],axis=2)).max())for s in sections)
    if gauge>1e-12:raise ValueError('Original source section calibre changed')
    full=m.vertices.copy();full[raw]=q;attrs={k:v.copy()for k,v in m.attributes.items()}
    n=attrs['normals'][raw].astype(float)*attrs['normalSigns'][raw,None];n=np.einsum('nij,nj->ni',R,n);n/=np.linalg.norm(n,axis=1)[:,None]
    attrs['normals'][raw]=np.rint(n*32767).astype('<i2');attrs['normalSigns'][raw]=1
    attrs['normals'][port]=m.attributes['normals'][port];attrs['normalSigns'][port]=m.attributes['normalSigns'][port]
    if not np.array_equal(full[port],m.vertices[port]):raise ValueError('Protected source mouth moved')
    candidate=Mesh(full,m.faces,m.frame,attrs,m.source);out.mkdir();file=out/'phased-part.npz';candidate.save(file)
    report=copy.deepcopy(data['candidate']);r=next(r for r in report['rows']if r['partId']==args.part)
    r.update(path=str(file),sha256=sha(file),geometrySha256=candidate.geometry_sha)
    report.update(phaseProposal=proof,sourceSectionGaugeError_m=gauge,complete=True,installed=False,installationAllowed=False,
                  inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},sourceCodeSha256=sha(__file__),
                  algorithmSha256=sha(Path(__file__).parents[1]/'algorithms/section_axis_phase.py'),elapsedSeconds=time.monotonic()-start)
    if report['elapsedSeconds']>18:raise TimeoutError('18-second source phase stage')
    write_json(out/'report.json',report);print(json.dumps(dict(phaseDegrees=args.angle_degrees,sourceSectionGaugeError_m=gauge,elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
