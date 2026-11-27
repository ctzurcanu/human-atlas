"""Scoped source proposal using explicit paired named digit bones, not a single artery guide."""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,complete_source_sections
from geometry_tools.algorithms.paired_bone_transport import paired_bone_transport


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['candidate','profiles','bone-guides','part','output']:p.add_argument('--'+name,required=True)
    p.add_argument('--digit',type=int,required=True,choices=range(1,6))
    p.add_argument('--surface-widths',help='Hash-bound measured source/recipient named-bone shaft gauges; changes carriage centres, never vessel calibre')
    p.add_argument('--surface-frame-transfer',action='store_true',help='Use actual named-bone surface axes and lengths instead of the bone-centroid curve')
    p.add_argument('--paired-vascular-reference',help='Explicit source and recipient paired palmar-vessel cross-section frames')
    args=p.parse_args();start=time.monotonic()
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in ['candidate','profiles','bone-guides']}
    data={k:json.loads(v.read_text())for k,v in paths.items()};out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new named-bone transport stage')
    report=copy.deepcopy(data['candidate']);row=next(r for r in report['rows']if r['partId']==args.part)
    original=Mesh.load(row['sourcePath'])
    if sha(row['sourcePath'])!=row['sourceSha256']:raise ValueError('Original source changed')
    current=Mesh.load(row['path'])
    if sha(row['path'])!=row['sha256']or not np.array_equal(original.faces,current.faces):raise ValueError('Candidate source ancestry changed')
    profile=next(r for r in data['profiles']['rows']if r['partId']==args.part)
    if sha(profile['path'])!=profile['sha256']:raise ValueError('Original section pairing changed')
    with np.load(profile['path'],allow_pickle=False)as z:pairing=z['retained_source_triangle_pairs'].copy()
    graph=source_profile_graph(original.vertices,original.faces,triangle_pairing=pairing,extend_annuli=True)
    _,physical=np.unique(original.vertices,axis=0,return_inverse=True)
    sections=complete_source_sections(graph,physical)
    guides=[g for g in data['bone-guides']['rows']if g['side']==row['side']]
    bones=next(g['namedBoneCorrespondences']for g in guides if g['digit']==args.digit)
    normals=[]
    for key in ['sourceCentre_m','targetCentre_m']:
        palm=np.array([g['namedBoneCorrespondences'][0][key]for g in guides]);_,_,vh=np.linalg.svd(palm-palm.mean(0))
        n=vh[-1];normals.append(n if n[2]>=0 else -n)
    width_scales=None
    if args.surface_widths:
        path=Path(args.surface_widths).resolve();widths=json.loads(path.read_text());paths['surface-widths']=path
        chosen=[next(r for r in widths['rows']if r['segment']==g['segment'])for g in bones]
        for fixed in chosen:
            for name in ['source','target']:
                if sha(fixed[name+'Path'])!=fixed[name+'Sha256']:raise ValueError('Measured named-bone surface changed')
        width_scales=[g['transverseWidthRatios']for g in chosen]
    if args.paired_vascular_reference:
        if args.surface_frame_transfer or args.surface_widths:raise ValueError('Choose one explicitly defined carriage model')
        from geometry_tools.algorithms.paired_bone_transport import paired_vascular_frame_carriage
        path=Path(args.paired_vascular_reference).resolve();reference=json.loads(path.read_text());paths['paired-vascular-reference']=path
        if reference['sourcePartId']!=args.part or reference['sourceMeshSha256']!=row['sourceSha256']:
            raise ValueError('Paired vascular reference belongs to another source')
        for name in ['sourceInventory','nativeGuide']:
            r=reference[name]
            if sha(r['path'])!=r['sha256']:raise ValueError('Paired source/reference course changed')
        if 'sourceCurves_m'in reference:
            from geometry_tools.algorithms.paired_bone_transport import paired_curve_carriage
            q,R,proof=paired_curve_carriage(original.vertices,sections,physical,
                    reference['sourceCurves_m'],reference['targetCurves_m'],reference['sourceAxis'],reference['targetAxis'],
                    np.mean(reference['sourcePair_m'],axis=0),np.mean(reference['targetPair_m'],axis=0),reference['axialScale'])
        else:
            q,R,proof=paired_vascular_frame_carriage(original.vertices,sections,physical,
                    reference['sourcePair_m'],reference['targetPair_m'],reference['sourceAxis'],reference['targetAxis'],reference['axialScale'])
    elif args.surface_frame_transfer:
        if not args.surface_widths:raise ValueError('Surface frames require explicit measured paired bone gauges')
        from geometry_tools.algorithms.paired_bone_transport import surface_frame_carriage
        from geometry_tools.open3d_backend import SurfaceIndex
        carriers=[]
        for fixed in chosen:
            a,b=Mesh.load(fixed['sourcePath']),Mesh.load(fixed['targetPath']);fa,fb=np.array(fixed['sourceFrame']),np.array(fixed['targetFrame'])
            axial=np.ptp(b.vertices@fb[:,2])/np.ptp(a.vertices@fa[:,2])
            carriers.append(dict(sourceFrame=fa,targetFrame=fb,sourceCentre=a.vertices.mean(0),targetCentre=b.vertices.mean(0),
                                 scales=[*fixed['transverseWidthRatios'],float(axial)],index=SurfaceIndex(a)))
        q,R,proof=surface_frame_carriage(original.vertices,sections,carriers,physical)
    else:
        q,R,proof=paired_bone_transport(original.vertices,sections,
                   [g['sourceCentre_m']for g in bones],[g['targetCentre_m']for g in bones],*normals,physical_ids=physical,
                   transverse_station_scales=width_scales)
    # Literal source aliases share the one named-bone field and section pose.
    for identity in np.unique(physical):
        ids=np.flatnonzero(physical==identity)
        if not np.all(q[ids]==q[ids[0]]):raise ValueError('Original source aliases diverged')
    gauge=max(float(np.abs(np.linalg.norm(q[s,None]-q[None,s],axis=2)-
                          np.linalg.norm(original.vertices[s,None]-original.vertices[None,s],axis=2)).max())for s in sections)
    if gauge>1e-12:raise ValueError('Original source section calibre changed')
    attrs={k:v.copy()for k,v in original.attributes.items()}
    n=attrs['normals'].astype(float)*attrs.get('normalSigns',np.ones(len(q)))[:,None]
    n=np.einsum('nij,nj->ni',R,n);n/=np.linalg.norm(n,axis=1)[:,None]
    attrs.update(normals=np.rint(n*32767).astype('<i2'),normalSigns=np.ones(len(q),'i1'))
    moved=Mesh(q,original.faces,original.frame,attrs,original.source);out.mkdir();file=out/'posed-part.npz';moved.save(file)
    row.update(path=str(file),sha256=sha(file),geometrySha256=moved.geometry_sha)
    for pair in report['pairs']:
        if args.part in [pair['firstPart'],pair['secondPart']]:
            pair.update(literalCoordinatesAndEffectiveNormals=False,packedCoordinatesAndEffectiveNormalWordsLiteral=False,requiresRecoupling=True)
    report.update(namedBoneTransport=proof,sourceAnnularPropagation=True,sourceFacesLiteral=True,
                  sourceSectionGaugeError_m=gauge,
                  complete=True,installed=False,installationAllowed=False,portsRequireRecoupling=True,
                  inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},
                  sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-start)
    if report['elapsedSeconds']>18:raise TimeoutError('18-second named-bone proposal deadline')
    write_json(out/'report.json',report);print(json.dumps(dict(sections=len(sections),sourceFacesRetained=True,elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
