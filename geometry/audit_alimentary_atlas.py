"""Read-only oesophagus-to-rectum lumen audit of a local TA98 catalogue.

python3 -m geometry.audit_alimentary_atlas --catalogue .local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json --report reports/male-lumen.json

Extracts the actual innermost tissue faces, excluding finite annular rims. Does
not substitute route polylines, whole-organ bounds or generic exteriors for a
lumen. Exact-coordinate topology is supplemented by labelled, measured mouth
gaps; proximity is diagnostic and never establishes an anatomical join.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from .tube_continuity import continuity_audit,boundary_loops
from .tube_sections import TubeGeometryError


def _asset(root,url):
    if url.startswith('/local-models/'):p=root/'.local-models'/url[len('/local-models/'):]
    elif url.startswith('/'):p=root/'public'/url[1:]
    else:raise TubeGeometryError('expected a local model asset URL')
    p=p.resolve()
    if not p.is_relative_to(root.resolve()):raise TubeGeometryError('asset outside the atlas workspace')
    return p


def audit_catalogue(catalogue, root, maximum_gap_examples=16):
    path=Path(catalogue);root=Path(root);cat=json.loads(path.read_text());buffers={}
    vertices=[];faces=[];faceowners=[];sources=[];offset=0;unrepresented=[]
    for part in cat['parts']:
        pid=part['id']
        if part.get('suppressed'):continue
        wanted=(pid.startswith('OESOPHWALL:') and pid.endswith(':epithelium') or
            pid.startswith('GASTRIC:') and ':027:' in pid or
            ('ANALWALL:' in pid or 'ANALCANALWALL:' in pid) and pid.endswith(':epithelium') or
            ':BOWELWALL:' in pid and (pid.endswith(':epithelium') or pid.endswith(':circular-folds')) or
            part.get('caecocolicLayer')=='epithelium' or
            pid.startswith('GASTRIC:') and pid.endswith(':pyloric-connection'))
        if not wanted:continue
        chunk=cat['chunks'][part['chunk']]['url']
        if chunk not in buffers:buffers[chunk]=_asset(root,chunk).read_bytes()
        blob=buffers[chunk]
        p=np.frombuffer(blob,'<f4',part['vertexCount']*3,part['positions']).reshape(-1,3)
        f=np.frombuffer(blob,'<u4',part['indexCount'],part['indices']).reshape(-1,3)
        if 'normalSigns' in part:
            signs=np.frombuffer(blob,'i1',part['vertexCount'],part['normalSigns'])
            f=f[np.all(signs[f]<0,axis=1)]
        elif pid.endswith(':pyloric-connection') and part.get('pairedTissue'):
            # Explicit legacy connector schema: outer faces, reversed inner
            # faces, then separately grouped annular tissue rims (installer).
            g=part['tissueDrawGroups'][0]
            if g['start']!=0 or g['count']%6:raise TubeGeometryError('unrecognized legacy connector layout')
            half=g['count']//6;f=f[half:2*half]
        else:
            unrepresented.append(dict(part=pid,reason='no documented innermost-face ownership'));continue
        if not len(f):
            unrepresented.append(dict(part=pid,reason='no owned lumen triangles'));continue
        used=np.unique(f);remap=np.full(len(p),-1,int);remap[used]=np.arange(len(used))
        origin=np.array(part.get('geometryOrigin',part.get('sourceOffset',[0,0,0])),float)
        vertices.append(p[used].astype(float)+origin);faces.append(remap[f]+offset);faceowners.extend([pid]*len(f));offset+=len(used)
        sources.append(dict(part=pid,asset=chunk,asset_sha256=hashlib.sha256(blob).hexdigest(),lumen_faces=len(f)))
    if not vertices:raise TubeGeometryError('catalogue contains no identifiable finite-wall lumens')
    v=np.concatenate(vertices);f=np.concatenate(faces);owners=np.array(faceowners)
    report=continuity_audit(v,f)
    unique,inverse=np.unique(v,axis=0,return_inverse=True);welded=inverse[f]
    mouthrows=[];pairs=[]
    try:
        rings=boundary_loops(welded)
        for ring in rings:
            touching=np.any(np.isin(welded,ring),axis=1)
            mouthrows.append(dict(parts=sorted(set(owners[touching])),controls=len(ring),center_m=unique[ring].mean(axis=0).tolist()))
        for i in range(len(rings)):
            for j in range(i+1,len(rings)):
                a,b=unique[rings[i]],unique[rings[j]]
                d=np.linalg.norm(a[:,None]-b[None],axis=2)
                pairs.append(dict(first=i,second=j,center_gap_m=float(np.linalg.norm(a.mean(0)-b.mean(0))),
                    symmetric_boundary_gap_m=float(max(d.min(axis=0).max(),d.min(axis=1).max())),
                    first_parts=mouthrows[i]['parts'],second_parts=mouthrows[j]['parts'],establishes_join=False))
    except TubeGeometryError as err:report['mouth_localization_error']=str(err)
    # No finite anal-canal lining is present in the currently installed recipes.
    # Naming anus or an external muscle does not supply the missing conduit.
    anal=[p['part'] for p in sources if ('ANALWALL:' in p['part'] or 'ANALCANALWALL:' in p['part']) and p['part'].endswith(':epithelium')]
    report.update(catalogue_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),catalogue=str(path),
        lumen_sources=sources,labelled_mouths=mouthrows,nearest_mouth_pairs=sorted(pairs,key=lambda x:x['symmetric_boundary_gap_m'])[:maximum_gap_examples],
        unrepresented_sources=unrepresented,finite_anal_canal_lining_ids=anal,
        oesophagus_to_anus_pass=report['topological_continuity_pass'] and bool(anal),
        atlas_modified=False,diagnostic_scope='exact installed lumen topology; geometric intersections, sections and named end-collar review also required')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--catalogue',required=True);parser.add_argument('--root',default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--report',required=True);args=parser.parse_args()
    if Path(args.report).resolve()==Path(args.catalogue).resolve():parser.error('report must differ from catalogue')
    report=audit_catalogue(args.catalogue,args.root);out=Path(args.report);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(components=report['components'],mouths=report['mouths'],oesophagus_to_anus_pass=report['oesophagus_to_anus_pass'],
        finite_anal_canal_lining_ids=report['finite_anal_canal_lining_ids'],report=str(out),atlas_modified=False)))
    return 0 if report['oesophagus_to_anus_pass'] else 2


if __name__=='__main__':raise SystemExit(main())
