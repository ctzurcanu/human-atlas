"""Fixed gallbladder free-surface witnesses for the first donor placement."""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.review import prepare

WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
TARGET=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'

def main():
    began=time.monotonic();audit=json.loads((WORK/'source-audit.json').read_text());atlas=json.loads(TARGET.read_text());fit=WORK/'first-similarity';summary=[];surfaces={}
    for layer in ['covering','mucosa']:
        vertices=[];faces=[];owners=[];offset=0
        for r in audit['rows']:
            if not r['id'].startswith('TA98PART:ZA:Gallbladder:GALLWALL:'+layer+':'):continue
            part=next(p for p in atlas['parts'] if p['id']==r['id']);mesh=Mesh.load(ROOT/r['path']);chosen=[]
            for group in part['tissueDrawGroups']:
                if not group.get('neighbor'):chosen.extend(range(group['start']//3,(group['start']+group['count'])//3))
            assert chosen
            vertices.append(mesh.vertices);faces.append(mesh.faces[chosen]+offset);offset+=len(mesh.vertices);owners.extend([part['id']]*len(chosen))
        surface=Mesh(np.concatenate(vertices),np.concatenate(faces),'model:local-ta98:assembled-metres',source=dict(surfaceOnly=True,kind='gallbladder '+('external serosal surface'if layer=='covering'else'innermost luminal mucosal surface'),parentCatalogueSha256=digest_file(TARGET)))
        surface.save(fit/('fixed-'+layer+'-free-surface.npz'));surfaces[layer]=surface
        for donor in ['artery','vein']:
            candidate=Mesh.load(fit/(donor+'.npz'));arrays,proof=pair_contacts(candidate,surface,0.,max_pairs=250000,max_tests=2000000)
            atomic_npz(fit/(donor+'-'+layer+'-witness.npz'),**arrays)
            atomic_json(fit/(donor+'-'+layer+'-contacts.json'),proof)
            summary.append(dict(donor=donor,context=layer,proof=proof))
    atomic_json(fit/'context-checks.json',dict(schema='cystic-blood-first-fit-context/v1',checks=summary,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-began))
    spec=dict(frame='model:local-ta98:assembled-metres',meshes=[dict(mesh=str(fit/'artery.npz'),name='Donor cystic artery · first similarity only',color=[.82,.12,.1,1]),dict(mesh=str(fit/'vein.npz'),name='Donor cystic vein · first similarity only',color=[.13,.27,.8,1]),dict(mesh=str(fit/'fixed-covering-free-surface.npz'),name='Fixed male gallbladder outer surface',color=[.45,.56,.24,.55]),dict(mesh=str(WORK/'fixed-target-22.npz'),name='Fixed right hepatic arterial region',color=[.8,.19,.13,1]),dict(mesh=str(fit/'fixed-right-portal-region.npz'),name='Fixed right portal venous region',color=[.2,.27,.66,1])],views=[dict(name='front',direction=[0,0,1]),dict(name='right-oblique',direction=[1,.2,1]),dict(name='section',direction=[0,0,1],section=dict(normal=[0,0,1],point=surfaces['covering'].bounds.mean(0).tolist(),clear_outer=True))],resolution=640,wall_seconds=45)
    atomic_json(fit/'review-spec.json',spec);review=prepare(None,spec);atomic_json(fit/'review-queued.json',review)
    print(json.dumps(dict(output=str(fit),checks=[dict(donor=r['donor'],context=r['context'],proof=r['proof'])for r in summary],review=review,elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
