"""Bounded, saved assessment of the installed liver and adjoining conduits."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sex', choices=['male', 'female'], required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); start = time.monotonic()
    out = args.output.resolve(); out.mkdir(exist_ok=False)
    cat = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas' / args.sex / 'atlas.json'
    a = json.loads(cat.read_text()); parts = {p['id']: p for p in a['parts']}
    frame = 'model:local-' + ('female-' if args.sex == 'female' else '') + 'ta98:assembled-metres'
    reader = AtlasReader(); meshes = {}
    def read(pid):
        if pid not in meshes:
            meshes[pid] = reader.read_part(a, parts[pid], frame=frame, catalogue_path=cat)
        return meshes[pid]
    host = parts['TA98PART:ZA:Liver:capsule']; full = read(host['id'])
    outer = Mesh(full.vertices, full.faces[:host['tissueDrawGroups'][0]['count']//3], frame)
    outer.save(out/'liver-exterior.npz'); bvh = BVH(outer)
    stomach = [p for p in a['parts'] if p['id'].startswith('GASTRIC:'+args.sex+':019:') and not p.get('suppressed')]
    contacts = []
    for p in stomach:
        m = read(p['id']); m.save(out/(p['id'].replace(':','_')+'.npz'))
        _, r = pair_contacts(outer, m, 0, first_bvh=bvh, max_pairs=200000, max_tests=4000000)
        contacts.append(dict(part=p['id'], **r))
        print(p['id'], r['crossing_face_pairs'], r['coplanar_overlap_face_pairs'], r['truncated'], flush=True)
    terms = {'A05.8.01.'+str(i).zfill(3) for i in [28,31,32,52,53,68]}
    possible = []
    for p in a['parts']:
        if p.get('ta98Term') in terms or any(w in p.get('name','').lower() for w in ['papillary process','caudate process','fibrous appendix']):
            possible.append({k:p.get(k) for k in ['id','name','ta98Term','suppressed','bounds']})
    roots = ['TA98HEPATIC:'+args.sex+':'+r for r in ['artery','portal','hepatic_vein','biliary']]
    adjacent = [p for p in a['parts'] if not p.get('suppressed') and not p.get('hepaticRegionOf') and p['id'] not in roots and any(w in p.get('name','').lower() for w in ['common hepatic artery','proper hepatic artery','coeliac trunk','celiac trunk','vena cava','cystic duct','bile duct','portal vein']) and not p['id'].startswith('HRA:')]
    summaries = []
    for pid in roots + [p['id'] for p in adjacent]:
        m=read(pid); m.save(out/(pid.replace(':','_').replace(' ','_')+'.npz'))
        summaries.append(dict(part=pid, name=parts[pid]['name'], vertices=len(m.vertices), faces=len(m.faces), bounds=m.bounds.tolist(), geometry_sha=m.geometry_sha))
    atomic_json(out/'assessment.json', dict(sex=args.sex, source_catalogue=str(cat), source_sha256=digest_file(cat), stomach_contacts=contacts, missing_source_candidates=possible, conduits=summaries, elapsed_seconds=time.monotonic()-start))
    print(json.dumps(dict(sex=args.sex, stomach_crossings=sum(r['crossing_face_pairs'] for r in contacts), candidates=possible, elapsed_seconds=time.monotonic()-start)),flush=True)
if __name__=='__main__':main()
