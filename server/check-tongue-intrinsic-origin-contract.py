"""Small regression using the actual saved, hash-bound curved TV source."""
from __future__ import annotations
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'server')]
from geometry_tools.common import atomic_json,digest_file
from tongue_basal_contract import corner_normal_mesh,require

def main():
    started=time.monotonic()
    spec=importlib.util.spec_from_file_location('origin_contract_packer',ROOT/'server/build-tongue-coupled-candidate.py')
    p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
    stage=ROOT/'.local-models/geometry-workbench/tongue-completion/curved-intrinsic-packing-inputs-v2-normal-field'
    construction=stage/'construction.json'
    row=next(r for r in json.loads(construction.read_text())['meshes'] if r['key']=='transverse-vertical-shared')
    mesh=corner_normal_mesh(row['mesh'],row['geometry_sha256'])
    origin,literal=p.intrinsic_storage_origin(row,mesh)
    require(literal and np.array_equal((mesh.vertices-origin).astype('<f4').astype(float)+origin,mesh.vertices),'Actual reviewed source failed literal replay')
    bank=bytearray();packed=p.pack_record(bank,mesh,dict(id='actual-reviewed-TV'),origin,literal=True)
    words,signs=p.old.normal_words(mesh)
    require(np.array_equal(np.frombuffer(bank,dtype='<i2',count=words.size,offset=packed['normals']).reshape(words.shape),words)
            and np.array_equal(np.frombuffer(bank,dtype='i1',count=signs.size,offset=packed['normalSigns']),signs)
            and np.array_equal(np.frombuffer(bank,dtype='<u4',count=mesh.faces.size,offset=packed['indices']).reshape(mesh.faces.shape),mesh.faces),
            'Actual packed buffer changed face indices or encoded normal words/signs')
    legacy,strict=p.intrinsic_storage_origin({},mesh)
    require(not strict and np.array_equal(legacy,p.old.ORIGIN),'Legacy implicit-origin behavior changed')
    require(not np.array_equal((mesh.vertices-legacy).astype('<f4').astype(float)+legacy,mesh.vertices),'Regression source no longer witnesses old-origin coordinate change')
    cases={}
    for name,field,value in [('wrong-origin','geometryOrigin',[0.,1.5,0.]),
                             ('mutated-mesh-hash','sha256','0'*64),
                             ('mutated-geometry-hash','geometry_sha256','0'*64),
                             ('mutated-proof-hash','storageOriginProof',dict(row['storageOriginProof'],sha256='0'*64))]:
        changed=copy.deepcopy(row);changed[field]=value
        try:p.intrinsic_storage_origin(changed,mesh)
        except (ValueError,AssertionError) as exc:cases[name]=str(exc)
        else:raise ValueError('Invalid origin binding was accepted: '+name)
    try:p.pack_record(bytearray(),mesh,dict(id='actual-TV-legacy-literal-test'),legacy,literal=True)
    except ValueError as exc:cases['literal-pack-rejects-legacy-origin']=str(exc)
    else:raise ValueError('Literal packer accepted changed coordinates')
    report=dict(schema='tongue-intrinsic-origin-contract-regression/v1',passed=True,
                source_construction=dict(path=str(construction),sha256=digest_file(construction)),
                reviewed_origin_literal=True,literal_face_indices_normal_words_and_signs=True,legacy_implicit_behavior_preserved=True,
                rejects=cases,elapsed_seconds=time.monotonic()-started)
    out=stage/'origin-contract-regression.json';atomic_json(out,report)
    print(json.dumps(report))

if __name__=='__main__':main()
