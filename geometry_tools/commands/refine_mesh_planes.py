"""Exact shared-edge refinement without moving a mesh or inventing caps.

Usage: python -m geometry_tools.commands.refine_mesh_planes SOURCE.npz --levels 0.1,0.2
       --axis 1 --output NEW_DIRECTORY
The output includes source-edge maps for UV/normal interpolation; it is not a
runtime replacement bank. Original point IDs and original surface are retained.
"""
from pathlib import Path
import argparse,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from geometry_tools.course_io import sha,load_arrays,write_json
from geometry_tools.algorithms.parallel_plane_refinement import refine_parallel_planes


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('source');parser.add_argument('--levels',required=True);parser.add_argument('--axis',type=int,default=1);parser.add_argument('--output',required=True);args=parser.parse_args();start=time.monotonic();out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new output directory')
    source=load_arrays(args.source);levels=np.array([float(s)for s in args.levels.split(',')]);r=refine_parallel_planes(source['vertices'],source['faces'],levels,axis=args.axis,time_limit_seconds=10.);out.mkdir();np.savez_compressed(out/'refinement.npz',**{k:v for k,v in r.items()if isinstance(v,np.ndarray)});write_json(out/'report.json',dict(schema='exact-plane-mesh-refinement/v1',sourcePath=str(Path(args.source).resolve()),sourceSha256=sha(args.source),arraysSha256=sha(out/'refinement.npz'),sourceCodeSha256=sha(__file__),originalVertices=len(source['vertices']),vertices=len(r['vertices']),originalFaces=len(source['faces']),faces=len(r['faces']),sourceFaceAreaAndOrientationReplayed=True,originalPointsLiteral=bool(np.array_equal(r['vertices'][:len(source['vertices'])],source['vertices'])),capsAdded=0,sourceAttributeMapsRetained=True,runtimeBankProduced=False,installed=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start));print(str(out/'report.json'))


if __name__=='__main__':main()
