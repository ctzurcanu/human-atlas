"""One bounded volume/topology feasibility measurement; never an install producer."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit,topology
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def read(p):return json.loads(Path(p).read_text())
def require(x,msg):
    if not x:raise ValueError(msg)
def main(out):
    began=time.monotonic();fp=W/'female-axis-tongue-fit-v1/fit.json';fit=read(fp)
    require(out.is_relative_to(W)and not out.exists(),'Use a new feasibility directory');out.mkdir(parents=True)
    origin=np.array([0.,1.5,0.]);clearance=.0001;hosts=[];solids=[]
    for row in fit['target_context']:
        if not any(s in row['partId']for s in [':hyoid:','mylohyoid','geniohyoid']):continue
        require(time.monotonic()-began<45,'Feasibility stage deadline approached')
        require(digest_file(row['mesh'])==row['sha256'],'Fixed host changed');m=Mesh.load(row['mesh'])
        check=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000)
        require(check['passed'],'A proposed exclusion host is not a valid closed self-free solid: '+row['partId'])
        s,_=make_solid(m.vertices-origin,m.faces);require(s.status()==md.Error.NoError,'Invalid host solid')
        hosts.append(dict(partId=row['partId'],mesh=row['mesh'],sha256=row['sha256'],geometry_sha256=m.geometry_sha,audit=check))
        solids.append(s)
    require(len(solids)==4,'Require actual hyoid and three closed oral-floor hosts')
    occupied=md.Manifold.batch_boolean(solids,md.OpType.Add)
    require(occupied.status()==md.Error.NoError,'Fixed-host union failed')
    cutter=occupied.minkowski_sum(md.Manifold.cube([2*clearance]*3,center=True))
    require(cutter.status()==md.Error.NoError,'Declared numerical clearance construction failed')
    raw=cutter.to_mesh64();atomic_npz(out/'classifier-only-cutter.npz',vertices=np.asarray(raw.vert_properties)[:,:3]+origin,faces=np.asarray(raw.tri_verts,int))
    results=[]
    for i,row in enumerate(fit['parts']):
        if not any(s in row['partId']for s in ['Genioglossus','Hyoglossus']):continue
        require(time.monotonic()-began<45,'Feasibility stage deadline approached');require(digest_file(row['mesh'])==row['sha256'],'Axis source changed')
        m=Mesh.load(row['mesh']);check=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000);require(check['passed'],'Axis source not valid')
        solid,_=make_solid(m.vertices-origin,m.faces);result=solid-cutter;require(result.status()==md.Error.NoError and not result.is_empty(),'Exclusion failed')
        before=float(solid.volume());after=float(result.volume());loss=(before-after)/before
        raw=result.to_mesh64();v=np.asarray(raw.vert_properties)[:,:3]+origin;f=np.asarray(raw.tri_verts,int);probe=Mesh(v,f,m.frame,{},dict(kind='unqualified-exclusion-feasibility-only'))
        topo=topology(probe);path=out/('probe-'+str(i)+'.npz');probe.save(path)
        rowout=dict(partId=row['partId'],source=row['mesh'],source_sha256=row['sha256'],source_geometry_sha256=m.geometry_sha,source_audit=check,probe=str(path),probe_sha256=digest_file(path),probe_geometry_sha256=probe.geometry_sha,before_m3=before,after_m3=after,removed_fraction=loss,removed_percent=loss*100,topology=topo,volume_loss_within_five_percent=loss<=.05,topology_unchanged=topo['closed_oriented_manifold']and len(topo['components'])==1 and topo['components'][0]['closed_orientable_genus']==0)
        results.append(rowout);atomic_json(out/('progress-'+str(i)+'.json'),rowout)
    acceptable=all(r['volume_loss_within_five_percent']and r['topology_unchanged']for r in results)
    report=dict(schema='female-tongue-single-exclusion-feasibility/v1',source_fit=str(fp),source_fit_sha256=digest_file(fp),construction_script_sha256=digest_file(__file__),hosts=hosts,parts=results,clearance_m=clearance,clearance_metric='L-infinity cube dilation; minimum Euclidean separation 0.1mm and maximum directional dilation sqrt(3)*0.1mm',maximum_directional_dilation_m=float(np.sqrt(3)*clearance),geometry_reconstruction_qualified=False,normal_charts_created=False,source_attachment_claim='No shared cortical face or histological footprint inferred',preliminary_volume_topology_feasible=acceptable,full_crossing_self_cover_and_render_required=True,fixed_skeleton=True,qualifying_candidate=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'evidence.json',report);print(json.dumps(dict(output=str(out),parts=[{k:r[k]for k in ['partId','removed_percent','volume_loss_within_five_percent','topology_unchanged']}for r in results],feasible=acceptable,elapsed_seconds=report['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
