"""Finite, source-qualified cardiac reference geometry in one canonical frame.

Native tissue regions remain literal intersections of their existing tissue.
Pericardial strata share corresponding boundaries and annular opening rims.
These are reference geometries, not replacements for the installed whole heart.
"""
from pathlib import Path
import argparse,json,runpy,sys,time
import numpy as np
import manifold3d as md
import igl
from scipy.spatial import ConvexHull
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.geometry_instance import instance_pose
types=runpy.run_path(str(ROOT/'server/heart-reference-types.py'))
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--guides',required=True);args=ap.parse_args()
started=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';out=base/args.stage;out.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};groups={g['key']:g for g in a['finiteCardiacIntegration']['termGroups']}
bundle=json.loads(path.with_name('cardiac-label-sites.json').read_text());sites={g['term']:g for g in bundle['integration']['termGroups']};reader=AtlasReader();cache={};solids={}
ep=parts['TA98HEART:male:epicardium'];R,scale=instance_pose(ep);origin=np.asarray(ep['geometryOrigin'])
local=lambda q:(np.asarray(q)-origin)@R/scale
world=lambda q:(np.asarray(q)@R.T)*scale+origin
def native(pid):
    if pid not in cache:
        m=reader.read_part(a,parts[pid],frame='assembled-metres',catalogue_path=path);cache[pid]=Mesh(local(m.vertices),m.faces,'heart-native-local-metres')
    return cache[pid]
def solid(pid):
    if pid not in solids:
        m=native(pid);v,remap=np.unique(m.vertices,axis=0,return_inverse=True);f=remap[m.faces]
        s=md.Manifold(md.Mesh64(np.array(v,dtype=float,order='C',copy=True),np.array(f,dtype='uint64',order='C',copy=True)))
        if str(s.status())!='Error.NoError':raise ValueError('Native tissue is not a finite Boolean carrier: '+pid)
        solids[pid]=s
    return solids[pid]
def q(term):return local(sites[term]['anchor']['point'])
def center(key):return np.mean(np.concatenate([native(pid).vertices for pid in groups[key]['partIds']]),axis=0)
def ellipsoid(center,radii,axes=None):
    transform=np.eye(4);transform[:3,:3]=(np.eye(3)if axes is None else axes.T)@np.diag(radii);transform[:3,3]=center
    return md.Manifold.sphere(1,24).transform(transform[:3])
def as_mesh(s):
    if str(s.status())!='Error.NoError' or s.is_empty():raise ValueError('Empty or invalid finite reference')
    m=s.to_mesh64();return Mesh(np.asarray(m.vert_properties[:,:3]),np.asarray(m.tri_verts,dtype=int),'heart-native-local-metres')
def from_arrays(v,f):
    volume=np.einsum('ij,ij->i',v[f[:,0]],np.cross(v[f[:,1]],v[f[:,2]])).sum()/6
    if volume<0:f=f[:,::-1]
    result=md.Manifold(md.Mesh64(np.array(v,dtype=float,order='C',copy=True),np.array(f,dtype='uint64',order='C',copy=True)))
    if str(result.status())!='Error.NoError':raise ValueError('Constructed finite reference has invalid topology')
    return result
def tube(points,radius,segments=10,samples=24):
    points=np.asarray(points);length=np.linalg.norm(np.diff(points,axis=0),axis=1);t=np.r_[0,np.cumsum(length)];t/=t[-1];u=np.linspace(0,1,samples);curve=CubicSpline(t,points,bc_type='natural')(u);tangent=np.gradient(curve,axis=0);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    frames=[];previous=None
    for direction in tangent:
        x=np.cross(direction,[1,0,0]if abs(direction[0])<.8 else[0,1,0])if previous is None else previous-direction*(previous@direction)
        x/=np.linalg.norm(x);y=np.cross(direction,x);frames.append((x,y));previous=x
    theta=np.arange(segments)*2*np.pi/segments;v=np.array([p+radius*(np.cos(theta)[:,None]*x+np.sin(theta)[:,None]*y)for p,(x,y)in zip(curve,frames)]).reshape(-1,3);f=[]
    for i in range(samples-1):
        for j in range(segments):
            aa=i*segments+j;bb=i*segments+(j+1)%segments;cc=(i+1)*segments+(j+1)%segments;dd=(i+1)*segments+j;f.extend([[aa,bb,cc],[aa,cc,dd]])
    for j in range(1,segments-1):f.extend([[0,j+1,j],[(samples-1)*segments,(samples-1)*segments+j,(samples-1)*segments+j+1]])
    return from_arrays(v,np.asarray(f))
def finite_sheet(mid,faces,thickness):
    mid=np.asarray(mid);faces=np.asarray(faces);norm=np.zeros_like(mid);n=np.cross(mid[faces[:,1]]-mid[faces[:,0]],mid[faces[:,2]]-mid[faces[:,0]])
    for k in range(3):np.add.at(norm,faces[:,k],n)
    norm/=np.maximum(np.linalg.norm(norm,axis=1)[:,None],1e-15)
    return layered_sheet(mid,faces,norm,[0,thickness])[0]
def layered_sheet(mid,faces,norm,depths):
    count=len(mid);fields=[mid+d*norm for d in depths];edges=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]);_,inverse,uses=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True);boundary=edges[uses[inverse]==1];result=[]
    for outer,inner in zip(fields[:-1],fields[1:]):
        v=np.r_[outer,inner];f=np.r_[faces[:,::-1],faces+count]
        side=np.array([[x,y,y+count]for x,y in boundary]+[[x,y+count,x+count]for x,y in boundary],dtype=int).reshape(-1,3)
        result.append(from_arrays(v,np.r_[f,side]))
    return result
def leaflet(center,normal,radius,thickness,span=1.9):
    normal=np.asarray(normal);normal/=np.linalg.norm(normal);x=np.cross(normal,[1,0,0]if abs(normal[0])<.8 else[0,1,0]);x/=np.linalg.norm(x);y=np.cross(normal,x)
    angles=np.linspace(-span/2,span/2,15);rim=np.asarray(center)+radius*(np.cos(angles)[:,None]*x+np.sin(angles)[:,None]*y);tip=np.asarray(center)+radius*.25*x+normal*radius*.2;mid=np.r_[rim,tip[None]];f=np.array([[i,i+1,len(rim)]for i in range(len(rim)-1)])
    return finite_sheet(mid,f,thickness)
def nearest(pid,points,peer=None):
    m=native(pid);faces=m.faces
    if peer:
        spans=[np.arange(d['start']//3,(d['start']+d['count'])//3)for d in parts[pid].get('tissueDrawGroups',[])if d.get('neighbor')==peer]
        if spans:faces=faces[np.concatenate(spans)]
    return igl.point_mesh_squared_distance(np.asarray(points).reshape(-1,3),m.vertices,faces)[2]
landmarks=json.loads((base/'cardiac-label-landmarks-v604/report.json').read_text());stats={r['cid']:r for r in landmarks['landmarks']};guide_meta=json.loads((base/args.guides/'construction.json').read_text());guide_by_term={r['term']:r for r in guide_meta['records']}
results={};issues=[]
def save(term,s,method,**extra):
    # Boolean slivers below the Float32 coordinate precision must be removed
    # before conversion. The one-micrometre bound is recorded and audited.
    mesh=as_mesh(s.as_original().simplify(.000001));mesh.save(out/(term+'.npz'));results[term]=dict(type='finite-mesh',file=term+'.npz',vertices=len(mesh.vertices),triangles=len(mesh.faces),method=method,bufferSimplificationToleranceMetres=.000001,**extra)
    print(term,len(mesh.faces),'triangles',round(time.monotonic()-started,2),'s',flush=True)

# Source-derived material networks are constructed independently of donor
# triangulation. Their anatomical extent still needs placement review.
for term,record in guide_by_term.items():
    if record.get('status')!='constructed-guide':continue
    try:
        if record['triangles']<128:raise ValueError('Coarse sampling lost the fine anatomical network')
        guide=Mesh.load(base/args.guides/(term+'.npz'));g=sites[term]
        # One proper source-to-heart registration is shared by the complete
        # conduction cohort; never translate its branches independently.
        s=from_arrays(guide.vertices,guide.faces)
        if term.startswith('A12.1.06.')or term in ['A12.1.00.025','A12.1.00.026']:
            inside=s^solid(g['anchor']['partId'])
            if inside.is_empty()or inside.volume()<s.volume()*.25:raise ValueError('Registered guide does not sufficiently occupy its correct native tissue')
            s=inside
        save(term,s,'independent anatomical material guide',referenceCid=record['cid'],guidePitchMetres=record['pitch_metres'])
    except ValueError as e:issues.append(dict(term=term,reason=str(e)))

# Region controls use actual finite native tissue, not their outer faces alone.
for term,cid in [('A12.1.01.002',20454),('A12.1.03.002',23809),('A12.1.02.008',16459),('A12.1.00.025',20644),('A12.1.00.026',23777)]:
    if term in results:continue
    try:
        st=stats[cid];extent=np.ptp(np.asarray(st['principal_extent_90_m']),axis=0)/(2*scale);axes=np.asarray(st['principal_axes'])@R
        # The guide bounds a region of the actual carrier; it cannot introduce
        # a solid appendage into an atrial lumen or manufacture fibrous tissue.
        roi=ellipsoid(q(term),np.maximum(extent,.001),axes)
        s=solid(sites[term]['anchor']['partId'])^roi
        save(term,s,'finite native tissue region bounded by source landmark extent',referenceCid=cid)
    except ValueError as e:issues.append(dict(term=term,reason=str(e)))

# The anterior/posterior LV papillary references already own actual finite
# stored source tissue; no new proxy or duplicate papillary body is added.

ra='TA98HEART:male:myocardium:right-atrium';rv='TA98HEART:male:myocardium:right-ventricle';la='TA98HEART:male:myocardium:left-atrium';lv='TA98HEART:male:myocardium:left-ventricle';ias='TA98HEART:male:myocardium:interatrial-septum';ivs='TA98HEART:male:myocardium:interventricular-muscular-septum';endo='TA98HEART:male:endocardium'
tv=center('tricuspid-valve');mv=center('mitral-valve');av=center('aortic-valve');pv=center('pulmonary-valve')
for term in [t for t in sites if types['required_type'](t)=='finite-mesh' and t not in results and t not in types['EXISTING_MESHES'] and not t.startswith('A12.1.08.')]:
    try:
        anchor=q(term);carrier=sites[term]['anchor']['partId'];method='anatomically scoped finite reference construction';extra={}
        if term=='A12.1.00.021':
            s=solid(lv)^ellipsoid(anchor,[.009/scale,.007/scale,.009/scale]);method='literal apical myocardial volume; fibre organization remains illustrative'
        elif term=='A12.1.00.024':
            pap=native(a['finiteCardiacIntegration']['papillaryPartIds'][-1]).vertices.mean(0);septal=nearest(ivs,[pap+[-.012,.004,0]])[0];s=tube([septal,(septal+pap)/2+[0,.002,0],pap],.00035/scale);extra['variableExample']=True;method='representative false tendon between septum and papillary tissue; no leaflet attachment'
        elif term in ['A12.1.00.028','A12.1.00.029']:
            points=np.array([tv,.6*tv+.4*pv,pv])if term.endswith('.028')else np.array([q('A12.1.01.015'),anchor,q('A12.1.06.004')]);s=tube(nearest(carrier,points),.00045/scale);method='finite fibrous band following its native attachment neighborhood'
        elif term in ['A12.1.01.003','A12.1.02.007']:
            points=np.array([q('A12.1.01.011'),anchor,q('A12.1.01.010')])if term.endswith('.003')else np.array([tv,anchor,pv]);s=tube(nearest(carrier,points,endo),.0012/scale);method='finite crest along the corresponding native endocardial boundary'
        elif term=='A12.1.01.014':
            s=solid(ias)^ellipsoid(anchor,[.005/scale,.005/scale,.002/scale]);method='literal superior interatrial-fold tissue; no invented discrete tubercle'
        elif term in ['A12.1.01.015','A12.1.01.016']:
            normal=q('A12.1.01.005')-anchor;s=leaflet(anchor,normal,.0035/scale if term.endswith('.015')else .0015/scale,.00015/scale);method='finite partial crescent valve leaflet; the remainder of the ostium stays uncovered'
        elif term in ['A12.1.02.014','A12.1.04.016']:
            key='pulmonary-valve'if term.startswith('A12.1.02')else'aortic-valve';c=center(key);pieces=[]
            for pid in groups[key]['partIds']:
                m=native(pid);node=m.vertices[np.linalg.norm(m.vertices-c,axis=1).argmin()];pieces.append(ellipsoid(node,[.00045/scale,.00035/scale,.00045/scale]))
            s=md.Manifold.batch_boolean(pieces,md.OpType.Add);method='three finite leaflet-margin nodules, individually located on their actual cusps'
        elif term=='A12.1.03.003':
            anchor=q('A12.1.03.002');pieces=[]
            for i in range(4):
                pp=nearest(la,[anchor+[i*.0012/scale,-.004/scale,0],anchor+[i*.0012/scale,0,.002/scale],anchor+[i*.0012/scale,.004/scale,0]],endo);pieces.append(tube(pp,.00045/scale))
            s=md.Manifold.batch_boolean(pieces,md.OpType.Add);method='four finite pectinate ridges in the native left-auricular endocardial neighborhood'
        elif term=='A12.1.03.005':
            s=solid(ias)^ellipsoid(anchor,[.006/scale,.005/scale,.0008/scale]);method='finite existing septal flap/closure tissue; no patent adult foramen'
        elif term=='A12.1.04.006':
            pieces=[]
            for pid in groups['mitral-valve']['partIds']:
                m=native(pid);v=m.vertices;_,_,axes=np.linalg.svd(v-v.mean(0),full_matrices=False);extent=np.ptp((v-v.mean(0))@axes.T,axis=0);target=v[np.argmax((v-mv)@axes[0])];pieces.append(solid(pid)^ellipsoid(target,np.maximum(extent*.18,.0006/scale),axes))
            s=md.Manifold.batch_boolean(pieces,md.OpType.Add);method='finite native mitral commissural leaflet regions; no invented accessory leaflet'
        elif term in ['A12.1.06.003','A12.1.06.004']:
            cid=types['SOURCE_CIDS'][term];st=stats[cid];extent=np.ptp(np.asarray(st['principal_extent_90_m']),axis=0)/(2*scale)
            if term.endswith('.004'):
                carrier=ias;target=q(term);m=native(carrier);tri=m.triangles;centers=tri.mean(1);order=np.argsort(np.linalg.norm(centers-target,axis=1))[:12];options=[]
                for index in order:
                    tt=tri[index];normal=np.cross(tt[1]-tt[0],tt[2]-tt[0]);normal/=max(np.linalg.norm(normal),1e-15);long=tv-centers[index];long-=normal*(long@normal);long/=max(np.linalg.norm(long),1e-15);axes=np.array([long,np.cross(normal,long),normal])
                    for sign in [-1,1]:
                        seed=centers[index]+sign*normal*.00025/scale;candidate=solid(carrier)^ellipsoid(seed,np.array([.0025,.0013,.0006])/scale,axes)
                        if not candidate.is_empty():options.append((candidate.volume(),candidate))
                if not options:raise ValueError('No finite native AV-nodal tissue region')
                s=max(options,key=lambda x:x[0])[1]
                if s.volume()*scale**3<.0000000005:raise ValueError('AV-nodal reference collapses below 0.5 cubic mm')
                method='finite AV-nodal reference in native inferior atrial septal tissue at the Koch-apex neighborhood; histological extent is illustrative'
            else:
                toward=native(carrier).vertices.mean(0)-anchor;toward/=max(np.linalg.norm(toward),1e-15);anchor=anchor+toward*.0004/scale
                s=solid(carrier)^ellipsoid(anchor,np.maximum(extent,.0005/scale),np.asarray(st['principal_axes'])@R);method='source-sized nodal volume contained in its anatomically corresponding atrial tissue'
            extra['carrierPartId']=carrier
        elif term=='A12.1.06.008':
            pieces=[]
            for ventricle in [rv,lv]:
                host=native(ventricle);low=host.vertices[np.argmin(host.vertices[:,1])];start=nearest(ventricle,[q('A12.1.06.007'if ventricle==lv else'A12.1.06.006')],endo)[0]
                for i in range(6):
                    theta=(i-2.5)*np.pi/8;end=nearest(ventricle,[low+[np.cos(theta)*.010/scale,.010/scale,np.sin(theta)*.010/scale]],endo)[0];midpoint=nearest(ventricle,[(start+end)/2+[np.sin(theta)*.006/scale,0,0]],endo)[0];pieces.append(tube([start,midpoint,end],.00018/scale,segments=8,samples=18)^solid(ventricle))
            s=md.Manifold.batch_boolean(pieces,md.OpType.Add);method='finite branching subendocardial reference network following both native ventricular walls; terminal density is illustrative'
        elif term in ['A12.1.06.005','A12.1.06.006','A12.1.06.007']:
            start=q('A12.1.06.004')if term.endswith('.005')else q('A12.1.06.005');end=anchor;points=nearest(carrier,[start,(start+end)/2,end]);s=tube(points,.0003/scale)^solid(carrier);method='finite intramural conducting band clipped to native septal tissue; source distribution remains qualified'
        elif term in ['A12.1.01.008','A12.1.02.020','A12.1.02.021','A12.1.04.011']:
            pieces=[];extent=.006/scale
            for i in range(5):
                pp=nearest(carrier,[anchor+[-extent,i*.001/scale,0],anchor+[0,i*.001/scale,.002/scale],anchor+[extent,i*.001/scale,0]],endo);pieces.append(tube(pp,.00065/scale))
            s=md.Manifold.batch_boolean(pieces,md.OpType.Add);method='finite source-neighborhood myocardial ridge references'
        else:raise ValueError('No anatomically defined finite construction')
        save(term,s,method,**extra)
    except (ValueError,RuntimeError)as e:issues.append(dict(term=term,reason=str(e)))

# A single source-bound pericardial carrier, with shared corresponding layers.
# The fixed heart is contained; actual outgoing root courses define open ports.
try:
    outer=native('TA98HEART:male:epicardium');free=np.concatenate([np.arange(d['start']//3,(d['start']+d['count'])//3)for d in ep['tissueDrawGroups']if not d.get('neighbor')]);cloud=outer.triangles[free].reshape(-1,3)
    hull=md.Manifold.hull_points(cloud).refine(2).to_mesh64();mid=np.asarray(hull.vert_properties[:,:3]);faces=np.asarray(hull.tri_verts,dtype=int)
    norm=np.zeros_like(mid);fn=np.cross(mid[faces[:,1]]-mid[faces[:,0]],mid[faces[:,2]]-mid[faces[:,0]])
    for k in range(3):np.add.at(norm,faces[:,k],fn)
    norm/=np.linalg.norm(norm,axis=1)[:,None];centers=mid[faces].mean(1);keep=np.ones(len(faces),dtype=bool);ports=[]
    for pid in a['finiteCardiacIntegration']['rootPartIds']:
        m=native(pid);c=m.vertices.mean(0);_,_,axes=np.linalg.svd(m.vertices-c,full_matrices=False);axis=axes[0];t=(m.vertices-c)@axis;radial=(m.vertices-c)-t[:,None]*axis;radius=float(np.quantile(np.linalg.norm(radial,axis=1),.85));x=(centers-c)@axis;rad=np.linalg.norm((centers-c)-x[:,None]*axis,axis=1);mask=(rad<radius+.0015/scale)&(x>=t.min()-.001/scale)&(x<=t.max()+.001/scale)
        keep&=~mask;ports.append(dict(partId=pid,removedCarrierFaces=int(mask.sum())))
    faces=faces[keep]
    # Adjacent root ports can meet at a single carrier vertex. Expand those
    # local ports by their incident facets to remove pinched boundary links.
    for _ in range(8):
        edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1);keys,counts=np.unique(edges,axis=0,return_counts=True);boundary=keys[counts==1];ids,degree=np.unique(boundary,return_counts=True);bad=ids[degree!=2]
        if not len(bad):break
        faces=faces[~np.any(np.isin(faces,bad),axis=1)]
    assert not len(bad),'Pinched pericardial port boundaries'
    cells=layered_sheet(mid,faces,norm,[.001/scale,.00103/scale,.0012/scale,.0018/scale])
    save('A12.1.08.008',cells[0],'finite parietal serosal stratum on a shared pericardial carrier',thicknessMetres=.00003)
    save('A12.1.08.009',cells[1],'finite parietal subserosal stratum sharing the serosal boundary',thicknessMetres=.00017)
    save('A12.1.08.002',cells[2],'finite fibrous pericardial wall with native-root opening ports',thicknessMetres=.0006,ports=ports)
    save('A12.1.08.006',cells[0]+cells[1],'finite parietal serous wall assembled from its two corresponding strata')
    # Thin attachment sheets use actual root/heart landmarks in the same frame.
    left=q('A12.1.03.002');veins=[native(pid).vertices.mean(0)for pid in a['finiteCardiacIntegration']['rootPartIds']if 'Left superior pulmonary vein'in parts[pid]['name']];v=veins[0];fold=np.array([left,v,(left+v)/2+[0,.004/scale,-.002/scale]])
    save('A12.1.08.010',finite_sheet(fold,np.array([[0,1,2]]),.00015/scale),'finite triangular Marshall fold between the left auricle and left superior pulmonary-venous neighborhood')
    p0=q('A12.1.08.004');trachea=next(p for p in a['parts']if p['name']=='Trachea'and not p.get('suppressed'));t=native(trachea['id']).vertices;carina=t[np.argmin(t[:,1])];diaphragm=next(p for p in a['parts']if p['id']=='TA98PART:ZA:Diaphragm:Tendon:central');insertion=nearest(diaphragm['id'],[p0])[0];sheet=np.array([carina+[-.003/scale,0,0],carina+[.003/scale,0,0],insertion+[-.010/scale,0,0],insertion+[.010/scale,0,0]])
    save('A12.1.08.004',finite_sheet(sheet,np.array([[0,2,1],[1,2,3]]),.0002/scale),'finite dorsal bronchopericardial sheet bounded by actual tracheal-carinal and central diaphragmatic neighborhoods; lateral reflections remain illustrative',fixedAttachmentPartIds=[trachea['id'],diaphragm['id']])
    anchor=q('A12.1.08.003');pieces=[];sternal=[p for p in a['parts']if not p.get('suppressed')and p.get('system')=='skeletal'and any(k in p['name'].lower()for k in ['manubrium','body of sternum','xiphoid'])]
    for sign in [-1,1]:
        heart_point=nearest('TA98HEART:male:epicardium',[anchor+[0,sign*.015/scale,0]])[0]
        if sternal:
            candidates=[nearest(p['id'],[heart_point])[0]for p in sternal];attachment=min(candidates,key=lambda x:np.linalg.norm(x-heart_point))
        else:raise ValueError('No fixed sternal attachment source')
        pieces.append(tube([heart_point,(heart_point+attachment)/2+[0,.001/scale,0],attachment],.00065/scale))
    save('A12.1.08.003',md.Manifold.batch_boolean(pieces,md.OpType.Add),'two finite sternopericardial bands with endpoints on actual fixed sternal and pericardial surfaces')
except (ValueError,RuntimeError)as e:issues.append(dict(term='pericardial-assembly',reason=str(e)))

# Empty reference spaces are clipped against actual native tissue. They are
# rendered only as unfilled boundary cages; their crop faces are not tissue.
for term in sorted(types['SPACES']):
    try:
        center_point=q(term);radii=np.array([.012,.006,.006])/scale
        if term=='A12.1.04.009':center_point=.7*av+.3*center('left-ventricle');radii=np.array([.008,.012,.008])/scale
        roi=ellipsoid(center_point,radii)
        relevant=[]
        for pid in set(sum([a['finiteCardiacIntegration'].get(k,[])for k in ['wallPartIds','valvePartIds','papillaryPartIds','chordalPartIds','rootPartIds']],[])):
            m=native(pid);lo=m.vertices.min(0);hi=m.vertices.max(0)
            if np.all(hi>=center_point-radii)&np.all(lo<=center_point+radii):relevant.append(pid)
        conservative=[]
        for pid in relevant:
            if term=='A12.1.04.009'and pid in a['finiteCardiacIntegration']['rootPartIds']:continue
            try:obstacle=solid(pid)
            except ValueError:
                # A pericardial recess is outside the arterial/caval lumen as
                # well as its wall. A native root envelope is conservative;
                # it is never asserted as a newly built arterial wall.
                obstacle=md.Manifold.hull_points(native(pid).vertices);conservative.append(pid)
            roi=roi-obstacle
        if roi.is_empty():raise ValueError('No demonstrated empty region at the scoped anatomical site')
        m=as_mesh(roi);m.save(out/(term+'.npz'));results[term]=dict(type='space-region',file=term+'.npz',vertices=len(m.vertices),triangles=len(m.faces),filled=False,method='bounded empty region excluding actual native cardiac tissue',excludedTissuePartIds=relevant,conservativeOpenRootEnvelopes=conservative)
        print(term,len(m.faces),'unfilled space triangles',round(time.monotonic()-started,2),'s',flush=True)
    except (ValueError,RuntimeError)as e:issues.append(dict(term=term,reason=str(e)))

atomic_json(out/'construction.json',dict(schema='cardiac-finite-reference-geometries/v1',geometries=results,
    required_finite_terms=38,source_catalogue_sha256=digest_file(path),source_bank_sha256=a['finiteCardiacIntegration']['coupledBankSha256'],
    original_tissue_geometry_unchanged=True,original_skeleton_unchanged=True,issues=issues,
    remaining=[t for t in sites if types['required_type'](t)=='finite-mesh'and t not in types['EXISTING_MESHES']and t not in results],
    elapsed_seconds=time.monotonic()-started,installed=False,anatomical_acceptance=False))
print('Finite reference constructions',len(results),'remaining', [t for t in sites if types['required_type'](t)=='finite-mesh'and t not in types['EXISTING_MESHES']and t not in results],'elapsed',round(time.monotonic()-started,2),'s',flush=True)
