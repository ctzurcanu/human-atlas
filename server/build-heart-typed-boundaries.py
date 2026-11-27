"""Native facets, transported surface curves and literal anatomical mouth rims.

This stage constructs no tissue and never spans an opening with a filled disk.
It corrects the dimension of the reference independently of later tissue jobs.
"""
from pathlib import Path
import argparse,copy,json,runpy,sys,time
import numpy as np
import igl
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra,connected_components
from scipy.spatial import ConvexHull,cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.geometry_instance import instance_pose
types=runpy.run_path(str(ROOT/'server/heart-reference-types.py'))
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--volume-stage',default='typed-final-buffers-v665');args=ap.parse_args()
started=time.monotonic();out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
atlas=json.loads(path.read_text());bundle=json.loads(path.with_name('cardiac-label-sites.json').read_text())
parts={p['id']:p for p in atlas['parts']};groups={g['key']:g for g in atlas['finiteCardiacIntegration']['termGroups']}
sites={g['term']:g for g in bundle['integration']['termGroups']};reader=AtlasReader();meshes={};graphs={};welds={}
ep=parts['TA98HEART:male:epicardium'];R,scale=instance_pose(ep);origin=np.asarray(ep['geometryOrigin'])
local=lambda q:((np.asarray(q)-origin)@R/scale).tolist()
def mesh(pid):
    if pid not in meshes:meshes[pid]=reader.read_part(atlas,parts[pid],frame='assembled-metres',catalogue_path=path)
    return meshes[pid]
def weld(pid):
    if pid not in welds:welds[pid]=np.unique(mesh(pid).vertices,axis=0,return_index=True,return_inverse=True)
    return welds[pid]
def chart(pid,peer=None,free=False):
    ds=parts[pid].get('tissueDrawGroups',[])
    selected=[np.arange(d['start']//3,(d['start']+d['count'])//3)for d in ds if
        (not d.get('neighbor')if free else d.get('neighbor')==peer)]
    return np.concatenate(selected)if selected else np.array([],dtype=int)if ds else np.arange(len(mesh(pid).faces))
def projected(pid,points,faces=None):
    m=mesh(pid);f=m.faces if faces is None else m.faces[faces]
    ds,ii,pp=igl.point_mesh_squared_distance(np.asarray(points,float).reshape(-1,3),m.vertices,f)
    return pp,np.asarray(ii if faces is None else faces[ii]),np.sqrt(ds)
def surface_path(pid,points,faces=None):
    m=mesh(pid);v,first,remap=weld(pid);selected=remap[m.faces if faces is None else m.faces[faces]]
    key=(pid,None if faces is None else digest_file(path)+str(len(faces))+str(int(faces[0])))
    if key not in graphs:
        edges=np.unique(np.sort(np.concatenate([selected[:,[0,1]],selected[:,[1,2]],selected[:,[2,0]]]),axis=1),axis=0)
        weights=np.linalg.norm(v[edges[:,0]]-v[edges[:,1]],axis=1)
        graphs[key]=coo_matrix((np.r_[weights,weights],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(v),)*2).tocsr()
    pp,ii,_=projected(pid,points,faces);vertices=[int(remap[m.faces[f][np.argmin(np.linalg.norm(m.vertices[m.faces[f]]-p,axis=1))]])for p,f in zip(pp,ii)]
    path_ids=[]
    for a,b in zip(vertices[:-1],vertices[1:]):
        dist,pre=dijkstra(graphs[key],directed=False,indices=a,return_predecessors=True)
        if not np.isfinite(dist[b]):raise ValueError(f'Disconnected surface route: {pid}')
        route=[b]
        while route[-1]!=a:route.append(int(pre[route[-1]]))
        path_ids.extend(route[::-1][:-1])
    path_ids.append(vertices[-1]);return v[path_ids],first[path_ids].tolist()
def loops(faces):
    edges=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]);keys,inv,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    boundary=keys[counts==1];adj={}
    for a,b in boundary:adj.setdefault(int(a),[]).append(int(b));adj.setdefault(int(b),[]).append(int(a))
    if any(len(ns)!=2 for ns in adj.values()):return []
    unused=set(adj);result=[]
    while unused:
        start=min(unused);route=[start];previous=None;here=start
        while True:
            nxt=next(n for n in adj[here]if n!=previous)
            if nxt==start:break
            route.append(nxt);previous,here=here,nxt
            if len(route)>len(adj):raise ValueError('Invalid boundary traversal')
        unused.difference_update(route);result.append(route)
    return result
def area(q):
    center=q.mean(0);_,_,basis=np.linalg.svd(q-center,full_matrices=False);uv=(q-center)@basis[:2].T
    return abs(float(np.sum(uv[:,0]*np.roll(uv[:,1],-1)-uv[:,1]*np.roll(uv[:,0],-1)))/2)
def rim_for(pid,peer):
    m=mesh(pid);v,first,remap=weld(pid);faces=chart(pid,peer);rings=loops(remap[m.faces[faces]])
    if not rings:return None
    # Inner/outer annular boundaries are nested; select the smaller enclosed
    # contour, not a disk, nor the remote distal end of a vessel.
    ring=min(rings,key=lambda r:area(v[r]));return dict(points=local(v[ring]),closed=True,sourcePartId=pid,sourceVertexIds=first[ring].tolist(),areaMetres2=area(v[ring]))
def sectional_rim(pid,target):
    m=mesh(pid);v,first,remap=weld(pid);ff=remap[m.faces];c=v.mean(0);_,_,axes=np.linalg.svd(v-c,full_matrices=False);normal=axes[0];projection=(v-c)@normal;target_t=(np.asarray(target)-c)@normal;candidates=[]
    # Every contour comes from literal triangle/plane intersections. Choose a
    # complete annular station close to the recipient; never synthesize a disk.
    stations=np.unique(np.clip(target_t+np.linspace(-.002,.002,9),projection.min()+.0001,projection.max()-.0001))
    for station in stations:
        d=(v-c)@normal-station;segments=[]
        for face in ff[np.min(d[ff],axis=1)*np.max(d[ff],axis=1)<0]:
            pp=[]
            for aa,bb in zip(face,np.roll(face,-1)):
                if d[aa]*d[bb]<0:pp.append(v[aa]+d[aa]/(d[aa]-d[bb])*(v[bb]-v[aa]))
            if len(pp)==2:segments.append(pp)
        if not segments:continue
        p=np.asarray(segments).reshape(-1,3);keys,first_p,inv=np.unique(np.rint(p/1e-9).astype('int64'),axis=0,return_index=True,return_inverse=True);points=p[first_p];edges=inv.reshape(-1,2);adj={}
        for aa,bb in edges:adj.setdefault(int(aa),[]).append(int(bb));adj.setdefault(int(bb),[]).append(int(aa))
        if any(len(ns)!=2 for ns in adj.values()):continue
        unused=set(adj);rings=[]
        while unused:
            start=min(unused);route=[start];prev=None;here=start
            while True:
                nxt=next(x for x in adj[here]if x!=prev)
                if nxt==start:break
                route.append(nxt);prev,here=here,nxt
                if len(route)>len(adj):break
            unused.difference_update(route)
            if len(route)>2:rings.append(points[route])
        if len(rings)<2:continue
        ring=min(rings,key=area);candidates.append((abs(station-target_t),dict(points=local(ring),closed=True,sourcePartId=pid,areaMetres2=area(ring),sectionPlanePoint=(c+station*normal).tolist(),sectionPlaneNormal=normal.tolist(),recipientStationDistanceMetres=abs(station-target_t))))
    return min(candidates,key=lambda x:x[0])[1]if candidates else None
def footprint(pid,q,radius,faces=None,normal=None):
    m=mesh(pid);faces=np.arange(len(m.faces))if faces is None else faces
    pp,fi,_=projected(pid,[q],faces);tri=m.triangles[faces];centers=tri.mean(1)
    keep=np.linalg.norm(centers-pp[0],axis=1)<=radius
    if normal is not None:
        n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);keep&=(n@normal)>0
    if not np.any(keep):keep[np.argmin(np.linalg.norm(centers-pp[0],axis=1))]=True
    # Boolean tissue charts can have coincident T-junctions. Their topological
    # facet connectivity must not reduce a material surface to one tiny face.
    # Keep the spatial footprint on one anatomically facing side instead.
    picked=faces[keep]
    t=m.triangles[picked];size=float(np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1).sum()/2)
    assert size>0
    return dict(partId=pid,faces=picked.tolist(),areaMetres2=size,sourceBankSha256=atlas['finiteCardiacIntegration']['coupledBankSha256'])

result={};issues=[];endo='TA98HEART:male:endocardium';ias='TA98HEART:male:myocardium:interatrial-septum';ra='TA98HEART:male:myocardium:right-atrium'
for term in sorted(types['CURVES']):
    g=sites[term];pid=g['anchor']['partId'];points=g.get('points')
    if term=='A12.1.01.006':continue
    if not points or len(points)<2:
        q=np.asarray(g['anchor']['point']);m=mesh(pid);v=m.vertices;_,_,axes=np.linalg.svd(v-v.mean(0),full_matrices=False)
        points=[q-axes[0]*.003,q,q+axes[0]*.003]
    if term=='A12.1.00.008':
        q=np.asarray(g['anchor']['point']);points=[q+[-.001,.002,0],q,q+[.001,.002,0]]
    if term=='A12.1.00.006':
        m=mesh(pid);outer_ids=np.unique(m.faces[chart(pid,'TA98HEART:male:epicardium')]);v=m.vertices[outer_ids];points=[]
        for yy in np.linspace(np.quantile(v[:,1],.08),np.quantile(v[:,1],.92),9):
            band=v[np.abs(v[:,1]-yy)<.004];points.append(band[np.argmin(band[:,0])])
    try:
        coords,vids=surface_path(pid,points,chart(pid,free=True)if pid.endswith(':epicardium')else None)
        result[term]=dict(type='boundary-curve',paths=[dict(points=local(coords),closed=False,partId=pid,sourceVertexIds=vids)],sourceBound=True)
    except ValueError as e:issues.append(dict(term=term,reason=str(e)))

# Actual finite septal/endocardial facets replace point markers of the fossa,
# smooth caval sinus and Koch triangle. No extra coincident tissue is created.
for term,radius in [('A12.1.01.005',.005),('A12.1.01.007',.005),('A12.1.00.030',.006),('A12.1.01.012',.008)]:
    g=sites[term];pid=endo
    facing=mesh(ra).vertices.mean(0)-np.asarray(g['anchor']['point']);facing/=np.linalg.norm(facing);patch=footprint(pid,g['anchor']['point'],radius,chart(pid,free=True),normal=facing)
    result[term]=dict(type='surface-region',surfaces=[patch],sourceBound=True,
                      state='closed-adult-septal-site'if term=='A12.1.01.007'else'present',boundaryProvisional=True)
    if term=='A12.1.00.030':
        # A triangle is a bounded region, never a circular point halo. Its
        # three landmarks are the CS ostium, septal tricuspid hinge and the
        # Todaro/central-fibrous apex, all projected onto native RA facets.
        av_node=Mesh.load(ROOT/'.local-models/geometry-workbench/heart'/args.volume_stage/'A12.1.06.004.npz').vertices.mean(0)@R.T*scale+origin;cs=np.array(sites['A12.1.01.009']['anchor']['point']);septal=next(pid for pid in groups['tricuspid-valve']['partIds']if 'septal'in parts[pid]['name'].lower());vv=mesh(septal).vertices;center=np.concatenate([mesh(pid).vertices for pid in groups['tricuspid-valve']['partIds']]).mean(0);distance=np.linalg.norm(vv-center,axis=1);hinges=vv[distance>=np.quantile(distance,.7)];axis=av_node-cs;axis/=np.linalg.norm(axis);perp=(hinges-cs)-((hinges-cs)@axis)[:,None]*axis;hinge=hinges[np.argmax(np.linalg.norm(perp,axis=1))];witnesses=[cs,hinge,av_node];charts={};proposals=[]
        for owner in [endo]:
            ff=chart(owner,free=True);tt=mesh(owner).triangles[ff];nn=np.cross(tt[:,1]-tt[:,0],tt[:,2]-tt[:,0]);ff=ff[np.einsum('ij,ij->i',nn,mesh(ra).vertices.mean(0)-tt.mean(1))>0]
            if len(ff):charts[owner]=ff;proposals.append(projected(owner,witnesses,ff)[0])
        proposals=np.asarray(proposals);corners=np.array([proposals[np.argmin(np.linalg.norm(proposals[:,i]-point,axis=1)),i]for i,point in enumerate(witnesses)]);axes=np.linalg.svd(corners-corners.mean(0),full_matrices=False)[2];basis=axes[:2];uv=(corners-corners.mean(0))@basis.T;surfaces=[]
        for owner in [endo]:
            ff=chart(owner,free=True);triangles=mesh(owner).triangles[ff];centers=triangles.mean(1);relative=centers-corners.mean(0);cc=relative@basis.T;bary=np.linalg.solve(np.c_[uv[1]-uv[0],uv[2]-uv[0]],(cc-uv[0]).T).T;n=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);towards=mesh(ra).vertices.mean(0)-centers;chosen=ff[(bary[:,0]>=0)&(bary[:,1]>=0)&(bary.sum(1)<=1)&(np.abs(relative@axes[2])<=.005)&(np.einsum('ij,ij->i',n,towards)>0)]
            if len(chosen):
                tt=mesh(owner).triangles[chosen];surfaces.append(dict(partId=owner,faces=chosen.tolist(),areaMetres2=float(np.linalg.norm(np.cross(tt[:,1]-tt[:,0],tt[:,2]-tt[:,0]),axis=1).sum()/2),sourceBankSha256=atlas['finiteCardiacIntegration']['coupledBankSha256']))
        assert surfaces,'Koch region has no native atrial/septal facets'
        result[term].update(surfaces=surfaces,boundaryLandmarks=local(corners))
fossa=result['A12.1.01.005']['surfaces'][0];tt=mesh(endo).triangles[fossa['faces']];center=tt.mean((0,1));normal=mesh(ra).vertices.mean(0)-center;normal/=np.linalg.norm(normal);axis=np.cross(normal,[1,0,0]if abs(normal[0])<.8 else[0,1,0]);axis/=np.linalg.norm(axis);other=np.cross(normal,axis);theta=np.linspace(0,2*np.pi,64,endpoint=False);witness=center+.005*(np.cos(theta)[:,None]*axis+np.sin(theta)[:,None]*other);ff=chart(endo,free=True);triangles=mesh(endo).triangles[ff];n=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);ff=ff[n@normal>0];rim,face,_=projected(endo,witness,ff);result['A12.1.01.006']=dict(type='boundary-curve',paths=[dict(points=local(rim),closed=True,partId=endo,sourceFaceIds=face.tolist())],sourceBound=True,boundaryProvisional=True)

for term,key in [('A12.1.02.015','pulmonary-valve'),('A12.1.04.017','aortic-valve')]:
    surfaces=[]
    valve_center=np.concatenate([mesh(pid).vertices for pid in groups[key]['partIds']]).mean(0)
    for pid in groups[key]['partIds']:
        m=mesh(pid);v=m.vertices;c=v.mean(0);basis=np.linalg.svd(v-c,full_matrices=False)[2];u=(v-c)@basis[0];free=[]
        for low,high in zip(np.linspace(u.min(),u.max(),13)[:-1],np.linspace(u.min(),u.max(),13)[1:]):
            band=v[(u>=low)&(u<=high)]
            if len(band):free.append(band[np.argmin(np.linalg.norm(band-valve_center,axis=1))])
        centers=m.triangles.mean(1);distance=np.min(np.linalg.norm(centers[:,None]-np.asarray(free)[None],axis=2),axis=1);chosen=np.flatnonzero(distance<=.0012);tri=m.triangles[chosen];surfaces.append(dict(partId=pid,faces=chosen.tolist(),areaMetres2=float(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1).sum()/2),sourceBankSha256=atlas['finiteCardiacIntegration']['coupledBankSha256']))
    result[term]=dict(type='surface-region',surfaces=surfaces,sourceBound=True,boundaryProvisional=True)

# Literal root/tissue interfaces supply the caval and four pulmonary rims.
for term,pids,peer in [
    ('A12.1.01.011',['ZA:Superior vena cava'],endo),
    ('A12.1.01.010',['TA98HEPATIC:male:region:A12.3.09.001'],ra),
    ('A12.1.03.004',[pid for pid in atlas['finiteCardiacIntegration']['rootPartIds']if 'pulmonary vein'in parts[pid]['name'].lower()],endo),
]:
    rims=[rim_for(pid,peer)for pid in pids]
    if not all(rims):rims=[r or sectional_rim(pid,sites[term]['anchor']['point'])for pid,r in zip(pids,rims)]
    if all(rims):result[term]=dict(type='opening-rim',paths=rims,sourceBound=True,filled=False)
    else:issues.append(dict(term=term,reason='Recipient interface has no simple native mouth ring'))

cs='ZA:Coronary sinus'
if cs in parts:
    v,first,remap=weld(cs);rings=loops(remap[mesh(cs).faces])
    if rings:
        target=np.asarray(sites['A12.1.01.009']['anchor']['point']);ring=min(rings,key=lambda ids:np.linalg.norm(v[ids].mean(0)-target));result['A12.1.01.009']=dict(type='opening-rim',paths=[dict(points=local(v[ring]),closed=True,sourcePartId=cs,sourceVertexIds=first[ring].tolist(),areaMetres2=area(v[ring]))],sourceBound=True,filled=False,recipientJoinVerified=False)

# A representative smallest-vein foraminal set is a boundary annotation on a
# documented endocardial distribution, not invented filled plugs or a claim
# that individual microscopic channels were present in this source specimen.
term='A12.1.01.004';m=mesh(ra);ff=chart(ra,endo);target=np.asarray(sites[term]['anchor']['point']);pp,face,_=projected(ra,[target],ff);tri=m.triangles[face[0]];normal=np.cross(tri[1]-tri[0],tri[2]-tri[0]);normal/=np.linalg.norm(normal);axis=np.cross(normal,[1,0,0]if abs(normal[0])<.8 else[0,1,0]);axis/=np.linalg.norm(axis);other=np.cross(normal,axis);paths=[]
for offset in [-.001,0,.001]:
    center=projected(ra,[pp[0]+offset*axis],ff)[0][0];theta=np.linspace(0,2*np.pi,16,endpoint=False);ring=center+.0002*(np.cos(theta)[:,None]*axis+np.sin(theta)[:,None]*other);ring=projected(ra,ring,ff)[0];paths.append(dict(points=local(ring),closed=True))
result[term]=dict(type='opening-rim',paths=paths,sourceBound=False,filled=False,representativeMicroscopicDistribution=True,individualSourceChannelsVerified=False)

# Valve ostia and commissures use the actual native leaflet apparatus. The
# hinge hull bounds an empty orifice and does not add an opening disk.
for key,opening,commissures in [('tricuspid-valve','A12.1.02.002',None),('mitral-valve','A12.1.04.002',None),('pulmonary-valve','A12.1.02.009','A12.1.02.016'),('aortic-valve','A12.1.04.010','A12.1.04.018')]:
    pids=groups[key]['partIds'];allv=np.concatenate([mesh(pid).vertices for pid in pids]);center=allv.mean(0);_,_,basis=np.linalg.svd(allv-center,full_matrices=False);uv=(allv-center)@basis[:2].T;hull=ConvexHull(uv);rim=allv[hull.vertices]
    result[opening]=dict(type='opening-rim',paths=[dict(points=local(rim),closed=True,sourcePartIds=pids,areaMetres2=float(hull.volume))],sourceBound=True,filled=False,boundaryProvisional=True)
    if commissures:
        points=[];gaps=[]
        for i in range(3):
            a=mesh(pids[i]).vertices;b=mesh(pids[(i+1)%3]).vertices;radius=np.linalg.norm((a-center)@basis[:2].T,axis=1);a=a[radius>=np.quantile(radius,.7)];dist,idx=cKDTree(b).query(a);j=int(np.argmin(dist));points.append((a[j]+b[idx[j]])/2);gaps.append(float(dist[j]))
        result[commissures]=dict(type='point-set',points=local(points),sourceBound=True,sourcePartIds=pids,attachmentGapMetres=gaps)

for term in sorted(types['EXISTING_MESHES']):result[term]=dict(type='finite-mesh',partIds=sites[term]['partIds'],sourceBound=True,existingGeometry=True)
for term,children in types['COLLECTIONS'].items():result[term]=dict(type='collection',terms=children,sourceBound=True)
atomic_json(out/'typed-geometries.json',dict(schema='cardiac-native-typed-geometries/v1',geometries=result,
    source_catalogue_sha256=digest_file(path),source_bank_sha256=atlas['finiteCardiacIntegration']['coupledBankSha256'],
    frame='heart-native-local-metres',issues=issues,constructed_terms=len(result),new_tissue_meshes=0,
    elapsed_seconds=time.monotonic()-started,anatomical_acceptance=False))
print('Native typed definitions',len(result),'unresolved',issues,'elapsed',round(time.monotonic()-started,2),'s',flush=True)
