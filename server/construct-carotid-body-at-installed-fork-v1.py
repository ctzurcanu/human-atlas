"""One source-bound macro-body candidate; no fitting or atlas installation."""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.distance import point_triangle
from geometry_tools.audit import audit,topology
from geometry_tools.common import atomic_json

SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
EXPECTED='a9ae0fbaaf62fdf600a966e43641f54a4b20621944784e31ec43f43260428638'
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1'
O=B/'body-at-installed-fork-v1'
AP=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME='model:local-ta98:assembled-metres'

def winding(triangles,point):
    v=triangles-point;n=np.linalg.norm(v,axis=2)
    numerator=np.einsum('ij,ij->i',v[:,0],np.cross(v[:,1],v[:,2]))
    denominator=n.prod(1)+np.einsum('ij,ij->i',v[:,0],v[:,1])*n[:,2]+np.einsum('ij,ij->i',v[:,1],v[:,2])*n[:,0]+np.einsum('ij,ij->i',v[:,2],v[:,0])*n[:,1]
    return float(np.sum(2*np.arctan2(numerator,denominator))/(4*np.pi))

def main():
    began=time.monotonic();assert not O.exists();O.mkdir()
    assert SHA(AP)==EXPECTED
    E=B/'installed-body-sinus-landmarks-v1/evidence.json';e=json.loads(E.read_bytes())
    assert e['complete']and e['source_catalogue_sha256']==EXPECTED
    C=np.array(e['one_analytic_clearance_constraint']['query_centre_m'])
    Q=np.array(e['body_landmark_frame']['basis_columns_medial_cranial_posterior'])
    R=np.array(e['body_size_constraint']['illustrative_diameters_m'])/2
    assert np.max(abs(Q.T@Q-np.eye(3)))<1e-12 and abs(np.linalg.det(Q)-1)<1e-12
    inventory=B/'complete-authority-inventory.json';inv=json.loads(inventory.read_bytes())
    assert inv['term_count']==222 and any(t['authority']['id']=='A12.2.04.007'and t['representation_class']=='whole-organ'for t in inv['terms'])
    contract=dict(schema='one-installed-fork-carotid-body-contract/v1',source_catalogue_sha256=EXPECTED,
        full_222_term_inventory_sha256=SHA(inventory),landmark_evidence_sha256=SHA(E),
        whole_assembly='Fixed installed left CCA/bifurcation/ICA/ECA and six coupled lingual cells, plus one independently constructed chemoreceptor macro-body.',
        planned_part=dict(term='A12.2.04.007',parent='A12.2.04.006',representation='macroscopic-organ',existing_source_mesh=None,missing_before=True),
        centre_m=C.tolist(),basis_columns_medial_cranial_posterior=Q.tolist(),semiaxes_m=R.tolist(),
        volume_mm3=float(4*np.pi*np.prod(R)/3*1e9),size_qualification='Axial dimensions are reference-informed; 6mm long-axis length and smooth ellipsoid shape are illustrative, not patient-specific.',
        primary_reference='https://www.nature.com/articles/s41598-025-05974-2.pdf',
        adjoining_interfaces=['Adjacent to the actual posteromedial fork, not the historical coincident source roots.',
            'No attachment, glomic arterial mouth, capillary bed or carotid-sinus nerve is constructed or certified by this macro-body.',
            'No histological capsule is added; the closed surface bounds the organ only.'],
        fixed='Every current master record, skeleton, heart, bowel, arterial/lingual wall and female anatomy remains byte-identical.',
        gate='Actual current-forest bounds select all potentially contacting sources. Exact analytic ellipsoid metric excludes every selected active triangle; closed source centre containment is checked separately. Open/invalid source occupancy remains unknown.',
        maximum_stages=1,optimizer=False,wall_seconds=40,installed=False,anatomical_acceptance=False)
    atomic_json(O/'contract-before-geometry.json',contract)
    atlas=json.loads(AP.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    fp=ROOT/'reports/structural-dependencies/current/male-forest.json';forest=json.loads(fp.read_bytes())
    assert forest['sourceSha256']==EXPECTED
    sp=ROOT/forest['spatialIndex']['file'];assert SHA(sp)==forest['spatialIndex']['sha256']
    with np.load(sp,allow_pickle=False)as z:bounds=z['bounds'].copy();ids=z['partIds'].copy()
    assert len(ids)==len(parts)and set(ids)==set(parts)
    distance=np.maximum(np.maximum(bounds[:,0]-C,C-bounds[:,1]),0)
    take=np.linalg.norm(distance,axis=1)<=R.max()
    rows=[];violations=[];hidden=set(atlas.get('defaultHidden',[]))
    for pid in ids[take]:
        if time.monotonic()-began>34:raise TimeoutError('Native body evidence deadline before completion')
        p=parts[str(pid)];m=reader.read_part(atlas,p,frame=FRAME,catalogue_path=AP,catalogue_sha=EXPECTED)
        # A rigid ellipsoid-to-unit-sphere map gives an exact continuous
        # triangle-surface exclusion, rather than checking only vertices.
        t=((m.triangles-C)@Q)/R
        d,_,bary=point_triangle(np.zeros((len(t),3)),t);j=int(np.argmin(d));radius=float(np.sqrt(max(d[j],0)))
        active=not(p.get('suppressed')or p.get('defaultHidden')or pid in hidden)
        top=topology(m);closed=top['closed_oriented_manifold']
        w=winding(m.triangles,C)if closed else None
        envelope=p['system']in ['integumentary','fascia']
        row=dict(id=str(pid),active=active,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,
            minimum_normalized_surface_radius=radius,nearest_face=j,nearest_barycentric=bary[j].tolist(),
            no_surface_crossing=radius>1,conservative_gap_lower_bound_m=max(0,(radius-1)*R.min()),
            closed_oriented_source=closed,centre_winding=w,centre_inside=abs(w)>.5 if w is not None else None,
            envelope_system=envelope,occupied_volume_qualification='Enclosing skin/fascia is not interpreted as solid filled tissue.'if envelope else 'Open/invalid sources have unknown occupancy.'if not closed else 'Centre containment supplements analytic continuous surface exclusion; source self-intersections are not independently certified.')
        if active and(radius<=1 or closed and abs(w)>.5 and not envelope):violations.append(str(pid))
        rows.append(row)
    atomic_json(O/'complete-source-context.json',dict(complete=True,source_catalogue_sha256=EXPECTED,
        forest_sha256=SHA(fp),spatial_index_sha256=SHA(sp),all_master_parts=len(ids),
        broadphase_excluded_parts=int((~take).sum()),actual_narrowphase_parts=len(rows),contexts=rows,
        rejected_active_sources=violations,passed=not violations,geometry_created=False,
        elapsed_seconds=time.monotonic()-began))
    assert not violations,violations
    nlon=24;nrings=10;sphere=[[0,1,0]]
    for j in range(1,nrings+1):
        phi=np.pi*j/(nrings+1)
        for k in range(nlon):
            theta=2*np.pi*k/nlon;sphere.append([np.sin(phi)*np.cos(theta),np.cos(phi),np.sin(phi)*np.sin(theta)])
    sphere.append([0,-1,0]);sphere=np.array(sphere);v=(sphere*R)@Q.T+C;f=[]
    for k in range(nlon):f.append([0,1+k,1+(k+1)%nlon])
    for j in range(nrings-1):
        lo=1+j*nlon;hi=lo+nlon
        for k in range(nlon):
            kk=(k+1)%nlon;f.extend([[lo+k,hi+k,hi+kk],[lo+k,hi+kk,lo+kk]])
    south=len(v)-1;lo=1+(nrings-1)*nlon
    for k in range(nlon):f.append([lo+k,south,lo+(k+1)%nlon])
    f=np.array(f,np.int64);t=v[f];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);wrong=np.einsum('ij,ij->i',n,t.mean(1)-C)<0;f[wrong]=f[wrong][:,[0,2,1]]
    normals=(sphere/R)@Q.T;normals/=np.linalg.norm(normals,axis=1)[:,None]
    src=dict(kind='original-analytic-macroscopic-organ-candidate',term='A12.2.04.007',parent='A12.2.04.006',side='left',
        construction_contract_sha256=SHA(O/'contract-before-geometry.json'),capsule_represented=False,
        geometry_origin_policy='Explicit local bank origin; applied once',clinical_packing_accepted=False)
    origin=np.array([0,1.625,0.]);path=O/'left-carotid-body-macro-candidate.npz'
    np.savez_compressed(path,positions=v-origin,geometryOrigin=origin,faces=f,attribute_normals=normals,frame=np.array(FRAME),source_json=np.array(json.dumps(src,sort_keys=True)))
    m=Mesh.load(path);proof=audit(m,require_closed=True,expected_components=1,expected_genus=[0]);assert proof['passed'],proof['failures']
    assert np.max(abs(m.vertices-v))<1e-15 and SHA(AP)==EXPECTED
    atomic_json(O/'candidate.json',dict(complete=True,candidate_only=True,path=str(path),sha256=SHA(path),
        geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,triangles=len(f),vertices=len(v),
        analytic_volume_mm3=contract['volume_mm3'],mesh_volume_mm3=proof['topology']['oriented_triangle_integral_m3']*1e9,
        source_catalogue_unchanged=True,bytes=path.stat().st_size,images=0,UVs=0,numerical_audit=proof,
        context_sha256=SHA(O/'complete-source-context.json'),anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began))
    print(json.dumps(dict(complete=True,triangles=len(f),vertices=len(v),actual_context_parts=len(rows),active_context_parts=sum(r['active']for r in rows),candidate=str(path),bytes=path.stat().st_size,elapsed_seconds=time.monotonic()-began)))

if __name__=='__main__':main()
