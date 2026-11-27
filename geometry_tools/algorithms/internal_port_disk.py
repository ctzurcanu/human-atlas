"""Remove a certified blocking disk between two opposite source tube collars.

No coordinates, attributes, external mouths or surrounding faces change.
Terminal caps, partial leaflets and same-side branches are not inferred joins.
"""
import numpy as np
from geometry_tools.model import Mesh
from .three_way_source_junction import triple_ring_parts
from .source_root_sections import recover_root_sections
from .tube_continuity import boundary_loops, continuity_audit


def remove_internal_port_disk(mesh):
    root, groups = triple_ring_parts(mesh)
    vertices, representatives, physical = np.unique(
        mesh.vertices, axis=0, return_index=True, return_inverse=True)
    root_ids = physical[root]
    disk_groups = [ids for ids in groups
                   if set(np.unique(physical[mesh.faces[ids]])) == set(root_ids)]
    if len(disk_groups) != 1:
        raise ValueError('Exactly one complete ring-only source disk is required')
    disk = disk_groups[0]
    audit = continuity_audit(mesh.vertices, mesh.faces[disk],
                             expected_mouths=1, expected_components=1, expected_euler=1)
    if not audit['topological_continuity_pass']:
        raise ValueError('Source patch is not one complete disk with the shared ring as boundary')
    center = mesh.vertices[root].mean(0)
    _, _, basis = np.linalg.svd(mesh.vertices[root]-center, full_matrices=False)
    normal = basis[-1]
    plane_error = float(np.abs((mesh.vertices[root]-center) @ normal).max())
    sides, collars, ports = [], [], []
    for group in groups:
        if np.array_equal(group, disk):
            continue
        sections = recover_root_sections(mesh.vertices, mesh.faces[group], root, maximum_layers=1)
        if not sections:
            raise ValueError('Both opposite source collars must be independently certified')
        collar = group[np.asarray(sections[0]['sourceFaceIds'])]
        ring = np.asarray(sections[0]['rawIds'])
        outside = np.setdiff1d(np.unique(physical[mesh.faces[collar]]), root_ids)
        signed = (vertices[outside]-center) @ normal
        margin = plane_error + 1e-8
        side = 1 if np.all(signed > margin) else -1 if np.all(signed < -margin) else 0
        if side == 0:
            raise ValueError('A collar does not lie strictly beyond the shared source disk')
        sides.append(side)
        collars.append(collar)
        ports.append(ring)
    if sorted(sides) != [-1, 1]:
        raise ValueError('Source collars must lie on opposite sides; no same-side weld is inferred')
    retained = np.setdiff1d(np.arange(len(mesh.faces)), disk)
    candidate = Mesh(mesh.vertices.copy(), mesh.faces[retained].copy(), mesh.frame,
                     {k: v.copy() for k, v in mesh.attributes.items()}, mesh.source)
    carrier_ids = np.concatenate(collars)
    carrier_audit = continuity_audit(mesh.vertices, mesh.faces[carrier_ids],
                                    expected_mouths=2, expected_components=1, expected_euler=0)
    if not carrier_audit['topological_continuity_pass']:
        raise ValueError('Removing the disk does not expose one continuous two-port source carrier')
    proof = dict(complete=True, method='remove-internal-port-disk',
                 sourceSharedRawIds=root.tolist(), removedSourceFaceIds=disk.tolist(),
                 retainedSourceFaceIds=retained.tolist(), sourceCollarFaceIds=carrier_ids.tolist(),
                 sourcePortRawIds=[p.tolist() for p in ports], sourceRingPlaneError_m=plane_error,
                 oppositeSourceCollars=True, diskTopology=audit, carrierAudit=carrier_audit,
                 newFaceCount=0, newVertexCount=0, newFacesStart=len(retained),
                 newFaceSourceLineage=[], computationRootCapsRetained=0, emittedLumenDisks=0,
                 sourceControlsAttributesRetained=True, finiteVesselWallEstablished=False,
                 fullSurfacePackingPatencyAndAnatomicalAcceptancePending=True)
    return candidate, proof
