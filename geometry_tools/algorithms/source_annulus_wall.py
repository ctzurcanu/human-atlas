"""An explicitly declared missing lateral wall between existing source profiles.

The partial-annulus topology proof is reused; permission to emit this wall is
a separate caller-owned construction decision, never inferred from a hole.
"""
from .partial_annulus_closure import complete_five_quad_collar


def reconstruct_declared_sidewall(vertices,faces,cell_faces,root_ids):
    added,proof=complete_five_quad_collar(vertices,faces,cell_faces,root_ids)
    proof=dict(proof,computationOnly=False,emittedGapFacesPermitted=True,
        representationClass='explicit-source-profile-lateral-wall',
        sourcePointsMoved=0,rootOrEndCap=False,
        acceptance='Caller must verify source ownership, full union, profiles, fields, contacts and rendered review')
    return added,proof
