# Mesh audit

`npm run audit:meshes` regenerates one JSON report per model. Import commands run
the audit for their own output. Each report lists mesh boundary edges, edges
shared by more than two faces, inconsistent triangle winding, duplicate and
degenerate faces, source normals opposed to triangle winding, missing UVs for
mapped materials, exact matching bounds, and parts extending beyond a single
whole-body skin mesh when one exists.

These are **diagnostic flags**, not automatic anatomical corrections. Open
fascia, skin patches, and hollow structures legitimately have boundary edges.
The bounds test is only a placement warning. A texture map without UVs is a
build failure because the image cannot be sampled correctly.

The viewer also checks the intersected mesh topology at section time. It draws
a narrow cut edge for substantial open sheets and a filled cap for solid
contours. Hollow organs keep an open lumen.
