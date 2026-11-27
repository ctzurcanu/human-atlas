# Uterine vascular repair inventory

Authority: `.local-models/ta98-authority.sqlite`, SHA256 `b4c41fa4fe13f68e2735bc59c94a744d58957e1e0dcf4e2172415a400f61ff48`.
All descendants of uterine artery, uterine veins and uterine venous plexus are listed before changing geometry. No vascular geometry is accepted.

| TA98 | Parent | Latin / English | Existing source meshes | Missing representation |
| --- | --- | --- | --- | --- |
| A12.2.15.029F | A12.2.15.001 | a. uterina / uterine artery | LOCAL:female:l_uterine_artery:1704, LOCAL:female:r_uterine_artery:1702 | False |
| A12.2.15.030F | A12.2.15.029F | rr. helicini / helicine branches | None | True |
| A12.2.15.031F | A12.2.15.029F | rr. vaginales / vaginal branches | None | True |
| A12.2.15.032F | A12.2.15.031F | a. azygos vaginae / azygos artery of vagina | None | True |
| A12.2.15.033F | A12.2.15.029F | r. ovaricus / ovarian branches | None | True |
| A12.2.15.034F | A12.2.15.029F | r. tubarius / tubal branch | None | True |
| A12.3.10.015F | A12.3.10.004 | vv. uterinae / uterine veins | LOCAL:female:l_uterine_vein:1770, LOCAL:female:r_uterine_vein:1772 | False |
| A12.3.10.016F | A12.3.10.004 | plexus venosus uterinus / uterine venous plexus | LOCAL:female:uterine_venous_plexus:1778, LOCAL:female:uterine_venous_plexus:1779 | False |

The two named uterine vein trunks provide partial parent coverage of A12.3.10.015F. This exposes their existing meshes in the hierarchy; it does not establish every tributary or accept their geometry.

The unassigned local uterine arterial plexus is an additional source concept; exact coordinate joins to named uterine arteries do not automatically establish its TA98 identity.

Representation classes: the source objects are polygonal vessel/plexus exteriors. Open terminal loops are openings, not separate solids. Boundary rings are interface constraints. No source outer sheet alone establishes intima, media, adventitia, vasa vasorum or microscopic endothelium.

Assembly obligation: preserve the six source interfaces and 21 paired coordinates across the three plexus meshes and four uterine artery/vein trunks. All source topology defects, branch/end identities and vascular supply relationships require review; do not move a whole plexus independently or indiscriminately cap all open ends.

These female-specific roots have no corresponding uterine anatomy in the male model. The male pelvic vascular anatomy needs its own authoritative inventory.

## Source seam and hierarchy result

The two native uterine venous plexus meshes remain 113 and 46 disconnected coordinate-welded components. Their nearest cross-component gaps are 0.3543–6.2708 mm and 0.5447–9.1133 mm; none is within one Float32 coordinate ULP. A hypothetical 0.1 micrometre coordinate union joins no components. Thus a rounding-scale weld is not a justified repair. No coordinates or topology were changed. Source-bound diagnostic: `reports/ta98-female/uterine-venous-component-seam-probe.json`. Full branch/end interpretation remains before reconstruction.

Existing named left/right native uterine vein trunks now have reviewed, source-scoped partial parent mappings to A12.3.10.015F in `scripts/ta98-reviewed-mappings.json`. Metadata and the public TA98 hierarchy are regenerated; the uterine-veins selector resolves these existing source meshes in both local female and local female TA98 catalogues. This fixes a hierarchy omission without asserting complete uterine venous anatomy or accepting the source geometry. The authority inventory contains eight terms; five arterial branch terms still lack separate representations.

The hierarchy validator now checks the established colon-wall plan nodes as well as the small-intestine plan, instead of its obsolete pre-colon node total. It verifies all 7,474 canonical term selectors and full descendant mesh sets across eight installed catalogue contexts, plus explicit partial-parent coverage of the two uterine vein trunks. The camera-intent test passes; no viewer camera/geometry code changed. Native source terminals, all six vascular source joins, clinical flow, plexus defects and the candidate wall's surrounding-structure acceptance remain unresolved.
