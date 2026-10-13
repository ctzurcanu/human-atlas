# Anatomical geometry work

Before splitting, replacing, or reconstructing an anatomical object:

1. Inventory every applicable TA98 descendant of the whole organ against the authoritative terminology. Record all terms, their parents, existing source meshes, and missing representations in a checked-in report.
2. Distinguish tissue regions, wall layers, openings, boundary curves, and microscopic structures. An opening is not a separate solid. A generic outer surface does not establish a wall layer or microscopic feature.
3. Define the whole-organ assembly, each planned part, and every adjoining interface before modifying geometry. Include cardia when dissecting the stomach. Do not substitute one coarse source partition for the full terminology.
4. Repair the complete surface before partitioning it. Inspect the actual rendered geometry and cross-sections. Zero open edges alone is insufficient: check connected components and unintended tunnels, as well as self-intersections, lumen continuity, and shared coordinates/normals.
5. Partition shared tissue without duplicating overlapping faces. Keep missing geometry and provisional anatomical boundaries explicitly marked as such. Do not call a repair complete while visible defects remain.

Local review models have unverified redistribution rights. Keep them local unless their licenses are established and deployment is authorized.

Follow `ANATOMICAL-GEOMETRY-MANUAL.md` for organ packing, geometry acceptance and runtime efficiency. Run the applicable evidence scripts before accepting a new assembly. A script pass does not replace reference, render, cross-section or anatomical review.

# Exploded presentation invariant

Every exploded group must occupy the same shared layout plane, perpendicular to the current point of view. Orbit, camera roll, saved views and section-camera changes must not tilt that plane relative to the viewer. Preserve the relative anatomy within each group. Do not replace this rule while changing geometry, cameras, or models. Run `npm run test:explosion-plane` for relevant changes.

# Camera intent invariant

For every model, selecting or clearing anatomy, opening panels, showing components or relations, changing visibility, isolation, hierarchy, sections or depth must preserve camera position, orientation, orbit target and zoom. Only explicit camera gestures, Center, view presets, saved-view navigation, rotation or Reset authorize movement. Explicit Explode/Implode commands also authorize fitting while their animation runs: all moving anatomy must remain inside the usable viewport. Explode begins at the exact current view and fits only the current moving extent. Full Implode restores the exact pre-explosion position, target, orientation, zoom and projection offset; this takes precedence over centering or fitting a deliberately cropped original view. Explicit camera gestures or navigation establish a new view. Selection must never activate fitting or centering. The explosion plane remains perpendicular to POV. Run the camera-intent and explosion-plane tests for related changes.

The first frame of each Explode/Implode command preserves both the current camera and the current anatomical arrangement, including commands that reverse an ongoing animation. Never fit a future explosion stage or reset the view at this boundary. Run `npm run test:explosion-viewport` to verify the exact start and return as well as viewport containment.
