# Model visits

The Visit button sits between Find and Sections. Opening or closing its panel
preserves the camera. Forward and Reverse explicitly start at opposite ends of
a path. Pause leaves the camera at the current position; Resume continues there.
Dragging while playing changes the look direction independently of path position.
Up/Down change speed (0.1–4×); Space pauses/resumes when focus is outside a form
control. Reset, view presets, Center and Explode stop playback.

Around/front/back paths are generated from the current model bounds. They face
the model by default, with an option to face the travel tangent. Local male and
female TA98 models additionally offer the alimentary guide, facing its tangent
by default. This path uses an arc-length sampled centripetal 3D spline through the guide
landmarks. There are no straight-span/rounded-corner transitions.
A 12 mm spatial filter removes guide jitter; a 25 mm camera look-ahead and bounded
turn rate prevent abrupt heading changes during playback. OrbitControls and
covering ray scans are suspended while the Visit camera drives playback. Frames use anterior-up (dorsal-down), retaining the prior up direction at axial
poles. The entry orientation is reached during the exterior approach. After the anal outlet, the camera turns
toward it during the exterior departure. Loop repeats playback when enabled. One traversal
at 1× takes 120 seconds. No extra organ geometry or texture is allocated.

The Visit panel uses the same dimensions as item details. Its position slider
tracks playback and supports seeking. Each alimentary segment has a matching
slider color and list bullet; clicking an organ seeks to its beginning. The
route starts 60 cm anterior to the mouth and ends outside between the anal outlet and toes, looking toward the anus. Opening or closing the panel and choosing a path do not move the camera;
playing or seeking explicitly authorizes path navigation.

## Local alimentary route generation

Run `npm run build:alimentary-visits` after changing installed digestive anatomy.
The Python script requires NumPy, as do the existing anatomy construction tools.
It reads the installed ordinary local TA98 catalogues, packed positions, origins,
source offsets and canonical face/rim ownership. It does not modify them. Its
outputs stay under `.local-models/ta98-review/visits/`; provenance and interface
bridge distances are recorded in `reports/alimentary-visit-paths.json`.
Runtime checks part IDs, packed offsets, counts, chunk URLs and bounds before
making a guide available for the loaded model.

The route covers mouth, pharynx, esophagus, cardia, gastric body, pyloric part,
duodenum, jejunum, ileum, caecum, ascending/transverse/descending/sigmoid colon,
rectum and anal outlet. The appendix is a blind branch and is excluded from the
through route. Small intestine and colon guides average installed innermost
coordinates by surface-geodesic distance between source openings/region ends.
Upper tract and stomach use surface-derived anatomical guides; the male anal
outlet uses the native sphincter landmarks because a complete canal wall is
unrepresented. The mouth guide lies above the tongue surface, not inside it.

**This camera guide does not establish anatomical lumen continuity.** Existing
inter-organ gaps, unresolved source registration, self-crossings and incomplete
wall representations remain unchanged. Interpolated bridges are recorded rather
than being asserted to be actual open lumens. No geometry acceptance is granted
by successful playback. These locally licensed TA98 guides are not published.

Validation: `npm run test:model-visit`, `npm run test:camera-intent`,
`npm run test:explosion-plane`, `npm run test:explosion-viewport`, `npm run check`.

The orbit slider marks four 90-degree regions; orbit framing fits the whole model
inside the usable viewport. Front/back sweeps run continuously from head to feet,
with height-based centimeter ranges for head, torso, pelvis, legs and feet. Those
region ranges are navigation bands, not certified anatomical boundaries. The
active navigation list item scrolls toward the list midpoint on region changes.

Playback glides for 2.4 seconds from the exact requesting camera pose to the
route start, using quintic easing. After a non-looping route completes, it glides
back and restores that pose, including target, orientation, zoom and projection
offset, exactly. Orbit playback also stops after one traversal when Loop is off.
Pausing freezes either glide as well as route playback. Explicit seeks cancel a
pending glide and place the camera at the requested route position.
