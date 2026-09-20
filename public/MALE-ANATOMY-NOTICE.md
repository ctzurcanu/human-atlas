# Anatomy Atlas — credits and licenses

Viewer created by **Brian Pridgen, MD**. Original viewer software is governed by
[the Viewer Software License](LICENSE). You may use the authorized hosted app
and share its link. Its original code is not offered under an open-source license.
The anatomy and third-party components below have separate permissions.

## Anatomical geometry and anatomy-derived metadata

This free educational release includes separately identifiable components with
**different licenses**. It must not be described as uniformly CC BY-SA.

- The previously released adapted anatomy and catalogue remain under
  [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), including
  Brian Pridgen's copyrightable contributions to those adaptations.
- Both kidneys by **Lissie Cowley** are under
  [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).
- Both cochleae and vestibules from **University of Dundee School of Medicine**
  are under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
  Shared adaptations of these components retain that license.

The kidney and inner-ear components may be used only for noncommercial purposes
unless the relevant rights holders grant additional permission. Their inclusion
here does not grant permission to use them in commercial products or promotional
material. Other components retain their own permissions, including commercial
reuse where their licenses permit it. Viewer software terms do not restrict the
rights granted by these component licenses.

### Open 3D Model

Open 3D Model — departments of anatomy of Leiden University Medical Center,
University Medical Center Utrecht, Maastricht University Medical Center and
University of Leuven KULAK, with Amsterdam UMC, Radboud University Medical Center
Nijmegen and University of Ghent. Funded by the Dutch Ministry of Education,
Culture and Science.

- [Upper limb, English labels](https://anatomytool.org/content/open3dmodel-upper-limb-english-labels)
- [Brachial plexus and branches, English labels](https://anatomytool.org/content/open3dmodel-brachial-plexus-and-branches-english-labels)
- Reviewed input edition: July 2025; CC BY-SA 4.0.

### Z-Anatomy and upstream contributors

[Z-Anatomy — The libre 3D atlas of anatomy](https://github.com/Z-Anatomy/Models-of-human-anatomy),
CC BY-SA 4.0. The reviewed source notice is available at
[the pinned upstream license](https://github.com/Z-Anatomy/Models-of-human-anatomy/blob/bb293be456f6d2245191f37752d14ec76894740d/License.txt).
Its upstream attributions are retained here:

- **BodyParts3D — The Database Center for Life Science — CC-BY-SA 2.1 Japan**, as
  credited by Z-Anatomy. [Original model](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html);
  [license](https://creativecommons.org/licenses/by-sa/2.1/jp/).
  Kousaku OKUBO is credited for the original BodyParts3D model.
- Gauthier KERVYN: design, 3D and anatomy.
- “Brainder” and “White matter” from the University of Washington, identified
  by Z-Anatomy as reference/included and adapted material. Its notice does not
  specify separate license terms for these entries; this is not a claim of
  independent clearance of those underlying sources.
- Cranial Nerves and Foramina — University of Dundee, CAHID —
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- The upstream project also credits Marcin ZIELINSKI (Blender add-on), Lluis
  VINENT (Unity development), Ana Teresa BIGIO (Portuguese translation), and
  Carlos TORRES VILLAR (Spanish translation). Their application software is
  not bundled in this viewer.

BodyParts3D: Mitsuhashi N, Fujieda K, Tamura T, Kawamoto S, Takagi T, Okubo K.
BodyParts3D: 3D structure database for anatomical concepts. Nucleic Acids Res
2009;37:D782–5. [PMID 18835852](https://pubmed.ncbi.nlm.nih.gov/18835852/).
The upstream BodyParts3D credit above concerns model ancestry; it does not
assert that a historical measurement dataset is included.

### Changes made for this viewer

Adaptation and integration by Brian Pridgen. Released adaptations include
combining the whole-body and upper-limb sources, mirroring the grafted upper
limb for the opposite side, preparing regional geometry chunks and a derived
body surface, and adapting names, catalogue metadata and display organization.
The viewer also crops regional skin and adjusts display visibility at runtime.
These are adaptations, not unmodified upstream releases. Source and mirroring
information are recorded per structure. Neither a mirrored structure nor a
display layer establishes independent anatomical validation.

### Educational components restored in this release

The six kidney and inner-ear meshes were recovered from the recorded Z-Anatomy
source edition, preserving their source placement. Export changes comprise
applying source object transforms/modifiers, coordinate conversion to glTF,
and assigning separate named objects and license metadata. They are not fused
with other meshes. Existing released mesh geometry was retained unchanged.
Model conversion and integration by Brian Pridgen; the component licenses above
continue to apply. This restoration does not establish anatomical validation.

### Material not included

Historical BodyParts3D-derived calibre measurements
and literature measurement tables are not bundled. Z-Anatomy's Wikipedia
text definitions are not supplied as a definitions dataset in this release.

## Rendering libraries and fonts

- [three.js](https://github.com/mrdoob/three.js), r160: MIT; copyright 2010–2023 Three.js Authors.
- [meshoptimizer](https://github.com/zeux/meshoptimizer): MIT; copyright 2016–2026 Arseny Kapoulkine.
- [IBM Plex](https://github.com/IBM/plex): SIL Open Font License 1.1; copyright 2017 IBM Corp., reserved font name “Plex”.

Full bundled permission, copyright, warranty and font notices are in
[THIRD-PARTY-LICENSES.txt](THIRD-PARTY-LICENSES.txt).
The Open 3D Model GPL-3.0 viewer and Z-Anatomy application code are not used.

## Attribution when sharing anatomy

Retain the relevant upstream credits, license links and modification notices,
and identify your further changes. A reasonable credit for the combined
anatomy is “Anatomy adapted by Brian Pridgen from Open 3D Model and Z-Anatomy,
with BodyParts3D and other upstream contributors; component licenses include
CC BY-SA 4.0, CC BY-NC 4.0 (kidneys, Lissie Cowley), and CC BY-NC-SA 4.0
(inner ear, University of Dundee School of Medicine); see credits.” Include this document or a link to its hosted
copy. Videos, screenshots, downloads and other shared outputs containing the kidney
or inner-ear models must retain their attribution and applicable license
notices and comply with their NonCommercial terms.
This example does not replace the applicable license terms.

## Intended use and warranties

Educational and reference use; the viewer is not intended for diagnosis,
clinical decision-making or intraoperative guidance. Anatomy may be incomplete
or inaccurate. No upstream institution or contributor is represented as
endorsing this viewer. Material is provided as is, subject to the warranty
and liability provisions of its applicable license. This intended-use notice
adds no restriction to rights granted by the anatomy licenses.
