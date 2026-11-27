# Median nerve source interfaces

Authority SHA256: `b4c41fa4fe13f68e2735bc59c94a744d58957e1e0dcf4e2172415a400f61ff48`. All descendants of A14.2.03.031 are inventoried before any reconstruction. Existing identity mappings do not certify geometry.

| TA98 | Parent | Latin / English | Male source meshes | Female source meshes |
| --- | --- | --- | --- | --- |
| A14.2.03.031 | A14.2.03.020 | n. medianus / median nerve | ZA:Median nerve.l (exact), ZA:Median nerve.r (exact) | LOCAL:female:r_median_nerve:798 (exact), LOCAL:female:l_median_nerve:807 (exact) |
| A14.2.03.032 | A14.2.03.031 | radix medialis nervi mediani / medial root of median nerve | O3M:Medial root of median nerve.r (exact), O3M:Medial root of median nerve.l (exact) | No active mapped source |
| A14.2.03.033 | A14.2.03.031 | radix lateralis nervi mediani / lateral root of median nerve | O3M:Lateral root of median nerve.r (exact), O3M:Lateral root of median nerve.l (exact) | No active mapped source |
| A14.2.03.034 | A14.2.03.031 | n. interosseus antebrachii anterior / anterior interosseous nerve | ZA:Anterior interosseous nerve of forearm.l (exact), ZA:Anterior interosseous nerve of forearm.r (exact) | No active mapped source |
| A14.2.03.035 | A14.2.03.031 | rr. musculares / muscular branches | ZA:Muscular branches of median nerve.l (exact), ZA:Muscular branches of median nerve.r (exact) | No active mapped source |
| A14.2.03.036 | A14.2.03.031 | r. palmaris / palmar branch | ZA:Palmar branch of median nerve.l (exact), ZA:Palmar branch of median nerve.r (exact) | No active mapped source |
| A14.2.03.037 | A14.2.03.031 | r. communicans ulnaris / communicating branch with ulnar nerve | No active mapped source | No active mapped source |
| A14.2.03.038 | A14.2.03.031 | nn. digitales palmares communes / common palmar digital nerves | ZA:Common palmar digital branches of median nerve.l (exact), ZA:Common palmar digital branches of median nerve.r (exact) | LOCAL:female:l_common_palmar_digital_branches_median_nerve:782 (exact), LOCAL:female:r_common_palmar_digital_branches_median_nerve:789 (exact) |
| A14.2.03.039 | A14.2.03.031 | nn. digitales palmares proprii / proper palmar digital nerves | ZA:Proper palmar digital branches of median nerve.l (exact), ZA:Proper palmar digital branches of median nerve.r (exact) | LOCAL:female:l_proper_palmar_digital_branches_median_nerve:783 (exact), LOCAL:female:r_proper_palmar_digital_branches_median_nerve:790 (exact) |

The male ZA median trunks and O3M roots are from separate publishers. The O3M left roots are mirrored from the right source. The importer does not establish accepted cross-source terminal registration. No source geometry was changed.

- O3M:Lateral root of median nerve.l: 45 triangle intersections with the ZA trunk; zero exact shared used coordinates; nearest used vertex distance 0.1411 mm. The distance is not a surface-distance or anatomical-clearance measurement.
- O3M:Medial root of median nerve.l: 0 triangle intersections with the ZA trunk; zero exact shared used coordinates; nearest used vertex distance 7.8880 mm. The distance is not a surface-distance or anatomical-clearance measurement.
- O3M:Lateral root of median nerve.r: 46 triangle intersections with the ZA trunk; zero exact shared used coordinates; nearest used vertex distance 0.1411 mm. The distance is not a surface-distance or anatomical-clearance measurement.
- O3M:Medial root of median nerve.r: 0 triangle intersections with the ZA trunk; zero exact shared used coordinates; nearest used vertex distance 7.8880 mm. The distance is not a surface-distance or anatomical-clearance measurement.

Do not restore a named cross-source branching link as if the mixed assembly were accepted. Native same-frame hand branches remain available in both models. The broad relation validator now checks the excluded roots reciprocally while retaining their existing selectable meshes.

Before root reconstruction: define both source trunks, both roots, proximal cords, distal branches and adjoining interfaces; review actual source surfaces in Blender and cross-sections; inventory shared tissue and source ownership; then repair the complete assembly before splitting it. Nerve exterior meshes do not establish epineurium, perineurium, endoneurium or fascicles. An open end is an interface opening, not an additional solid.

Blender inspection was attempted but the Mac was locked. Render and anatomical acceptance remain pending; the source-interface report is quantitative diagnosis only.
