# Hepatic fitting: adjoining TA98 inventory

Additional authority descendants used alongside [TA98-HEPATIC-PATHS.md](TA98-HEPATIC-PATHS.md). A generic exterior does not establish wall layers or microscopic features. Openings are not solids.

| TA98 | Name | Parent | Existing male source |
|---|---|---|---|
| A05.8.02.001 | gallbladder | A05.0.00.000 | ZA:Gallbladder |
| A05.8.02.002 | fundus of gallbladder | A05.8.02.001 | TA98PART:ZA:Gallbladder:A05.8.02.002 |
| A05.8.02.003 | infundibulum of gallbladder | A05.8.02.001 | TA98PART:ZA:Gallbladder:A05.8.02.003 |
| A05.8.02.004 | body of gallbladder | A05.8.02.001 | TA98PART:ZA:Gallbladder:A05.8.02.004 |
| A05.8.02.005 | neck of gallbladder | A05.8.02.001 | TA98PART:ZA:Gallbladder:A05.8.02.005 |
| A05.8.02.006 | serosa | A05.8.02.001 | missing |
| A05.8.02.007 | subserosa | A05.8.02.001 | missing |
| A05.8.02.008 | muscular layer | A05.8.02.001 | missing |
| A05.8.02.009 | mucosa | A05.8.02.001 | missing |
| A05.8.02.010 | mucosal folds | A05.8.02.009 | missing |
| A05.8.02.011 | cystic duct | A05.8.02.001 | FJ3080 |
| A05.8.02.012 | spiral fold | A05.8.02.011 | missing |
| A05.8.02.013 | bile duct | A05.8.02.001 | ZA:Bile duct; TA98PART:FJ3079:A05.8.02.013 |
| A05.8.02.014 | sphincter of bile duct | A05.8.02.013 | missing |
| A05.8.02.015 | superior sphincter | A05.8.02.014 | missing |
| A05.8.02.016 | inferior sphincter | A05.8.02.014 | missing |
| A05.8.02.017 | hepatopancreatic ampulla | A05.8.02.013 | missing |
| A05.8.02.018 | sphincter of ampulla | A05.8.02.017 | missing |
| A05.8.02.019 | glands of bile duct | A05.8.02.017 | missing |
| A10.1.02.101 | lesser omentum | A10.1.02.002 | ZA:Lesser omentum |
| A10.1.02.102 | hepatophrenic ligament | A10.1.02.101 | missing |
| A10.1.02.103 | hepato-oesophageal ligament | A10.1.02.101 | missing |
| A10.1.02.104 | hepatogastric ligament | A10.1.02.101 | missing |
| A10.1.02.105 | hepatoduodenal ligament | A10.1.02.101 | missing |
| A10.1.02.106 | hepatocolic ligament | A10.1.02.101 | missing |
| A12.2.12.015 | common hepatic artery | A12.2.12.012 | ZA:Common hepatic artery |
| A12.2.12.016 | gastroduodenal artery | A12.2.12.015 | ZA:Gastroduodenal artery |
| A12.2.12.017 | supraduodenal artery | A12.2.12.016 | missing |
| A12.2.12.018 | posterior superior pancreaticoduodenal artery | A12.2.12.016 | FJ3557 |
| A12.2.12.019 | pancreatic branches | A12.2.12.018 | missing |
| A12.2.12.020 | duodenal branches | A12.2.12.018 | missing |
| A12.2.12.021 | retroduodenal arteries | A12.2.12.016 | missing |
| A12.2.12.022 | right gastro-omental artery | A12.2.12.016 | FJ3596 |
| A12.2.12.023 | gastric branches | A12.2.12.022 | missing |
| A12.2.12.024 | omental branches | A12.2.12.022 | missing |
| A12.2.12.025 | anterior superior pancreaticoduodenal artery | A12.2.12.016 | FJ3409 |
| A12.2.12.026 | pancreatic branches | A12.2.12.025 | missing |
| A12.2.12.027 | duodenal branches | A12.2.12.025 | missing |
| A12.2.12.028 | right gastric artery | A12.2.12.015 | FJ3594 |
| A12.2.12.029 | hepatic artery proper | A12.2.12.015 | ZA:Proper hepatic artery |
| A12.2.12.030 | right branch | A12.2.12.029 | FJ3075; FJ3076; FJ3114; FJ3115; FJ3116; FJ3117 |
| A12.2.12.031 | cystic artery | A12.2.12.030 | missing |
| A12.2.12.032 | artery of caudate lobe | A12.2.12.030 | FJ1916; FJ2386 |
| A12.2.12.033 | anterior segmental artery | A12.2.12.030 | missing |
| A12.2.12.034 | posterior segmental artery | A12.2.12.030 | missing |
| A12.2.12.035 | left branch | A12.2.12.029 | FJ3091; FJ3092; FJ3093; FJ3095; FJ3106; FJ3107 |
| A12.2.12.036 | artery of caudate lobe | A12.2.12.035 | FJ3077 |
| A12.2.12.037 | medial segmental artery | A12.2.12.035 | missing |
| A12.2.12.038 | lateral segmental artery | A12.2.12.035 | missing |
| A12.2.12.039 | intermediate branch | A12.2.12.029 | missing |
| A13.3.03.031 | cystic node | A13.3.03.030 | ZA:Cystic node |

## Trial assembly and interfaces

Preserve the whole liver exterior, ligament material and eight overlapping source territories; they are not new exclusive tissue partitions. Apply the same spatial fitting function to all twenty inventoried native hepatic dependants, preserving their suppression status. Exact source-coordinate constraints must remain exact. Ports, microscopic features and source open rims remain unaccepted; no artificial closure disks are created.

Probe an anterior Z displacement of at most 15 mm, with smoothstep support X=[0,0.020] m and Y=[1.230,1.270] m. X/Y stay fixed. There is no Z dependence, so the continuous Jacobian determinant is one. Preserve UVs, indices and materials; transport normals by the inverse transpose. Vessel/omental endpoint courses remain unverified and must be reviewed before acceptance.

Gates: source replay; shared coordinates; unchanged topology; triangle inversion/degeneration; whole liver self-crossings; gastric and fixed-neighbor surface crossings; actual render. Trial files remain local and uninstalled.

Authority SHA256: b4c41fa4fe13f68e2735bc59c94a744d58957e1e0dcf4e2172415a400f61ff48
