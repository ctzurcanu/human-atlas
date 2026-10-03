/** TA98/source grouping is not always tissue containment. These reviewed
 * boundary pathways remain selectable neighbors, outside component assemblies.
 * Evidence: human hypothalamus atlas (PMC3654779), optic pathway sections
 * (PMID19837645), human stria medullaris tractography (PMC5952041).
 */
export const containmentExceptions=[
 {parent:'hypothalamus',child:'optic chiasm',note:'The optic chiasm is a visual-pathway structure at the inferior/anterior hypothalamic boundary; TA98 grouping does not establish hypothalamic tissue containment.'},
 {parent:'hypothalamus',child:'optic tract',note:'The optic tract is a visual-pathway structure along the hypothalamic boundary; it is not a subdivision of hypothalamic gray matter.'},
 {parent:'thalamus',child:'stria medullaris of thalamus',note:'The stria medullaris is an epithalamic white-matter pathway along the dorsal/medial thalamic surface; its name does not establish dorsal-thalamic tissue containment.'},
 {parent:'epiglottis',child:'hyoepiglottic ligament',note:'The hyoepiglottic ligament attaches the epiglottis to the hyoid; source/TA98 grouping does not establish a tissue subdivision contained within the epiglottis. Attachment interfaces remain unverified.'},
 {parent:'epiglottis',child:'thyroepiglottic ligament',note:'The thyroepiglottic ligament attaches the epiglottis to thyroid cartilage; source/TA98 grouping does not establish a tissue subdivision contained within the epiglottis. Attachment interfaces remain unverified.'},
] as const;
const base=(name:string|undefined)=>(name??'').toLowerCase().replace(/\s*\((?:left|right)\)\s*$/,'').replace(/^stria medullaris thalami$/,'stria medullaris of thalamus').replace(/^(hyo|thyro)-epiglottic ligament$/,'$1epiglottic ligament').trim();
export function isAorticBranchGrouping(child:string|undefined,parent:string|undefined):boolean{
 const c=base(child).replace(/^(?:left|right)\s+/,''),p=base(parent).replace(/^(?:left|right)\s+/,'').replace(/^inferior vena cava \((?:abdominal|thoracic) part\)$/,'inferior vena cava');
 if(p==='dorsal metatarsal arteries'&&['dorsal digital arteries','dorsal digital arteries of foot'].includes(c))return true;
 if(['azygos vein','hemi-azygos vein','hemiazygos vein'].includes(p)&&c==='ascending lumbar vein')return true;
 if(['inferior vena cava','ascending lumbar vein'].includes(p)&&['lumbar veins','ascending lumbar vein'].includes(c)&&c!==p)return true;
 // Source chapters place the entire downstream vascular tree under the
 // common carotid. Those distant branches are not contained in its tissue.
 // Preserve named tissue subdivisions, such as the carotid sinus or wall.
 if(p==='common carotid artery'&&c!==p&&!/^(?:wall|cervical part|thoracic part) of common carotid artery$/.test(c)&&(/\b(?:artery|arteries|arterial|branch|branches)\b/.test(c)||c.startsWith('?')))return true;
 return ['aortic arch','arch of aorta'].includes(p)&&['brachiocephalic trunk','common carotid artery','subclavian artery','proximal segment of subclavian artery'].includes(c)||p==='brachiocephalic trunk'&&['common carotid artery','subclavian artery'].includes(c)||p==='common carotid artery'&&['internal carotid artery','external carotid artery'].includes(c);
}
export function nonComponentRelationship(child:string|undefined,parent:string|undefined):string|undefined{
 if(base(child)==='proximal segment of subclavian artery'&&base(parent)==='subclavian artery')return 'These are adjoining source partitions. The named subclavian source remainder is not a verified whole assembly containing the proximal segment.';
 if(isAorticBranchGrouping(child,parent))return 'A vascular branch or tributary is an adjoining vessel, not a tissue component of the named vessel. Source chapter grouping does not establish physical containment.';
 const c=base(child),p=base(parent);return containmentExceptions.find(rule=>rule.child===c&&rule.parent===p)?.note;
}
