import type {Part} from './anatomy';

/** Reviewed native source extent; this does not certify a complete nerve. */
export function nerveSourceCoverageLimitation(part:Part):string|undefined{
 if(part.system==='nervous'&&/^LOCAL:female:[lr]_palmar_branch_ulnar_nerve:\d+$/.test(part.id))return 'The source labeled “Palmar branch ulnar nerve” is an unpartitioned distal digital network alongside the ring and little fingers. It is provisionally grouped under the superficial ulnar branch, not identified as the palmar cutaneous branch. Separate common/proper digital meshes, full fascicles and physical nerve junctions remain unverified.';
}
