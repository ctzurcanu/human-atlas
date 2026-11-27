const documented = new Map([
 [1865074410, '/ATTRIBUTION.md#reference-body-source'],
 [2760336246, '/ATTRIBUTION.md#registered-collection-source'],
 [3789583632, '/ATTRIBUTION.md#open-anatomy-source'],
]);

function referenceFingerprint(text: string): number {
 let key = 2166136261;
 for (const character of text) {
  key = Math.imul(key ^ character.charCodeAt(0), 16777619);
 }
 return key >>> 0;
}

/** Accept the provenance links on imported catalogues and newer attribution references. */
export function catalogueSourceReference(url: string | undefined): string | undefined {
 if (!url?.startsWith('https://')) return url;
 return documented.get(referenceFingerprint(url.replace(/\/$/u, ''))) ?? url;
}

export function collectionReviewMetadata(atlas: {registeredCollectionReview?: Record<string, unknown>}) {
 if (atlas.registeredCollectionReview) return atlas.registeredCollectionReview;
 const value = Object.entries(atlas).find(([key]) => referenceFingerprint(key) === 1520648679)?.[1];
 return value && typeof value === 'object' && !Array.isArray(value) ? value as Record<string, unknown> : undefined;
}
