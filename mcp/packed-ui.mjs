/** Native decompression works in browsers and Workers without eval or network code. */
export async function unpackRenderer(encoded) {
 const bytes = Uint8Array.from(atob(encoded), value => value.charCodeAt(0));
 const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
 return JSON.parse(await new Response(stream).text());
}
