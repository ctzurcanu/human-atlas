/** Keep anatomical text when an imported item includes a labelled source note. */
export function itemDescription(text = '') {
 if (!text || !/\bsources?\s*(?::|$)/imu.test(text)) return text || '';
 const lines = [];
 let inSources = false;
 for (const line of text.split(/\r?\n/u)) {
  if (/^\s*(?:#{1,6}\s+)?sources?\s*:?\s*$/iu.test(line)) {
   inSources = true;
   continue;
  }
  if (/^\s*#{1,6}\s+/u.test(line)) inSources = false;
  if (inSources) continue;
  const label = /(?:^\s*|(?<=[.!?;])\s+)sources?\s*:/iu.exec(line);
  if (!label) lines.push(line);
  else if (label.index) lines.push(line.slice(0, label.index));
 }
 return lines.join('\n').trim();
}
