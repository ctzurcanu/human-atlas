/** Keep canvas renderers, captures and UI appearance on the same system theme. */
export function followSystemTheme(document: Document, media: MediaQueryList) {
 const apply=()=>{
  const dark=media.matches;
  document.documentElement.dataset.theme=dark?'dark':'light';
  document.documentElement.classList.toggle('dark',dark);
  document.querySelector('meta[name="theme-color"]')?.setAttribute('content',dark?'#141b23':'#f2f3f3');
 };
 apply();
 media.addEventListener('change',apply);
 return ()=>media.removeEventListener('change',apply);
}
