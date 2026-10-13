import {createRoot} from 'react-dom/client';
import Home from '../app/page';
import {initializeAssetCache} from '../app/asset-cache';
import {assetUrl} from '../app/asset-url';
import {setAppViewerUrl} from '../app/viewer-location';
import {followSystemTheme} from '../app/theme';
import {connectUiHost,openUiLink} from './ui-bridge.mjs';
import '../app/globals.css';
import './theme.css';

const root = createRoot(document.getElementById('root')!);
const status = document.getElementById('status')!;
const stopTheme = followSystemTheme(document, window.matchMedia('(prefers-color-scheme: dark)'));
let generation = 0;
const ready = initializeAssetCache(assetUrl(''));
const error = (text: string) => {status.textContent = text; status.hidden = false;};
document.addEventListener('click', event => {
 const link = event.target instanceof Element ? event.target.closest<HTMLAnchorElement>('a[target="_blank"]') : null;
 if (!link) return;
 event.preventDefault();
 openUiLink(window, link.href);
});
connectUiHost(window, {
 async show(url: string) {
  const current = ++generation;
  // Local review meshes stay in the local atlas, outside remote host sandboxes.
  if (new URL(url).origin === 'http://localhost:3016') {
   root.render(null);
   const link = document.createElement('a');
   link.href = url; link.target = '_blank'; link.rel = 'noopener noreferrer';
   link.textContent = 'Open local anatomy view';
   status.replaceChildren(link); status.hidden = false;
   return;
  }
  setAppViewerUrl(url);
  await ready;
  if (current !== generation) return;
  status.hidden = true;
  root.render(<Home key={url}/>);
 },
 error,
 context(value: {theme?: string}) {
  if (value.theme !== 'dark' && value.theme !== 'light') return;
  stopTheme();
  document.documentElement.dataset.theme = value.theme;
  document.documentElement.classList.toggle('dark', value.theme === 'dark');
 },
 teardown() {generation++; stopTheme(); root.unmount();},
});
