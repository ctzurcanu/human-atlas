import {connectUiHost} from './ui-bridge.mjs';
const viewer = document.getElementById('viewer');
const status = document.getElementById('status');
const open = document.getElementById('open');
connectUiHost(window, {
 show(url) {
  viewer.src = url; viewer.hidden = false; status.hidden = true;
  open.href = url; open.hidden = false;
 },
 error(text) {status.textContent = text; status.hidden = false;},
 context(value) {
  if (value.theme === 'dark' || value.theme === 'light') document.documentElement.style.colorScheme = value.theme;
 },
 teardown() {viewer.src = 'about:blank';},
});
