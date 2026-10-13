import {createRoot} from 'react-dom/client';
import Home from '../app/page';
import {assetUrl} from '../app/asset-url';
import {initializeAssetCache} from '../app/asset-cache';
import '../app/globals.css';
import {followSystemTheme} from '../app/theme';
followSystemTheme(document,window.matchMedia('(prefers-color-scheme: dark)'));
void initializeAssetCache(assetUrl('')).finally(()=>createRoot(document.getElementById('root')!).render(<Home/>));
