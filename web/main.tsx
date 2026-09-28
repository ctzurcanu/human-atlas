import {createRoot} from 'react-dom/client';
import Home from '../app/page';
import {assetUrl} from '../app/asset-url';
import {initializeAssetCache} from '../app/asset-cache';
import '../app/globals.css';
void initializeAssetCache(assetUrl('')).finally(()=>createRoot(document.getElementById('root')!).render(<Home/>));
