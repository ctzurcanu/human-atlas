import {createPortal} from 'react-dom';
import {useEffect,useRef,useState,type RefObject} from 'react';
import {Download,ExternalLink,Square,Video,X} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {VIDEO_FPS,captureScale,drawPageFrame,renderAdvancedOverlay,renderPageOverlay,renderSceneLabels,videoBitrate,type SceneFrame,type SceneFrameCapture} from './page-capture';

export interface SceneCaptureOptions {labels:boolean;background:boolean}
export type SceneCapture=(options?:SceneCaptureOptions)=>Promise<Blob>;
type Dimensions={width:number;height:number};
interface Props {
 studioRef:RefObject<HTMLElement|null>;
 captureRef:RefObject<SceneCapture|null>;
 frameCaptureRef:RefObject<SceneFrameCapture|null>;
 currentViewUrl:()=>string;
 model:string;
 ready:boolean;
}
interface ActiveRecording {
 recorder:MediaRecorder;
 output:MediaStream;
 frame:number;
 cleanup:()=>void;
}

const VIDEO_SIZES=[
 {value:'1920x1080',label:'1920 × 1080 · Full HD',width:1920,height:1080},
 {value:'3840x2160',label:'3840 × 2160 · 4K',width:3840,height:2160},
 {value:'1280x720',label:'1280 × 720 · HD',width:1280,height:720},
 {value:'1440x900',label:'1440 × 900 · Laptop',width:1440,height:900},
 {value:'1080x1080',label:'1080 × 1080 · Square',width:1080,height:1080},
 {value:'1080x1920',label:'1080 × 1920 · Portrait',width:1080,height:1920},
 {value:'720x1280',label:'720 × 1280 · Portrait HD',width:720,height:1280},
] as const;
const initialVideoSize=()=>VIDEO_SIZES.some(item=>item.value===new URLSearchParams(location.search).get('recsize'))?new URLSearchParams(location.search).get('recsize')!:'current';
const windowSize=():Dimensions=>({width:window.innerWidth,height:window.innerHeight});
const fileName=(model:string,extension:string)=>`human-atlas-${model}-${new Date().toISOString().replaceAll(':','-')}.${extension}`;
export function downloadBlob(blob:Blob,name:string){
 const url=URL.createObjectURL(blob),link=document.createElement('a');
 link.href=url;link.download=name;document.body.appendChild(link);link.click();link.remove();
 window.setTimeout(()=>URL.revokeObjectURL(url),60_000);
}
const waitFrame=()=>new Promise<void>(resolve=>requestAnimationFrame(()=>requestAnimationFrame(()=>resolve())));
const pngBlob=(canvas:HTMLCanvasElement)=>new Promise<Blob>((resolve,reject)=>canvas.toBlob(blob=>blob?resolve(blob):reject(new Error('The PNG could not be created.')),'image/png'));
function recordingMime(){
 if(typeof MediaRecorder==='undefined')throw new Error('Video recording is unavailable in this browser.');
 for(const mime of ['video/webm;codecs=vp8','video/webm;codecs=vp9','video/webm','video/mp4'])if(MediaRecorder.isTypeSupported(mime))return mime;
 return '';
}

export default function ExportTools({studioRef,captureRef,frameCaptureRef,currentViewUrl,model,ready}:Props){
 const [open,setOpen]=useState(false),[source,setSource]=useState<'view'|'screen'>('view');
 const [labels,setLabels]=useState(true),[background,setBackground]=useState(true),[seconds,setSeconds]=useState(0);
 const [viewport,setViewport]=useState<Dimensions>(windowSize),[videoSize,setVideoSize]=useState(initialVideoSize);
 const [recording,setRecording]=useState(false),[busy,setBusy]=useState<'png'|'video'|null>(null),[status,setStatus]=useState(''),[videoFile,setVideoFile]=useState<{url:string;name:string}|null>(null),[videoData,setVideoData]=useState('');
 const pending=useRef<AbortController|null>(null),statusTimer=useRef<number|undefined>(undefined),recorderRef=useRef<ActiveRecording|null>(null);
 useEffect(()=>{const update=()=>setViewport(windowSize());window.addEventListener('resize',update);window.visualViewport?.addEventListener('resize',update);return()=>{window.removeEventListener('resize',update);window.visualViewport?.removeEventListener('resize',update);};},[]);
 const notice=(message:string)=>{setStatus(message);window.clearTimeout(statusTimer.current);if(message)statusTimer.current=window.setTimeout(()=>setStatus(''),5000);};
 const cleanRecording=(active:ActiveRecording)=>{cancelAnimationFrame(active.frame);active.cleanup();active.output.getTracks().forEach(track=>track.stop());};
 const stopRecording=()=>{
  const active=recorderRef.current;if(!active)return;
  recorderRef.current=null;
  if(active.recorder.state!=='inactive')active.recorder.stop();
  cleanRecording(active);setRecording(false);notice('Finishing video…');
 };
 useEffect(()=>()=>{pending.current?.abort();const active=recorderRef.current;if(active){if(active.recorder.state!=='inactive')active.recorder.stop();cleanRecording(active);}window.clearTimeout(statusTimer.current);},[]);
 const selectedVideoSize=():Dimensions=>VIDEO_SIZES.find(item=>item.value===videoSize)??viewport;
 const openSizedWindow=()=>{
  const size=VIDEO_SIZES.find(item=>item.value===videoSize);if(!size)return;
  const url=new URL(currentViewUrl());url.searchParams.set('recsize',size.value);
  const chromeWidth=window.top===window?Math.max(0,outerWidth-innerWidth):0,chromeHeight=window.top===window?Math.max(0,outerHeight-innerHeight):0;
  const popup=window.open(url.href,'human-atlas-recording',`popup=yes,resizable=yes,width=${size.width+chromeWidth},height=${size.height+chromeHeight}`);
  if(!popup){notice('The browser blocked the recording window.');return;}
  const adjust=()=>{try{popup.resizeTo(size.width+Math.max(0,popup.outerWidth-popup.innerWidth),size.height+Math.max(0,popup.outerHeight-popup.innerHeight));}catch{}};
  popup.addEventListener('load',adjust,{once:true});window.setTimeout(()=>{adjust();try{if(Math.abs(popup.innerWidth-size.width)>16||Math.abs(popup.innerHeight-size.height)>16)notice('The browser limited the window size; the video preset still applies.');}catch{}},700);popup.focus();setOpen(false);
 };
 const pagePng=async()=>{
  const studio=studioRef.current,frameCapture=frameCaptureRef.current;
  if(!studio||!frameCapture)throw new Error('The page is still loading.');
  try{const frame=frameCapture(),width=studio.clientWidth,height=studio.clientHeight;
   const [pageLabels,overlay]=await Promise.all([renderSceneLabels(frame.labels,width,height),renderPageOverlay(studio)]);
   const canvas=document.createElement('canvas');canvas.width=width;canvas.height=height;
   drawPageFrame(canvas,frameCapture().canvas,pageLabels,overlay,width,height);
   return pngBlob(canvas);
  }finally{frameCapture.release();}
 };
 const takePng=async()=>{
  if(!ready||!captureRef.current){notice('The anatomy is still loading.');return;}
  const capture=captureRef.current,task=new AbortController();pending.current=task;setBusy('png');
  try{
   setOpen(false);await waitFrame();
   const delay=Math.min(300,Math.max(0,Number(seconds)||0));
   if(delay){notice(`PNG in ${delay} second${delay===1?'':'s'}…`);await new Promise<void>((resolve,reject)=>{const timer=window.setTimeout(resolve,delay*1000);task.signal.addEventListener('abort',()=>{window.clearTimeout(timer);reject(new DOMException('Cancelled','AbortError'));},{once:true});});}
   if(task.signal.aborted)return;
   const blob=source==='view'?await capture({labels,background}):await pagePng();
   if(task.signal.aborted)return;
   downloadBlob(blob,fileName(model,'png'));notice('PNG downloaded.');
  }catch(error){if(!task.signal.aborted)notice(error instanceof Error?error.message:'The PNG could not be created.');}
  finally{if(pending.current===task)pending.current=null;setBusy(null);}
 };
 const startVideo=async()=>{
  const studio=studioRef.current,frameCapture=frameCaptureRef.current;
  if(!ready||!studio||!frameCapture){notice('The anatomy is still loading.');return;}
  if(recorderRef.current||busy)return;
  const size=selectedVideoSize();setBusy('video');setOpen(false);
  let output:MediaStream|undefined,frame=0,cleanup=()=>{};
  try{
   const mime=recordingMime(),canvas=document.createElement('canvas');canvas.width=size.width;canvas.height=size.height;
   await waitFrame();
   const sourceWidth=studio.clientWidth,sourceHeight=studio.clientHeight,scale=captureScale(size,sourceWidth,sourceHeight),initial=frameCapture(size);
   let [overlay,advancedOverlay]=await Promise.all([renderPageOverlay(studio,scale,true),renderAdvancedOverlay(studio,scale)]);
   drawPageFrame(canvas,frameCapture(size).canvas,initial.labels,overlay,sourceWidth,sourceHeight,advancedOverlay);
   const manualFrames=typeof CanvasCaptureMediaStreamTrack!=='undefined'&&typeof CanvasCaptureMediaStreamTrack.prototype.requestFrame==='function';
   output=canvas.captureStream(manualFrames?0:VIDEO_FPS);
   const videoTrack=output.getVideoTracks()[0] as CanvasCaptureMediaStreamTrack;
   const recorder=new MediaRecorder(output,{...(mime?{mimeType:mime}:{}),videoBitsPerSecond:videoBitrate(size)}),pieces:Blob[]=[];
   recorder.ondataavailable=event=>{if(event.data.size)pieces.push(event.data);};
   recorder.onerror=()=>{stopRecording();notice('The video encoder stopped recording.');};
   recorder.onstop=()=>{if(pieces.length){const type=recorder.mimeType||mime,extension=type.includes('mp4')?'mp4':'webm',blob=new Blob(pieces,{type}),name=fileName(model,extension),url=URL.createObjectURL(blob);setVideoFile(current=>{if(current)URL.revokeObjectURL(current.url);return {url,name};});void blob.arrayBuffer().then(buffer=>{let binary='';const bytes=new Uint8Array(buffer);for(let i=0;i<bytes.length;i+=0x8000)binary+=String.fromCharCode(...bytes.subarray(i,i+0x8000));setVideoData(btoa(binary));});if(location.hostname==='localhost')void fetch('/recording-upload',{method:'POST',headers:{'Content-Type':type},body:blob}).catch(()=>{});downloadBlob(blob,name);notice('Video downloaded.');}else notice('No video frames were recorded.');};
   let overlayDirty=false,overlayPending=false,advancedDirty=false,advancedPending=false,nextFrame=0,lastOverlay=0,currentFrame=initial,boundCapture:SceneFrameCapture|null=null,unsubscribe=()=>{};
   const present=(now:number,refresh=false)=>{
    const active=recorderRef.current;if(!active||active.recorder!==recorder||now+.5<nextFrame)return;
    nextFrame=Math.max(nextFrame+1000/VIDEO_FPS,now+500/VIDEO_FPS);
    // Re-submit a still frame without reading the WebGL buffer or redrawing
    // the anatomy. The buffer is valid only immediately after scene rendering.
    try{if(refresh)canvas.getContext('2d')?.drawImage(canvas,0,0,1,1,0,0,1,1);if(manualFrames)videoTrack.requestFrame();}
    catch(error){stopRecording();notice(error instanceof Error?error.message:'The page could not be recorded.');}
   };
   const paint=(current:SceneFrame,now:number)=>{try{drawPageFrame(canvas,current.canvas,current.labels,overlay,studio.clientWidth,studio.clientHeight,advancedOverlay);present(now);}catch(error){stopRecording();notice(error instanceof Error?error.message:'The page could not be recorded.');}};
   const bindCapture=(capture:SceneFrameCapture)=>{unsubscribe();boundCapture?.release();boundCapture=capture;currentFrame=capture(size);paint(currentFrame,performance.now());unsubscribe=capture.subscribe((current,now)=>{currentFrame=current;paint(current,now);},size);};
   const markOverlay=()=>{overlayDirty=true;advancedDirty=true;};
   const observer=new MutationObserver(mutations=>{for(const mutation of mutations){const element=mutation.target instanceof Element?mutation.target:mutation.target.parentElement;if(element?.closest('.anatomy-labels,[data-scene-canvas],.recording-outline'))continue;if(element?.closest('.advanced-panel'))advancedDirty=true;else markOverlay();}});
   observer.observe(studio.ownerDocument.body,{subtree:true,childList:true,characterData:true,attributes:true});
   document.addEventListener('scroll',markOverlay,true);document.addEventListener('input',markOverlay,true);document.addEventListener('change',markOverlay,true);
   window.addEventListener('resize',markOverlay);
   cleanup=()=>{unsubscribe();observer.disconnect();document.removeEventListener('scroll',markOverlay,true);document.removeEventListener('input',markOverlay,true);document.removeEventListener('change',markOverlay,true);window.removeEventListener('resize',markOverlay);boundCapture?.release();if(boundCapture!==frameCapture)frameCapture.release();};
   recorderRef.current={recorder,output,frame:0,cleanup};
   recorder.start(1000);setRecording(true);notice('');
   bindCapture(frameCapture);
   const draw=(now:number)=>{
    const active=recorderRef.current;if(!active||active.recorder!==recorder)return;
    try{if(frameCaptureRef.current){
     if(frameCaptureRef.current!==boundCapture)bindCapture(frameCaptureRef.current);
     present(now,true);
     const scale=captureScale(size,studio.clientWidth,studio.clientHeight);
     if(advancedDirty&&!advancedPending){advancedDirty=false;advancedPending=true;void renderAdvancedOverlay(studio,scale).then(image=>{advancedOverlay=image;boundCapture?.invalidate();}).catch(error=>{stopRecording();notice(error instanceof Error?error.message:'The presentation controls could not be rendered.');}).finally(()=>{advancedPending=false;});}
     if(overlayDirty&&!overlayPending&&!currentFrame.transitioning&&now-lastOverlay>=100){overlayDirty=false;overlayPending=true;lastOverlay=now;void renderPageOverlay(studio,scale,true).then(image=>{overlay=image;boundCapture?.invalidate();}).catch(error=>{stopRecording();notice(error instanceof Error?error.message:'The page controls could not be rendered.');}).finally(()=>{overlayPending=false;});}
    }}catch(error){stopRecording();notice(error instanceof Error?error.message:'The page could not be recorded.');return;}
    active.frame=requestAnimationFrame(draw);
   };
   frame=requestAnimationFrame(draw);recorderRef.current.frame=frame;
  }catch(error){cancelAnimationFrame(frame);cleanup();frameCapture.release();output?.getTracks().forEach(track=>track.stop());recorderRef.current=null;notice(error instanceof Error?error.message:'The video could not start.');}
  finally{setBusy(null);}
 };
 const pressButton=()=>{if(recording){stopRecording();return;}if(busy==='png'){pending.current?.abort();notice('PNG cancelled.');return;}if(!busy)setOpen(value=>!value);};
 return <>
  <Button variant="ghost" className={`icon-button ${open?'active':''} ${recording?'is-recording':''}`} onClick={pressButton} disabled={busy==='video'} aria-label={recording?'Stop recording':busy==='png'?'Cancel PNG':busy==='video'?'Preparing recording':'Download'} title={recording?'Stop recording':busy==='png'?'Cancel PNG':busy==='video'?'Preparing recording':'Download'}>{recording?<Square size={18} fill="currentColor"/>:<Download size={18}/>}</Button>
  {open&&createPortal(<aside className="download-panel glass" role="dialog" aria-label="Download">
   <div className="panel-heading download-heading"><strong>Download</strong><Button variant="ghost" className="icon-button" onClick={()=>setOpen(false)} aria-label="Close Download"><X size={18}/></Button></div>
   <div className="download-content">
    <fieldset><legend>PNG includes</legend><label><input type="radio" name="png-source" checked={source==='view'} onChange={()=>setSource('view')}/>Anatomy view</label><label><input type="radio" name="png-source" checked={source==='screen'} onChange={()=>setSource('screen')}/>Screen and controls</label></fieldset>
    {source==='view'&&<fieldset><legend>PNG view options</legend><label><input type="checkbox" checked={labels} onChange={event=>setLabels(event.target.checked)}/>Labels</label><label><input type="checkbox" checked={background} onChange={event=>setBackground(event.target.checked)}/>Background</label></fieldset>}
    <label className="download-delay">Take PNG after <span><input type="number" min="0" max="300" step="1" value={seconds} onChange={event=>setSeconds(Math.min(300,Math.max(0,Number(event.target.value)||0)))} aria-label="Seconds before PNG"/> seconds <small>(0 = now)</small></span></label>
    <Button variant="outline" onClick={()=>void takePng()} disabled={!ready||!!busy}><Download size={16}/> Download PNG</Button>
    <div className="download-divider"/>
    <div className="download-window-dimensions">Window <output>{viewport.width} × {viewport.height} px</output></div>
    <label className="download-video-size">Video size<select value={videoSize} onChange={event=>setVideoSize(event.target.value)} aria-label="Video size"><option value="current">Current window ({viewport.width} × {viewport.height})</option>{VIDEO_SIZES.map(size=><option key={size.value} value={size.value}>{size.label}</option>)}</select></label>
    {videoSize!=='current'&&<Button variant="outline" onClick={openSizedWindow}><ExternalLink size={16}/> Open window at this size</Button>}
    <Button variant="outline" onClick={()=>void startVideo()} disabled={!ready||!!busy}><Video size={16}/> Record video</Button>
   </div>
  </aside>,studioRef.current??document.body)}
  {recording&&createPortal(<div className="recording-outline" aria-hidden="true"/>,document.body)}
  {status&&createPortal(<p className="export-status glass" role="status">{status}</p>,studioRef.current??document.body)}
  {videoFile&&createPortal(<a className="export-status glass" role="link" aria-label="Download recorded video" href={videoFile.url} download={videoFile.name}>Download recorded video</a>,studioRef.current??document.body)}
  {videoData&&createPortal(<pre data-testid="recorded-video-base64" style={{position:'fixed',left:'-10000px',top:0,width:1,height:1,overflow:'hidden'}}>{videoData}</pre>,document.body)}
</>;
}
