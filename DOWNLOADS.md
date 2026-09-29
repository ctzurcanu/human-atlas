# Downloads

The **Download** panel exports the anatomy view or the visible atlas page as PNG. Video includes the anatomy, labels, and visible controls. Frames are composed locally in the webpage, so recording needs no screen-sharing permission or extension and does not resize the window. The 1 px red recording outline is excluded from PNG and video frames.

The video size selector fixes the saved file's pixel dimensions. When its aspect ratio differs from the window, the complete view fits inside the file with margins. **Open window at this size** requests a separate window at the chosen dimensions; the browser and screen may limit the actual window size. File dimensions still use the selected preset.

Video targets 60 fps and captures scene frames immediately after rendering, including intermediate camera and slide transitions. Labels are painted synchronously with the anatomy instead of being decoded as separate image snapshots. The actual frame rate depends on the device's rendering and encoding capacity.

Video uses a resolution-scaled encoding bitrate: approximately 50 Mbps for Full HD at 60 fps, with a 12 Mbps minimum and a 100 Mbps maximum. Actual bitrate depends on the browser encoder and scene motion. The anatomy, labels, and interface are rendered at the selected output density instead of enlarging a lower-resolution frame. Full HD and 4K presets provide more detail at the cost of larger files and more rendering work; the current-window default remains available.

Slide titles, counters, and advanced controls refresh separately from the rest of the interface during transitions. Larger interface snapshots update between camera transitions, including dropdowns and scrolled panels. Still frames reuse the composed output rather than repeatedly copying or rendering the anatomy. This captures the atlas's own page; it does not capture browser chrome, browser dialogs, or an embedding page's surrounding content. The renderer and technical sources are credited in [ATTRIBUTION.md](ATTRIBUTION.md).

## Video

In .webm format. To transform:

```
ffmpeg \
  -i human-atlas-male-detail-2026-09-29T10-08-16.486Z.webm \
  -map 0:v:0 \
  -vsync 0 \
  -c:v hevc_videotoolbox \
  -b:v 10M \
  -tag:v hvc1 \
  -color_range pc \
  -colorspace bt709 \
  -color_primaries bt709 \
  -color_trc bt709 \
  -video_track_timescale 1000 \
  human-atlas-male-detail.mp4
```

```
ffmpeg -i human-atlas-male-detail-2026-09-29T10-08-16.486Z.webm \
  -vf "minterpolate=fps=60:mi_mode=mci:mc_mode=aobmc:me_mode=bidir" \
  -c:v hevc_videotoolbox \
  -b:v 12M \
  -tag:v hvc1 \
  -color_range pc \
  -colorspace bt709 \
  -color_primaries bt709 \
  -color_trc bt709 \
  human-atlas-male-detail-60fps.mp4
```
