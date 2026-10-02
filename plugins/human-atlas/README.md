# Atlas

Atlas lets Claude search named anatomical structures and open them in an interactive 3D human atlas. Ask for an organ, a system or a region, and Atlas shows it inline with orbit, zoom, isolation, depth layers, explode and cross-sections.

## What it does

- Searches anatomical names and identifiers to find the exact structure you mean.
- Shows one or several structures in a focused 3D view.
- Offers detailed male, complete male, standard male, female, embryo and cell models. Coverage varies, and some reference models are partial.

Try: "Show the stomach in isolation." or "Show the stomach and pancreas together."

## Tools

All three tools are read-only.

- `get_anatomy_options` lists the models, systems, regions, layers and camera views.
- `search_anatomy` finds exact structure names and IDs.
- `show_anatomy` opens the interactive viewer.

## What it connects to

The plugin declares one remote MCP server in `.mcp.json`: `https://human-atlas-connect.ctzurcanu.workers.dev/mcp`, hosted on Cloudflare Workers. It needs no account, sign-in or credentials. Your search terms and view options are sent to that server to build a view. The viewer loads from `https://ctzurcanu.github.io/human-atlas/`. The plugin runs no local code, hooks or scripts, and stores nothing.

## Limits

Atlas provides general anatomy education. It does not access patient records, interpret medical scans, diagnose conditions or recommend treatment.

## Links

- Website: https://ctzurcanu.github.io/human-atlas/
- Support: https://ctzurcanu.github.io/human-atlas/support.html
- Privacy: https://ctzurcanu.github.io/human-atlas/privacy.html
- Terms: https://ctzurcanu.github.io/human-atlas/terms.html
- Source: https://github.com/ctzurcanu/human-atlas

## License

GPL-3.0-only. Vitruvian artwork is CC BY-NC 4.0; see `assets/NOTICE.txt`.
