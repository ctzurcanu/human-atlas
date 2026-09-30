# Human Atlas submission handoff

Publisher: Christian Tzurcanu. Support: GitHub Issues, confirmed by the publisher.

## Prepared files

- `plugins/human-atlas/`: portable public package source, with root `plugin.json`, root `mcp.json`, the publisher's icon and a live-viewer screenshot.
- `outputs/human-atlas-1.0.0.zip`: generated upload artifact; regenerate with `npm run plugin:package`.
- `public/privacy.html`, `public/terms.html`, `public/support.html`: publication pages for the three listing URLs. These need deployment to GitHub Pages before submission.
- `test-cases.json`: exactly five positive and three negative review cases, also embedded in the package manifest.
- `tool-checks.json`: successful direct checks against the deployed MCP endpoint. These are not a claim that all eight conversational prompts passed in ChatGPT.
- `walkthrough.md`: recording instructions covering all eight cases.

The user confirmed the isolated-stomach flow works in ChatGPT. The other conversational cases and mobile UI still need evaluation in ChatGPT. The listing screenshot was captured from the deployed viewer, not from a mockup.

## Video

The current walkthrough is https://www.youtube.com/watch?v=8A56yWDgX0w and is set in `review.demo_recording_url`. The local Resolve-compatible master is `outputs/atlas-demo-careful.mp4`; use that master if the portal requests a direct upload. Verify the YouTube video's visibility and coverage before submission.

## Portal steps

1. Review the prepared listing and publication pages. Publish the pages through the existing GitHub Pages workflow; verify the support, privacy and terms URLs return the intended pages, not a 404.
2. Open https://platform.openai.com/plugins under the organization and project that will own Human Atlas. Select the publisher's verified individual identity. A ChatGPT Pro subscription does not by itself establish developer identity verification.
3. Upload the ZIP. Confirm the category against the dashboard's available categories (the package currently proposes Education).
4. In MCP setup, connect `https://human-atlas-connect.ctzurcanu.workers.dev/mcp` with No Authentication. No reviewer login credentials are necessary.
5. The portal provides a domain-verification challenge. Supply its exact token for the MCP hostname. The Worker has a `/.well-known/openai-apps-challenge` route driven by `OPENAI_APPS_CHALLENGE`; when unset, it returns 404. Do not use a made-up token.
6. Configure the actual token with `npx wrangler secret put OPENAI_APPS_CHALLENGE --config cloudflare/wrangler.jsonc`, then deploy the tested Worker with `npm run connect:cloudflare:deploy`. The challenge URL must return only the exact token as plain text.
7. Review imported cases, release notes, video and automated findings. Run all eight prompts in ChatGPT and retain results. Test desktop and mobile rendering.
8. Complete review details and any required attestations in the portal. Submit for review. After approval, choose Publish plugin.

No portal submission, policy attestation, identity verification or public-directory publication has been performed by this preparation task.

## Iframe justification

The optional MCP UI embeds `https://ctzurcanu.github.io/human-atlas/`. Christian Tzurcanu controls the viewer repository and its GitHub Pages deployment. The same publisher controls the Cloudflare MCP endpoint. The embedded page supplies the existing interactive Three.js anatomy renderer, geometry, selection details, orbit/zoom, layers and sections; embedding it lets the user inspect the requested structures within the conversation. It is essential interactive content rather than a static promotional page. Its origin is explicitly declared in the MCP resource's frame CSP. GitHub Pages and Cloudflare are different hosting origins even though the publisher controls both. This explanation is for the portal's iframe review.

## Rebuild and checks

```sh
npm run test:plugin-review
npm run test:connect-cloudflare
npm run plugin:package
```

The ZIP contains only public package assets and configurations. It does not contain Cloudflare credentials, local anatomy data, `.app.json` references, lifecycle hooks, geometry bundles or the repository checkout.

## Sources used

- https://developers.openai.com/plugins/deploy/submission
- https://developers.openai.com/plugins/build/plugins
- https://developers.openai.com/plugins/plugin-guidelines
