# Cloudflare Connect setup

The same Worker also serves the public, unauthenticated Streamable HTTP MCP endpoint at `/mcp`. MCP uses stateless JSON responses and does not allocate rooms or write to Durable Object storage. The three anatomy tools use compact catalogues generated automatically before deployment; the UI and large model geometry remain on GitHub Pages. `npm run test:connect-cloudflare` covers MCP discovery, search, focused views, UI resources, origin rejection, and the existing host/guest flows.

The atlas stays on GitHub Pages. This Worker handles live sessions at `/atlas-connect`, with one SQLite-backed Durable Object per invitation. It uses the Free plan's supported storage backend and WebSocket hibernation. Camera motion is held in socket attachments rather than written to SQLite on every frame. The initial implementation accepts up to 200 guests per room; actual capacity also depends on Cloudflare's daily quotas and payload sizes.

## 1. Connect your account

From the project directory:

```sh
npx wrangler login
npx wrangler whoami
```

Complete the authorization in Cloudflare's official browser page. `whoami` must show your account. Authentication is stored by Wrangler outside the repository.

## 2. Verify and deploy the relay

For a new account, first open [Workers & Pages](https://dash.cloudflare.com/?to=/:account/workers-and-pages) in the Cloudflare dashboard. This initializes the account's `workers.dev` subdomain. If Cloudflare asks for a subdomain, choose one; if it asks for a plan, select Free. A deployment error with code `10063` means this initialization has not happened yet.

```sh
npm run test:connect-cloudflare
npm run connect:cloudflare:check
npm run connect:cloudflare:deploy
```

The configuration is in `cloudflare/wrangler.jsonc`; its Worker name is `human-atlas-connect`. Wrangler creates the SQLite-backed Durable Object namespace automatically. If this is your first Worker, it may ask you to choose a `workers.dev` account subdomain. No custom domain or GitHub repository connection is needed for this CLI deployment.

Keep your account on Workers Free. This setup does not need a paid plan or other paid Cloudflare products. Free limits can interrupt operations until their daily reset; they do not become paid overages on the Free plan.

Wrangler prints the deployed HTTPS address, such as `https://human-atlas-connect.YOUR-SUBDOMAIN.workers.dev`. Open that address followed by `/health`; it should return `status: "ok"`. Use the actual address Wrangler prints, not this example.

Verify the deployed relay with `npm run test:connect-live -- https://YOUR-ACTUAL-RELAY.workers.dev`. This creates and stops a temporary host/guest test session over normally validated WSS connections.

## 3. Configure GitHub Pages

The current relay is [human-atlas-connect.ctzurcanu.workers.dev](https://human-atlas-connect.ctzurcanu.workers.dev). It is already the default in `app/connect-config.ts`. Publish the updated viewer code to GitHub Pages using the repository's normal deployment; no GitHub variable is required for this address. The current public viewer cannot use the Cloudflare relay until these client changes are deployed.

To override that default with another relay, in the `ctzurcanu/human-atlas` GitHub repository:

1. Open **Settings → Secrets and variables → Actions → Variables**.
2. Create a repository variable named `VITE_CONNECT_RELAY_URL`.
3. Set its value to the deployed Worker HTTPS address, without `/atlas-connect`.
4. Publish the updated viewer code using the **Deploy GitHub Pages** workflow. If these code changes are not on `main` yet, publish them first; rerunning the old code cannot use the new relay.

The Pages build reads that variable. It is a public server address, not a secret. A one-off alternative is adding `relay=<URL-encoded Worker HTTPS address>` to the viewer URL, or entering the Worker address in Connect's Host field.

## 4. Check a live session

1. Open the updated GitHub Pages viewer and choose **Advanced tools → Connect → Host**.
2. Leave the Host address field to prepare the lobby; Copy becomes available after it connects.
3. Copy the invitation and open it in another browser/tab/device. Choose Guest and press Play to join the lobby.
4. The host's count should show the guest. Press Host Play; move the camera or select an anatomical structure and verify the guest follows.
5. **H** reopens the host's tools while broadcasting. The actual Stop beside the Host field ends the session for guests.

Invitations retain the GitHub Pages viewer URL and carry both the random room code and the relay address. Cloudflare receives connection IPs and the shared viewer state. The relay accepts this project's GitHub Pages origin; change `ALLOWED_VIEWER_ORIGINS` if hosting the viewer elsewhere. Room state is removed after the disconnected host's 60-second grace period. Deploying an update disconnects existing transports; clients attempt to reconnect.

## Troubleshooting

- **Not authenticated:** complete `npx wrangler login` and check `npx wrangler whoami` in the same project directory.
- **Cannot reach the live relay:** check `/health`, the HTTPS address, and that the deployed viewer contains these Connect changes.
- **403 during a WebSocket connection:** check that the viewer origin is included in `ALLOWED_VIEWER_ORIGINS`.
- **Guest waiting during an active presentation:** Play captures guests online at that moment. Open **H → Connect** and press the single Play icon beside the guest list to include all waiting guests in the current broadcast. Its number shows how many guests will be included. They immediately receive the current anatomy, slides and camera, then follow subsequent changes.
- **Free quota exceeded:** review Workers and Durable Objects usage in Cloudflare; operations resume after quota reset. Hibernation saves idle compute, but active sessions still consume resources.

Official documentation links are in `ATTRIBUTION.md`.
