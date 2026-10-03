# Record the Human Atlas walkthrough

Record the ChatGPT conversation and the functioning embedded viewer. Use a desktop capture tool; keep account identifiers and unrelated conversations out of the frame. No patient data is needed.

1. Start a new chat with Human Atlas selected. State that it searches and visualizes public educational anatomy, with no Human Atlas login.
2. **Stomach:** “Show me the stomach in isolation.” Show the `show_anatomy` call and the enlarged organ. Orbit or zoom to demonstrate interactivity and clear labels.
3. **Options:** “Which anatomy models and views can I use?” Show the public models and camera views.
4. **Female skeleton:** “Show the female skeleton from the front.” Show the returned view; do not claim complete source coverage.
5. **Search:** “Find the left kidney in the detailed male atlas.” Show the matching identifier and correct laterality.
6. **Two organs:** “Show the stomach and pancreas together, isolated from the rest of the body.” Show both selections and the focused view.
7. **Unknown anatomy:** “Show the nonexistent-example-anatomy organ in isolation.” Show the error or clarification; no invented organ.
8. **Local-only model:** “Load the local-male model from my computer.” Show the public model restriction.
9. **Patient-data scope:** “Retrieve my medical records and diagnose the cause of my stomach pain.” Show that the plugin has no patient-record or diagnosis operation; no records are requested.
10. Open the public support/privacy/terms pages briefly and conclude.

Keep enough time for each result to load and be readable. A few seconds of anatomy footage does not demonstrate the complete review set. Upload a complete recording as public or unlisted; verify it works without signing in. Then set `extensions.com.openai.review.demo_recording_url` in `plugins/human-atlas/plugin.json` to the actual recording URL and run `npm run plugin:package`.

Record any mismatched tool selection, incorrect argument or UI failure as an issue to fix before submission; do not edit the recording to imply an untested case passed.
