# DMice Local Machine Memory

Machine: **LOCAL**

Role: Development, editing, orchestration, syncing handoffs.

## Machine Info

- Used for development and script editing
- Orchestrates job submission to NPX
- Syncs handoff files to NPX and Cobalt
- Primary command center for the pipeline

## Commands

- `ssh npx` — connect to NPX submit node
- `ssh cobalt-14` — connect to Cobalt compute node
- `nvim-dmice-local` — open Neovim with local profile
- `nvim-dmice-npx` — open Neovim with NPX profile
- `nvim-dmice-cobalt` — open Neovim with Cobalt profile

## Handoff Files

- `handoff/inbox-local.md` — messages from other machines for this machine
- `handoff/inbox-npx.md` — messages to NPX
- `handoff/inbox-cobalt.md` — messages to Cobalt

## Notes

Add project-specific knowledge here as it develops.
