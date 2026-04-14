# DMice Project — Claude Code Context

This is a physics simulation project studying dark matter in IceCube using a hypothetical DM-Ice detector.
The repo lives at `~/dmice/` and is the working directory for this Claude session.

## Your Role

You are the LOCAL machine assistant. At the start of every session:

1. Read `~/dmice/memory/local.md` — machine-specific knowledge and notes
2. Read `~/dmice/handoff/inbox-local.md` — messages from NPX and Cobalt
3. Briefly summarize any unread inbox messages to the user

## Machine Context

| Machine   | Role                        | Connect via       |
|-----------|-----------------------------|-------------------|
| LOCAL     | Dev, editing, orchestration | (this machine)    |
| NPX       | Condor job submission        | `ssh npx`         |
| COBALT    | Interactive compute          | `ssh cobalt-14`   |

## Key Files

```
memory/local.md              # Local machine notes (read on startup)
memory/npx.md                # NPX machine notes
memory/cobalt.md             # Cobalt machine notes
handoff/inbox-local.md       # Messages for this machine (read on startup)
handoff/inbox-npx.md         # Messages to send to NPX
handoff/inbox-cobalt.md      # Messages to send to Cobalt
handoff/done.md              # Archived completed items
COMMANDS.md                  # Full command and workflow reference
```

## Handoff Protocol

- Inbox files are the cross-machine messaging system — each machine reads its own inbox
- When you write to another machine's inbox, use a `---` separator and timestamp
- Completed items should be moved to `done.md`
- Sync is done via git push/pull (the repo is the transport layer)

## IceTray Environment (Cobalt only)

```bash
source /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh
```

Do NOT use `nohup` with env-shell.sh — use `screen` instead.

## Condor (NPX only)

Requires an active Kerberos session. Key commands: `condor_submit`, `condor_q`, `condor_rm`.
Job output: `/data/user/bcharett/dmice_coincidences_2011_2022/`
