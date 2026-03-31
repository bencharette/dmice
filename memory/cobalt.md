# DMice Cobalt Memory

Machine: **COBALT** (cobalt-14)

Role: Interactive compute environment.

## Machine Info

- Cobalt compute node (cobalt-14 is primary)
- Interactive script execution
- Long-running sessions via screen/tmux
- Debugging and pipeline testing

## Key Paths

- Working directory: `~/dmice_work/`
- IceTray env-shell: `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh`
- Prometheus: `~/prometheus` (cloned from GitHub: Harvard-Neutrino/prometheus)

## Handoff Files

- `handoff/inbox-cobalt.md` — messages from other machines (LOCAL, NPX) for Cobalt
- `handoff/inbox-local.md` — write back to LOCAL
- `handoff/inbox-npx.md` — write to NPX

## Interactive Sessions

Use `screen` for long-running tasks:
```bash
screen -S <session-name>
# inside screen:
source /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh
python script.py
```

Do NOT use `nohup` with env-shell.sh — it forks and output is lost.

## Notes

Add Cobalt-specific knowledge here as needed.
