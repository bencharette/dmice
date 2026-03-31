# DMice NPX Memory

Machine: **NPX** (npx-submitter.icecube.wisc.edu)

Role: Condor submission and job monitoring.

## Machine Info

- Condor job submission node
- Pipeline orchestration
- Submits jobs to the HTCondor cluster
- Monitors job progress

## Key Paths

- Working directory: `~/dmice_work/`
- Job logs and output: `/data/user/bcharett/dmice_coincidences_2011_2022/`

## Handoff Files

- `handoff/inbox-npx.md` — messages from other machines (LOCAL, COBALT) for NPX
- `handoff/inbox-local.md` — write back to LOCAL
- `handoff/inbox-cobalt.md` — write to COBALT

## Condor Commands

Requires interactive Kerberos session:
- `condor_submit` — submit a job
- `condor_q` — check job queue
- `condor_rm` — remove a job

## Notes

Add NPX-specific knowledge here as needed.
