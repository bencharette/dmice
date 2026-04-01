# Agent Update Instructions

**Date:** 2026-04-01
**Purpose:** Replace `run_2020_2021_pipeline.sh` in the repo root with the updated version in this folder.

---

## What Changed

`run_2020_2021_pipeline.sh` was refactored so it runs **on NPX** instead of on the local machine.

**Before:** The script ran locally and used `ssh npx "..."` for every condor command and `ssh cobalt-14 "..."` for all file I/O on `/data/user/`.

**After:**
- All `ssh npx "..."` wrappers removed — condor commands (`condor_submit`, `condor_q`, `step3_submit.py`) now run directly since the script executes on NPX.
- `ssh cobalt-14 "..."` removed for file operations on `/data/user/` (shared filesystem, accessible from NPX).
- `ssh cobalt-14` **kept** for step 2 (`step2_run.sh`) and step 4 (IceTray merge) — these genuinely need cobalt's environment.
- Header comment updated to reflect the new run target (NPX, not local).

---

## Task: Replace the Script

1. Copy `updates/run_2020_2021_pipeline.sh` → repo root `run_2020_2021_pipeline.sh`
2. Verify the diff looks correct (no `ssh npx` remaining, `ssh cobalt-14` only on steps 2 and 4)
3. Commit and push

```bash
cp updates/run_2020_2021_pipeline.sh run_2020_2021_pipeline.sh
git diff run_2020_2021_pipeline.sh   # review
git add run_2020_2021_pipeline.sh
git commit -m "refactor: run pipeline on NPX, ssh only to cobalt for step2+merge"
git push
```

---

## Context

See `updates/context/` for full project memory files. Start with `project_dmice.md` for complete project context, then `MEMORY.md` as an index to the rest.
