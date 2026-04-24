# pipeline/ — Pipeline Orchestration Scripts

Shell scripts that run multi-step workflows across machines (Cobalt, NPX, WARD).

## Scripts

| File | Machine | Description |
|------|---------|-------------|
| `run_phase1_pipeline.sh` | NPX | Phase 1 validation: parquet → i3 → linefit comparison for 20 sim runs |
| `run_2020_2021_pipeline.sh` | NPX | Full 2020-2021 coincidence pipeline (step1 Condor → step2 → step3 Condor → merge) |
| `pipeline_2020_2021_followup.sh` | Cobalt | Step2 + step3 + merge after step1 Condor jobs finish |
| `phase34_run.sh` | Cobalt | Merge 2020-2021 and re-run all recos on full dataset |
| `overnight_run.sh` | Cobalt | Full overnight: recos on fixed file → wait for step3 → merge 2020-2021 → recos on full dataset |
| `run_prometheus_to_i3.sh` | NPX/local | Convert Prometheus parquet → npz → i3 (two-step wrapper) |
| `test_coinc_per_year.sh` | Cobalt | Quick per-year coincidence test on one representative subrun |

## Typical workflow order

```
NPX: run_2020_2021_pipeline.sh        # submits and monitors Condor
  or
Cobalt: overnight_run.sh              # runs recos + waits for Condor + merges
```
