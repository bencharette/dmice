---
name: Use dmice sim_linefit_comparison.py not dmice_results version
description: Always use sim_linefit_comparison.py from /home/bench/dmice/, not the copy in dmice_results
type: feedback
---

Use `/home/bench/dmice/sim_linefit_comparison.py` for simulation linefit comparisons, not the version in `/home/bench/dmice_results/`.

**Why:** The dmice/ version is the better/more current one. There are two copies and the user explicitly confirmed this preference.

**How to apply:** When writing scripts or pipeline steps that call sim_linefit_comparison.py, always reference the one in `/home/bench/dmice/`.
