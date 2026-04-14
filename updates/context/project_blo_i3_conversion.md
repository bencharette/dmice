---
name: BLO i3 conversion debugging session
description: Debugging session fixing blo_npz_to_i3.py so steamshovel displays tracks correctly
type: project
---

## Status (as of 2026-04-03)

**RESOLVED** — BLO muon simulations now display correctly in steamshovel with proper track direction and vertex position.

Test files on Cobalt: `~/dmice_work/output/blo_dmice_targeted_det1det2_both_1000events.i3.zst` (1000 events), `~/dmice_work/output/blo_muons_200hits.i3.zst` (2 events)
Local copies: `/home/ben/dmice/output/`

**How to convert NPZ to I3 on Cobalt:**
```bash
scp ~/dmice/BLO/blo_npz_to_i3.py ~/dmice/BLO/icecube_with_dmice.geo cobalt:~/
ssh cobalt "/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
  python ~/blo_npz_to_i3.py --geo ~/icecube_with_dmice.geo \
  ~/dmice_work/output/INPUT_repacked.npz \
  ~/dmice_work/output/OUTPUT.i3.zst"
```

**Steamshovel filter for 250+ DOM hits:**
Window → Frame Filter → `len(frame["InIcePulses"]) >= 250`

## Bugs fixed in this session

### 1. `omgeo.omtype` not set (both scripts)
- **File:** `BLO/blo_npz_to_i3.py` and `prometheus_to_i3.py`
- **Bug:** `I3OMGeo.omtype` left at default `Unknown` — steamshovel DOM artists skip unknown-type DOMs entirely, so nothing displayed.
- **Fix:** Set `omgeo.omtype = dataclasses.I3OMGeo.IceCube` in geometry loop.

### 2. `primary.time` not set (both scripts)
- **File:** `BLO/blo_npz_to_i3.py` and `prometheus_to_i3.py`
- **Bug:** `primary.time` left as NaN — steamshovel can't anchor the particle in time.
- **Fix:** Set to earliest DOM hit time: `t0 = np.min(dom_t[i])`.

### 3. Zenith direction wrong in `prometheus_to_i3.py`
- **File:** `prometheus_to_i3.py`
- **Bug:** Prometheus parquet `initial_state_zenith` is in anti-momentum (incoming) convention. Was used raw — needed π-flip to convert to I3Direction convention.
- **Fix:** `primary.dir = I3Direction(np.pi - zenith, azimuth)`.

### 4. Vertex hardcoded to (0, 0, 0) in `blo_npz_to_i3.py`
- **Bug:** BLO detector center is ~(-25, 38, -39) in IceCube coords — track at origin missed detector.
- **Fix:** Vertex set to 1500 m before the hit-DOM centroid along the momentum direction.

### 5. Zenith convention in `blo_npz_to_i3.py` (most complex bug)
- **Bug:** `batch_dm_ice_sim.py` stores `zenith_rad = arccos(dz)` where `dz > 0` = upgoing momentum direction from +z. `I3Direction` uses anti-momentum (incoming) convention. The two conventions require opposite treatment.
- **Fix:**
  - `primary.dir = I3Direction(np.pi - zen_mom, azi)` — π-flip for I3Direction
  - Vertex backstep uses raw momentum direction (no flip): `dz_mom = cos(zenith_rad[i])`
- **Why separate:** The vertex geometry must use the physical momentum direction; I3Direction must use incoming direction.

## Session 2026-04-03 — RESOLVED

### Bug 6: `dx_mom` used before assignment
- **File:** `BLO/blo_npz_to_i3.py`
- **Bug:** `primary.dir = I3Direction(-dx_mom, -dy_mom, -dz_mom)` was on line 135, but `dx_mom` etc. were calculated on lines 148-150.
- **Fix:** Moved momentum direction calculation earlier, right after extracting zenith/azimuth.

### Bug 7: Track not visible in steamshovel
- **Bug:** `primary.length` was not set — steamshovel needs explicit length to render tracks.
- **Fix:** Added `primary.length = 10000.0 * I3Units.m` (10km track).

### Bug 8: Track direction opposite to DOM timing
- **Bug:** Using anti-momentum direction `I3Direction(-dx_mom, -dy_mom, -dz_mom)` caused track to point backwards.
- **Fix:** Use momentum direction: `I3Direction(dx_mom, dy_mom, dz_mom)`.
- **Key insight:** For steamshovel visualization, use actual momentum direction, NOT the I3Direction anti-momentum convention.

### Bug 9: Track vertex position
- **Bug:** With momentum direction, vertex must be BEHIND hit centroid so track extends forward through detector.
- **Fix:** `vx = cx - BACKSTEP * dx_mom` (subtract, not add) with `BACKSTEP = 1500.0` metres.

## Final working configuration (blo_npz_to_i3.py)
```python
# Momentum direction components
dx_mom = np.sin(zen_mom) * np.cos(azi_ic)
dy_mom = np.sin(zen_mom) * np.sin(azi_ic)
dz_mom = np.cos(zen_mom)

primary.length = 10000.0 * I3Units.m
primary.dir = dataclasses.I3Direction(dx_mom, dy_mom, dz_mom)  # momentum direction

# Vertex 1.5km behind hit centroid
vx = cx - BACKSTEP * dx_mom
vy = cy - BACKSTEP * dy_mom
vz = cz - BACKSTEP * dz_mom
```

## Key data facts (blo_muons_200hits_repacked.npz)
- 2 events
- ev0: E=125 TeV, zen=1.116 rad (64°), azi=5.395 rad — upgoing, dz=0.44
- ev1: E=111 TeV, zen=0.666 rad (38°), azi=0.805 rad — upgoing, dz=0.79
- DOM hit times ev0: 2109–8822 ns, ev1: 236–8584 ns
- Hit centroid ev0: x=-25, y=38, z_ic=-39 m (≈ IceCube center)
- Vertex after backstep ev0: z_ic ≈ -699 m (below detector)

## Convention reference — CONFIRMED from dmice_muons_filtered.i3
- `I3Direction.zenith` = arccos(-dir.z) — IceCube convention, measured from above
- `I3Direction(zenith_ic, azi)` stores dir.z = -cos(zenith_ic)
- For downgoing muon (zenith_ic=87°): dir.z = -cos(87°) = -0.054 ✓ (confirmed from real data)
- `batch_dm_ice_sim.py` zenith: angle of momentum from +z (0=up, 90=horizontal), only upgoing [0,90°]
- Conversion to I3Direction: `I3Direction(π - zen_blo, azi)` — π-flip IS correct
- `prometheus_to_i3.py` zenith: standard Prometheus anti-momentum — needs same π-flip

## Reference data (dmice_muons_filtered.i3 on Cobalt)
Nearly-horizontal downgoing CORSIKA muons (zenith ~85-87°).
- MC vertex far outside detector: e.g., (-1106, -1403, -430) m at t=5586 ns
- MMC entry point at detector boundary: e.g., (-477, -642, -483) m
- SplineMPE fit vertex inside detector: e.g., (-10, -77, -520) m
- Key: real MC vertex is ~800m outside detector, not at the boundary
