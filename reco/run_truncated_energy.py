#!/usr/bin/env python3
"""
run_truncated_energy.py

Runs TruncatedEnergy reconstruction on the merged DM-Ice coincidence I3 file.
Uses MPEFit (already in file) as seed + OfflinePulses as hits.
Writes a CSV with energy, n_doms, n_hits per event.

Usage (on Cobalt with IceTray env):
    /cvmfs/.../env-shell.sh python3 run_truncated_energy.py
"""

import os
import csv

SPICEMIE_DIR    = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/SPICEMie"
SPICEMIE_DRIVER = os.path.join(SPICEMIE_DIR, "driverfiles")

GCD_FILE    = "/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2013.56429_V1.i3.gz"
DEFAULT_IN      = "/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022.i3.zst"
DEFAULT_OUT_I3  = "/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_with_energy.i3.zst"
DEFAULT_OUT_CSV = os.path.expanduser("~/dmice_work/output/comparison/real_hits_energy_v2.csv")

from icecube import icetray, dataio, dataclasses, recclasses
from icecube import truncated_energy
from icecube.icetray import I3Tray

os.makedirs(os.path.dirname(DEFAULT_OUT_CSV), exist_ok=True)

rows = []

ENERGY_OUT = "Reprocessed_TruncatedEnergy_SPICEMie"

@icetray.traysegment
def Truncated(tray, Name, Pulses="", Seed="", PhotonicsService=""):
    tray.AddModule(
        "I3TruncatedEnergy",
        RecoPulsesName     = Pulses,
        RecoParticleName   = Seed,
        ResultParticleName = ENERGY_OUT,
        I3PhotonicsServiceName = PhotonicsService,
        UseRDE             = True,
        # Only run if seed exists AND we haven't already got energy for this event
        If = lambda f: Seed in f and ENERGY_OUT + "_ORIG_Muon" not in f,
    )

def extract(frame):
    hdr  = frame["I3EventHeader"]
    year = str(hdr.start_time.utc_year)
    run  = hdr.run_id
    evt  = hdr.event_id

    energy_GeV = float("nan")
    # Try reprocessed key first, then fall back to original 2012 L2 key
    for key in [ENERGY_OUT + "_ORIG_Muon", "MPEFitTruncatedEnergy_SPICEMie_ORIG_Muon"]:
        try:
            energy_GeV = frame[key].energy
            if energy_GeV > 0:
                break
        except Exception:
            pass

    n_doms = 0
    n_hits = 0
    try:
        pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, UNIFIED_PULSES)
        n_doms = len(pulses)
        n_hits = sum(len(v) for v in pulses.values())
    except Exception:
        pass

    rows.append(dict(year=year, run=run, event=evt,
                     energy_GeV=energy_GeV, n_doms=n_doms, n_hits=n_hits))

tray = I3Tray()

tray.AddService("I3PhotonicsServiceFactory", "PhotonicsServiceMu",
    PhotonicsTopLevelDirectory = SPICEMIE_DIR,
    DriverFileDirectory        = SPICEMIE_DRIVER,
    PhotonicsLevel2DriverFile  = "mu_photorec.list",
    PhotonicsTableSelection    = 2,
    ServiceName                = "PhotonicsServiceMu",
)

tray.AddModule("I3Reader", FilenameList=[GCD_FILE, DEFAULT_IN])

# Normalise pulse key: create "UnifiedPulses" from year-appropriate source
PULSE_PRIORITY = ["OfflinePulses", "SRTInIcePulses", "ReextractedInIcePulses", "InIcePulses"]
UNIFIED_PULSES = "UnifiedPulses"

def unify_pulses(frame):
    for key in PULSE_PRIORITY:
        if key in frame:
            frame[UNIFIED_PULSES] = frame[key]
            return
tray.AddModule(unify_pulses, Streams=[icetray.I3Frame.Physics])

tray.AddSegment(Truncated,
    Pulses           = UNIFIED_PULSES,
    Seed             = "MPEFit",
    PhotonicsService = "PhotonicsServiceMu",
)

tray.AddModule(extract, Streams=[icetray.I3Frame.Physics])

tray.AddModule("I3Writer", Filename=DEFAULT_OUT_I3,
               Streams=[icetray.I3Frame.Physics])

tray.Execute()
tray.Finish()

with open(DEFAULT_OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["year","run","event","energy_GeV","n_doms","n_hits"])
    w.writeheader()
    w.writerows(rows)

n_with_e = sum(1 for r in rows if str(r["energy_GeV"]) not in ("nan",""))
print(f"Done: {len(rows)} events, {n_with_e} with energy")
print(f"CSV saved: {DEFAULT_OUT_CSV}")
print(f"I3  saved: {DEFAULT_OUT_I3}")
