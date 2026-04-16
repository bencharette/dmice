#!/usr/bin/env python3
"""
profile_era_keys.py

For each year, sample step3 coinc files and catalog:
  - sub_event_streams present
  - pulse keys present per stream
  - reconstruction keys present per stream
  - P-frames per DAQ event

Run on Cobalt inside IceTray environment:
  /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
      python3 ~/dmice/profile_era_keys.py
"""

import glob, os
from collections import defaultdict
from icecube import icetray, dataio

BASEDIR    = "/data/user/bcharett/dmice_coincidences_2011_2022/step3_coincidences"
YEARS      = list(range(2012, 2022))
FILES_PER_YEAR = 40   # sample size per year

PULSE_CANDIDATES = [
    "OfflinePulses", "SRTInIcePulses", "ReextractedInIcePulses",
    "InIcePulses", "InIceDSTPulses", "SplitInIcePulses",
    "SplitInIceDSTPulses", "TWOfflinePulses_FR_WIMP",
    "RTTWOfflinePulses_FR_WIMP",
]
RECO_CANDIDATES = [
    "LineFit", "PoleMuonLinefitCutsNanoDST", "PoleMuonLinefit",
    "MPEFit", "SPEFit2", "SPEFitSingle",
    "FiniteRecoFit", "I3DST", "I3DST13",
]

for year in YEARS:
    files = sorted(glob.glob(os.path.join(BASEDIR, str(year), "*", "*", "*.i3.zst")))
    if not files:
        print(f"\n{year}: no files found")
        continue

    # Sample files with actual content
    sample = []
    for f in files:
        if len(sample) >= FILES_PER_YEAR:
            break
        if os.path.getsize(f) > 5000:
            sample.append(f)

    # Per-stream stats
    stream_pframes    = defaultdict(int)   # stream → count of P-frames
    stream_pulse_keys = defaultdict(lambda: defaultdict(int))  # stream → key → count
    stream_reco_keys  = defaultdict(lambda: defaultdict(int))
    daq_events        = 0
    pframes_per_daq   = []
    current_daq_p     = 0

    for f in sample:
        try:
            i3f = dataio.I3File(f)
            while i3f.more():
                frame = i3f.pop_frame()
                if frame.Stop == icetray.I3Frame.DAQ:
                    if current_daq_p > 0:
                        pframes_per_daq.append(current_daq_p)
                    current_daq_p = 0
                    daq_events += 1
                elif frame.Stop == icetray.I3Frame.Physics:
                    hdr = frame["I3EventHeader"]
                    stream = getattr(hdr, "sub_event_stream", "") or "NoSplit"
                    stream_pframes[stream] += 1
                    current_daq_p += 1
                    for k in PULSE_CANDIDATES:
                        if k in frame:
                            stream_pulse_keys[stream][k] += 1
                    for k in RECO_CANDIDATES:
                        if k in frame:
                            stream_reco_keys[stream][k] += 1
            if current_daq_p > 0:
                pframes_per_daq.append(current_daq_p)
            i3f.close()
        except Exception as e:
            pass

    total_p = sum(stream_pframes.values())
    avg_p   = sum(pframes_per_daq) / len(pframes_per_daq) if pframes_per_daq else 0

    print(f"\n{'='*60}")
    print(f"YEAR {year}  ({len(sample)} files, {daq_events} DAQ events, "
          f"{total_p} P-frames, {avg_p:.1f} P/DAQ)")
    print(f"{'='*60}")

    for stream in sorted(stream_pframes, key=lambda s: -stream_pframes[s]):
        n = stream_pframes[stream]
        print(f"\n  Stream '{stream}'  ({n} frames, {n/max(total_p,1)*100:.0f}%)")

        pulses = {k: v for k, v in stream_pulse_keys[stream].items() if v > 0}
        if pulses:
            print(f"    Pulse keys:  " +
                  "  ".join(f"{k}({v})" for k, v in sorted(pulses.items(), key=lambda x: -x[1])))
        else:
            print(f"    Pulse keys:  (none of the candidates)")

        recos = {k: v for k, v in stream_reco_keys[stream].items() if v > 0}
        if recos:
            print(f"    Reco keys:   " +
                  "  ".join(f"{k}({v})" for k, v in sorted(recos.items(), key=lambda x: -x[1])))
        else:
            print(f"    Reco keys:   (none of the candidates)")
