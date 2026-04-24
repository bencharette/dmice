#!/bin/bash
# Wrapper script for running the Prometheus muon simulation via HTCondor.
# Condor sets CWD to initialdir (~/dmice), so ./output/ resolves correctly.
set -e

ICETRAY_ENV=/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh

exec "${ICETRAY_ENV}" python simulate_muons.py
