"""
load_dmice_artists.py — Startup script for Steamshovel

Registers the ICLineFitArtist and PivotLineFitArtist and adds both to the scene.

Usage:
    steamshovel -s ~/dmice/load_dmice_artists.py your_file.i3.zst

Or paste into the Steamshovel Python console:
    exec(open('/home/bcharett/dmice/load_dmice_artists.py').read())
"""

import sys
sys.path.insert(0, '/home/bcharett/dmice')

from steamshovel_artists import ICLineFitArtist, PivotLineFitArtist

scenario = window.gl.scenario
scenario.registerArtist(ICLineFitArtist)
scenario.registerArtist(PivotLineFitArtist)
scenario.add('ICLineFitArtist')
scenario.add('PivotLineFitArtist')

print("DMice artists loaded.")
print("  ICLineFitArtist  — blue arrow (IC-only LineFit)")
print("  PivotLineFitArtist — red arrow (DM-Ice Pivot LineFit)")
print("Both arrows are centred at the IC charge-weighted centroid.")
print("The red sphere marks the DM-Ice detector pivot point.")
