# tools/ — Utilities and Dev Tools

Helper scripts for visualization, workspace setup, and memory sync.

## Scripts

| File | Description |
|------|-------------|
| `steamshovel_artists.py` | Custom Steamshovel artists: `ICLineFitArtist` and `PivotLineFitArtist` for DM-Ice track display |
| `load_dmice_artists.py` | Convenience loader — registers steamshovel artists in the current session |
| `start-work.sh` | Launch the dmice tmux workspace (nvim + claude side-by-side) |
| `sync-memory.sh` | Sync local Claude memory files to Cobalt/NPX via scp |

## Notes

- `steamshovel_artists.py` is auto-registered by `startup.py` on Cobalt
- Direction convention: `ICLineFitArtist` uses the IceTray direction (points *toward* source); `PivotLineFitArtist` shows the pivot-shifted version
- `start-work.sh` usage: `bash ~/dmice/tools/start-work.sh [local|npx|cobalt]`
