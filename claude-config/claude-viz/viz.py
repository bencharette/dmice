#!/usr/bin/env python3
"""
Claude Workflow Visualizer — live TUI for Claude Code sessions.

Usage:
  python3 ~/claude-viz/viz.py           # live view, latest session
  python3 ~/claude-viz/viz.py --setup   # print hook config instructions
  python3 ~/claude-viz/viz.py --all     # show all sessions
  python3 ~/claude-viz/viz.py --once    # render once and exit
  python3 ~/claude-viz/viz.py --clear   # wipe events log and exit
"""
import json
import os
import sys
import time
import argparse
from pathlib import Path
from collections import defaultdict, OrderedDict

try:
    from rich.console import Console
    from rich.live import Live
    from rich.table import Table
    from rich.panel import Panel
    from rich.columns import Columns
    from rich.text import Text
    from rich.rule import Rule
    from rich import box
except ImportError:
    print("pip install rich")
    sys.exit(1)

EVENTS_FILE = Path.home() / ".claude-viz" / "events.jsonl"
REFRESH_HZ = 2
MAX_TAIL_LINES = 2000  # only read last N lines for performance

# ── tool display config ────────────────────────────────────────────────────────

TOOL_ICON = {
    "Bash": "⚡", "Read": "󰈙", "Edit": "✎", "Write": "󰏫",
    "Glob": "󰥨", "Grep": "󰍉", "Agent": "󰙳", "WebFetch": "󰖟",
    "WebSearch": "󰖟", "TodoWrite": "󰃯", "TodoRead": "󰃯",
    "Task": "󰄹", "NotebookEdit": "󱞁", "Skill": "󰠱",
}
TOOL_COLOR = {
    "Bash": "yellow", "Read": "cyan", "Edit": "green", "Write": "green",
    "Glob": "blue", "Grep": "blue", "Agent": "magenta", "WebFetch": "bright_blue",
    "WebSearch": "bright_blue", "TodoWrite": "white", "TodoRead": "white",
}

def tool_icon(name):
    return TOOL_ICON.get(name, "○")

def tool_color(name):
    return TOOL_COLOR.get(name, "white")

# ── helpers ────────────────────────────────────────────────────────────────────

def summarize(tool, inp):
    """Extract the most human-readable snippet from tool input."""
    if not inp:
        return ""
    if tool == "Bash":
        c = inp.get("command", "")
        # collapse whitespace / newlines
        c = " ".join(c.split())
        return (c[:80] + "…") if len(c) > 80 else c
    if tool in ("Read", "Write", "Edit"):
        p = inp.get("file_path", "")
        return os.path.basename(p) or p
    if tool == "Glob":
        return inp.get("pattern", "")
    if tool == "Grep":
        pat = inp.get("pattern", "")
        path = inp.get("path", "")
        return f"{pat}  in {os.path.basename(path)}" if path else pat
    if tool == "Agent":
        d = inp.get("description", inp.get("prompt", ""))
        return (d[:80] + "…") if len(d) > 80 else d
    if tool in ("WebFetch", "WebSearch"):
        return inp.get("url", inp.get("query", ""))[:80]
    vals = [str(v) for v in inp.values() if v]
    if vals:
        v = vals[0]
        return (v[:80] + "…") if len(v) > 80 else v
    return ""

def fmt_dur(seconds):
    if seconds < 0.001:
        return "<1ms"
    if seconds < 1:
        return f"{seconds*1000:.0f}ms"
    if seconds < 60:
        return f"{seconds:.1f}s"
    m = int(seconds // 60)
    s = seconds % 60
    return f"{m}m{s:.0f}s"

def is_error(response_str):
    if not response_str:
        return False
    low = response_str.lower()
    return any(k in low for k in ("error", "exception", "failed", "traceback", "not found"))

# ── data loading ───────────────────────────────────────────────────────────────

def load_events():
    if not EVENTS_FILE.exists():
        return []
    with open(EVENTS_FILE) as f:
        lines = f.readlines()
    # tail for performance
    lines = lines[-MAX_TAIL_LINES:]
    out = []
    for line in lines:
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return out

def pair_calls(events):
    """Pair PreToolUse with PostToolUse events."""
    calls = []
    for e in events:
        if e["type"] == "PreToolUse":
            calls.append({"pre": e, "post": None})
        elif e["type"] == "PostToolUse":
            # match to most recent unpaired pre with same tool_name
            for c in reversed(calls):
                if (c["pre"] and
                        c["pre"]["tool_name"] == e["tool_name"] and
                        c["post"] is None):
                    c["post"] = e
                    break
    return calls

# ── rendering ──────────────────────────────────────────────────────────────────

def make_session_table(calls, max_rows=50):
    t = Table(
        box=box.SIMPLE_HEAD,
        show_header=True,
        header_style="bold dim",
        expand=True,
        padding=(0, 1),
        show_edge=False,
    )
    t.add_column("#", width=4, justify="right", style="dim")
    t.add_column("", width=2, no_wrap=True)       # icon
    t.add_column("Tool", width=12, no_wrap=True)
    t.add_column("Details", ratio=1, no_wrap=True)
    t.add_column("Duration", width=8, justify="right", no_wrap=True)
    t.add_column("", width=2, no_wrap=True)        # status

    shown = calls[-max_rows:]
    offset = max(0, len(calls) - max_rows)

    for idx, c in enumerate(shown, start=offset + 1):
        pre, post = c["pre"], c["post"]
        if not pre:
            continue

        tool = pre["tool_name"]
        icon = tool_icon(tool)
        color = tool_color(tool)
        detail = summarize(tool, pre.get("tool_input") or {})

        if post:
            dur_s = post["ts"] - pre["ts"]
            dur = fmt_dur(dur_s)
            err = is_error(str(post.get("tool_response", "")))
            if err:
                status = "[red]✗[/red]"
                tool_str = f"[dim {color}]{tool}[/dim {color}]"
                detail_str = f"[dim]{detail}[/dim]"
                dur_str = f"[dim]{dur}[/dim]"
            else:
                status = "[green]✓[/green]"
                tool_str = f"[{color}]{tool}[/{color}]"
                detail_str = detail
                dur_str = f"[dim]{dur}[/dim]"
        else:
            elapsed = time.time() - pre["ts"]
            dur_str = f"[bold yellow]{fmt_dur(elapsed)}…[/bold yellow]"
            status = "[bold yellow]⟳[/bold yellow]"
            tool_str = f"[bold {color}]{tool}[/bold {color}]"
            detail_str = f"[bold]{detail}[/bold]"

        t.add_row(
            str(idx),
            icon,
            tool_str,
            detail_str,
            dur_str,
            status,
        )

    return t


def make_stats_bar(calls, session_id, start_ts):
    total = len([c for c in calls if c["pre"]])
    done = sum(1 for c in calls if c["post"] is not None)
    errors = sum(
        1 for c in calls
        if c["post"] and is_error(str(c["post"].get("tool_response", "")))
    )
    running = total - done
    elapsed = fmt_dur(time.time() - start_ts)
    sid = session_id[:12] if session_id else "?"

    parts = [
        f"[dim]session[/dim] [cyan]{sid}[/cyan]",
        f"[dim]tools[/dim] [white]{done}/{total}[/white]",
        f"[dim]elapsed[/dim] [white]{elapsed}[/white]",
    ]
    if running:
        parts.append(f"[yellow]{running} running[/yellow]")
    if errors:
        parts.append(f"[red]{errors} errors[/red]")
    return "  ".join(parts)


def make_all_sessions_table(by_session):
    t = Table(box=box.SIMPLE_HEAD, show_header=True, header_style="bold dim", expand=True, padding=(0,1))
    t.add_column("Session", style="cyan", width=14)
    t.add_column("Tools", width=7, justify="right")
    t.add_column("Errors", width=7, justify="right")
    t.add_column("Duration", width=10, justify="right")
    t.add_column("Last active", width=12, justify="right")

    for sid, evts in sorted(by_session.items(), key=lambda kv: -max(e["ts"] for e in kv[1])):
        calls = pair_calls(evts)
        total = len([c for c in calls if c["pre"]])
        errors = sum(1 for c in calls if c["post"] and is_error(str(c["post"].get("tool_response",""))))
        start = min(e["ts"] for e in evts)
        last = max(e["ts"] for e in evts)
        dur = fmt_dur(last - start)
        ago = fmt_dur(time.time() - last) + " ago"
        t.add_row(
            sid[:12],
            str(total),
            f"[red]{errors}[/red]" if errors else "[dim]0[/dim]",
            dur,
            ago,
        )
    return t


def render(show_all=False, session_filter=None):
    events = load_events()

    if not events:
        hook_path = Path.home() / "claude-viz" / "hook.py"
        return Panel(
            "[dim]No events yet.[/dim]\n\n"
            "Configure hooks in [cyan]~/.claude/settings.json[/cyan]:\n\n"
            f'  [white]python3 {hook_path} PreToolUse[/white]\n'
            f'  [white]python3 {hook_path} PostToolUse[/white]\n\n'
            "Then run Claude Code in another terminal.\n\n"
            f"Run [cyan]python3 ~/claude-viz/viz.py --setup[/cyan] for the full config snippet.",
            title="[bold blue]Claude Workflow Visualizer[/bold blue]",
            border_style="blue",
        )

    by_session = defaultdict(list)
    for e in events:
        by_session[e["session_id"]].append(e)

    if show_all:
        t = make_all_sessions_table(by_session)
        return Panel(
            t,
            title=f"[bold blue]Claude Workflow Visualizer[/bold blue]  [dim]({len(by_session)} sessions)[/dim]",
            border_style="blue",
        )

    # pick session
    if session_filter:
        sid = next((s for s in by_session if s.startswith(session_filter)), None)
        if sid is None:
            sid = max(by_session, key=lambda s: max(e["ts"] for e in by_session[s]))
    else:
        sid = max(by_session, key=lambda s: max(e["ts"] for e in by_session[s]))

    sevents = by_session[sid]
    calls = pair_calls(sevents)
    start_ts = min(e["ts"] for e in sevents)
    stats = make_stats_bar(calls, sid, start_ts)
    table = make_session_table(calls)

    n_sess = len(by_session)
    title = "[bold blue]Claude Workflow Visualizer[/bold blue]"
    if n_sess > 1:
        title += f"  [dim](session {list(by_session.keys()).index(sid)+1}/{n_sess} — use --all to list)[/dim]"

    return Panel(table, title=title, subtitle=stats, border_style="blue")

# ── setup instructions ─────────────────────────────────────────────────────────

SETTINGS_SNIPPET = """\
Add to [cyan]~/.claude/settings.json[/cyan] inside the top-level object:

[bold white]{
  "hooks": {
    "PreToolUse": [{
      "matcher": "",
      "hooks": [{"type": "command", "command": "python3 ~/claude-viz/hook.py PreToolUse"}]
    }],
    "PostToolUse": [{
      "matcher": "",
      "hooks": [{"type": "command", "command": "python3 ~/claude-viz/hook.py PostToolUse"}]
    }],
    "Stop": [{
      "matcher": "",
      "hooks": [{"type": "command", "command": "python3 ~/claude-viz/hook.py Stop"}]
    }]
  }
}[/bold white]

Then in a separate terminal run:
  [cyan]python3 ~/claude-viz/viz.py[/cyan]
"""

# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Claude Workflow Visualizer")
    ap.add_argument("--setup", action="store_true", help="Print hook config instructions")
    ap.add_argument("--all", action="store_true", help="Show all sessions overview")
    ap.add_argument("--session", metavar="ID", help="Filter to a specific session ID prefix")
    ap.add_argument("--once", action="store_true", help="Render once and exit (no live update)")
    ap.add_argument("--clear", action="store_true", help="Clear the events log and exit")
    args = ap.parse_args()

    console = Console()

    if args.setup:
        console.print(Panel(SETTINGS_SNIPPET, title="Setup", border_style="cyan"))
        return

    if args.clear:
        if EVENTS_FILE.exists():
            EVENTS_FILE.unlink()
            console.print(f"[green]Cleared[/green] {EVENTS_FILE}")
        else:
            console.print("[dim]No events file found.[/dim]")
        return

    if args.once:
        console.print(render(show_all=args.all, session_filter=args.session))
        return

    # live mode
    with Live(
        render(show_all=args.all, session_filter=args.session),
        console=console,
        refresh_per_second=REFRESH_HZ,
        screen=False,
        transient=False,
    ) as live:
        try:
            while True:
                live.update(render(show_all=args.all, session_filter=args.session))
                time.sleep(1.0 / REFRESH_HZ)
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
