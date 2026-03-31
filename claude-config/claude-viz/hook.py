#!/usr/bin/env python3
"""
Claude Code hook receiver.
Reads JSON payload from stdin, appends a structured event to ~/.claude-viz/events.jsonl.

Usage (in ~/.claude/settings.json):
  "command": "python3 ~/claude-viz/hook.py PreToolUse"
  "command": "python3 ~/claude-viz/hook.py PostToolUse"
  "command": "python3 ~/claude-viz/hook.py Stop"
"""
import sys
import json
import time
import os

EVENTS_FILE = os.path.expanduser("~/.claude-viz/events.jsonl")
MAX_RESPONSE_LEN = 500  # truncate large tool responses


def main():
    hook_type = sys.argv[1] if len(sys.argv) > 1 else "unknown"

    try:
        payload = json.load(sys.stdin)
    except Exception:
        sys.exit(0)  # never block Claude Code

    resp = payload.get("tool_response", "")
    if isinstance(resp, str) and len(resp) > MAX_RESPONSE_LEN:
        resp = resp[:MAX_RESPONSE_LEN] + "…"

    event = {
        "ts": time.time(),
        "type": hook_type,
        "session_id": payload.get("session_id", ""),
        "tool_name": payload.get("tool_name", ""),
        "tool_input": payload.get("tool_input") or {},
        "tool_response": resp,
    }

    os.makedirs(os.path.dirname(EVENTS_FILE), exist_ok=True)
    with open(EVENTS_FILE, "a") as f:
        f.write(json.dumps(event) + "\n")

    sys.exit(0)  # always exit 0 — never block tools


if __name__ == "__main__":
    main()
