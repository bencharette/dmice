# How to Create a `/note` Slash Command in Claude Code

This guide explains how to create a custom `/note` slash command that saves facts, decisions, or context to Claude's persistent memory for future sessions.

---

## What Is a Slash Command?

In Claude Code (the CLI), slash commands are shortcuts that invoke pre-written prompts. When you type `/note`, Claude expands it into a full instruction and executes it. Commands are defined as Markdown files in a special directory.

---

## Where to Put Command Files

Claude Code looks for custom slash commands in:

```
~/.claude/commands/          # personal commands (available in any project)
.claude/commands/            # project-scoped commands (available in this repo only)
```

Each `.md` file in those directories becomes a slash command named after the file.

---

## Step 1 — Create the Commands Directory

```bash
mkdir -p ~/.claude/commands
```

---

## Step 2 — Create the `note.md` File

```bash
cat > ~/.claude/commands/note.md << 'EOF'
Save the following as a persistent memory for future Claude sessions.

Determine the appropriate memory type:
- **user** — facts about the user's role, preferences, or expertise
- **feedback** — corrections or guidance on how Claude should behave
- **project** — ongoing work context, decisions, deadlines
- **reference** — pointers to external systems or resources

Write the memory to `/home/YOUR_USERNAME/.claude/projects/YOUR_PROJECT_PATH/memory/` as a `.md` file with this frontmatter:

```
---
name: <short descriptive name>
description: <one-line description>
type: <user|feedback|project|reference>
---

<memory content>
```

Then add a pointer line to `MEMORY.md` in the same directory.

Note to save: $ARGUMENTS
EOF
```

Replace `YOUR_USERNAME` and `YOUR_PROJECT_PATH` with your actual paths (e.g., `/home/bench/.claude/projects/-home-bench/memory/`).

---

## Step 3 — Use It

In any Claude Code session:

```
/note the pipeline now runs on NPX, not locally
```

Claude will determine the memory type, write a `.md` file, and update `MEMORY.md`.

You can also pass multi-line notes:

```
/note user prefers terse responses with no trailing summaries
```

---

## How the DMice Project Uses Memory

The DMice project stores memory at:
```
~/.claude/projects/-home-bench/memory/
```

Files:
| File | Type | Contents |
|------|------|----------|
| `project_dmice.md` | project | Full project context — machines, scripts, physics, known bugs |
| `user_role.md` | user | User is an IceCube physicist working on DM-Ice |
| `feedback_comparison_script.md` | feedback | Always use `sim_linefit_comparison.py` from `/home/bench/dmice/` |
| `reference_github.md` | reference | GitHub repo URLs |
| `reference_onboarding.md` | reference | How to onboard a new Claude instance |
| `MEMORY.md` | index | Index of all memory files (loaded every session) |

To onboard a new Claude instance, paste `project_dmice.md` at the start of the conversation — or rely on `MEMORY.md` being auto-loaded.

---

## Tip: Project-Scoped Commands

If you want `/note` only for this project, put it in the repo instead:

```bash
mkdir -p .claude/commands
# copy note.md here instead of ~/.claude/commands/
```

This keeps the command version-controlled alongside the project.

---

## Further Reading

- Claude Code docs: `claude --help` or `/help` inside a session
- Memory system: see `~/.claude/projects/.../memory/MEMORY.md`
- Command files support `$ARGUMENTS` for text passed after the command name
