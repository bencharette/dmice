---
name: DMice Claude onboarding prompt
description: What to paste into a new Claude conversation to give it full DMice project context
type: reference
---

To onboard a new Claude instance on the DMice project, paste the full contents of:

`/home/bench/.claude/projects/-home-bench/memory/project_dmice.md`

at the start of the conversation, prefixed with:

> "Here is context about the project I'm working on:"

This gives the new model:
- Full physics context (DM-Ice pivot fit, coordinate conventions, direction conventions)
- Machine layout (LOCAL / NPX / Cobalt roles)
- Script descriptions and known bugs/fixes
- Current project phase and status
