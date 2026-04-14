-- DMice Commands Plugin
-- Provides: DMiceMemory, DMiceInbox, DMiceSend, DMiceSyncPull, DMiceSyncPush, DMiceSyncAll, DMiceDashboard

local function get_machine()
  return vim.fn.stdpath("config"):match("nvim%-dmice%-(%a+)$") or "local"
end

local function get_memory_path(machine)
  return vim.fn.expand("~/dmice/memory/" .. (machine or get_machine()) .. ".md")
end

local function get_inbox_path(machine)
  return vim.fn.expand("~/dmice/handoff/inbox-" .. (machine or get_machine()) .. ".md")
end

local function get_dmice_root()
  return vim.fn.expand("~/dmice")
end

-- ── :DMiceMemory ─────────────────────────────────────────────────────────────

vim.api.nvim_create_user_command("DMiceMemory", function()
  vim.cmd("edit " .. get_memory_path())
end, { desc = "Open memory file for current machine" })

-- ── :DMiceInbox ──────────────────────────────────────────────────────────────

vim.api.nvim_create_user_command("DMiceInbox", function()
  vim.cmd("edit " .. get_inbox_path())
end, { desc = "Open inbox file for current machine" })

-- ── :DMiceSend <machine> ─────────────────────────────────────────────────────

vim.api.nvim_create_user_command("DMiceSend", function(opts)
  local target = opts.args
  if target == "" then
    vim.notify("Usage: :DMiceSend <machine>  (local | npx | cobalt)", vim.log.levels.WARN)
    return
  end
  local path = get_inbox_path(target)
  local timestamp = os.date("[%Y-%m-%d %H:%M]")
  local from = get_machine()
  local f = io.open(path, "a")
  if not f then
    vim.notify("Could not open " .. path, vim.log.levels.ERROR)
    return
  end
  f:write("\n---\n" .. timestamp .. " from " .. from .. ":\n")
  f:close()
  vim.cmd("edit " .. path)
  vim.cmd("normal! G")
  vim.cmd("startinsert!")
end, {
  nargs = 1,
  complete = function() return { "local", "npx", "cobalt" } end,
  desc = "Append timestamped message header to inbox-<machine>.md and open it",
})

-- ── Sync helpers ─────────────────────────────────────────────────────────────

local function git(args, on_done)
  local root = get_dmice_root()
  vim.fn.jobstart("git -C " .. root .. " " .. args, {
    on_exit = function(_, code)
      local verb = args:match("^%S+")
      if code == 0 then
        vim.notify("git " .. verb .. " OK", vim.log.levels.INFO)
      else
        vim.notify("git " .. verb .. " failed (exit " .. code .. ")", vim.log.levels.WARN)
      end
      if on_done then on_done(code) end
    end,
  })
end

local function do_push()
  local root  = get_dmice_root()
  local msg   = "handoff sync from " .. get_machine() .. " at " .. os.date("%Y-%m-%d %H:%M")
  local cmd   = string.format(
    "git -C %s add handoff/ memory/ && "
    .. "git -C %s diff --cached --quiet || git -C %s commit -m %q && "
    .. "git -C %s push",
    root, root, root, msg, root
  )
  vim.fn.jobstart(cmd, {
    on_exit = function(_, code)
      if code == 0 then
        vim.notify("Handoff sync pushed", vim.log.levels.INFO)
      else
        vim.notify("Push failed (exit " .. code .. ")", vim.log.levels.WARN)
      end
    end,
  })
end

-- ── :DMiceSyncPull ───────────────────────────────────────────────────────────

vim.api.nvim_create_user_command("DMiceSyncPull", function()
  git("pull --rebase")
end, { desc = "Pull latest handoff files from remote" })

-- ── :DMiceSyncPush ───────────────────────────────────────────────────────────

vim.api.nvim_create_user_command("DMiceSyncPush", function()
  do_push()
end, { desc = "Commit and push handoff + memory files" })

-- ── :DMiceSyncAll ────────────────────────────────────────────────────────────

vim.api.nvim_create_user_command("DMiceSyncAll", function()
  git("pull --rebase", function(code)
    if code ~= 0 then return end
    do_push()
    vim.defer_fn(function() vim.cmd("DMiceInbox") end, 500)
  end)
end, { desc = "Pull, push, then open inbox" })

-- ── :DMiceDashboard ──────────────────────────────────────────────────────────

vim.api.nvim_create_user_command("DMiceDashboard", function()
  vim.cmd("edit " .. vim.fn.expand("~/dmice/COMMANDS.md"))
end, { desc = "Open DMice command reference dashboard" })

-- ── Keybindings ──────────────────────────────────────────────────────────────

vim.keymap.set("n", "<leader>dm", "<cmd>DMiceMemory<cr>",    { desc = "Open machine memory",   silent = true })
vim.keymap.set("n", "<leader>di", "<cmd>DMiceInbox<cr>",     { desc = "Open machine inbox",    silent = true })
vim.keymap.set("n", "<leader>da", "<cmd>DMiceSyncAll<cr>",   { desc = "Sync all + open inbox", silent = true })
vim.keymap.set("n", "<leader>dd", "<cmd>DMiceDashboard<cr>", { desc = "Open DMice dashboard",  silent = true })
