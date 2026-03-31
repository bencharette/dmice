-- DMiceMarkRead command
-- Add this to all three Neovim profile configs:
--   ~/.config/nvim-dmice-local/init.lua
--   ~/.config/nvim-dmice-npx/init.lua
--   ~/.config/nvim-dmice-cobalt/init.lua
--
-- Usage:
--   :DMiceMarkRead        → inserts [READ 2026-03-22 14:30] above cursor line
--   :DMiceMarkRead done   → also moves the block to done.md
--   <leader>dr            → shortcut for :DMiceMarkRead
--   <leader>dR            → shortcut for :DMiceMarkRead done

local function get_timestamp()
  return os.date("[READ %Y-%m-%d %H:%M]")
end

-- Find the inbox file for the current machine (based on which profile is loaded)
local function get_inbox_path()
  local profile = vim.fn.stdpath("config")  -- e.g. ~/.config/nvim-dmice-cobalt
  local machine = profile:match("nvim%-dmice%-(%a+)$") or "local"
  return vim.fn.expand("~/dmice/handoff/inbox-" .. machine .. ".md")
end

local function get_done_path()
  return vim.fn.expand("~/dmice/handoff/done.md")
end

-- Insert a [READ timestamp] line above the current line
local function mark_read_at_cursor()
  local row = vim.api.nvim_win_get_cursor(0)[1]  -- 1-indexed
  local timestamp = get_timestamp()
  vim.api.nvim_buf_set_lines(0, row - 1, row - 1, false, { timestamp })
  vim.notify("Marked as read: " .. timestamp, vim.log.levels.INFO)
end

-- Move the current paragraph block to done.md with a timestamp header
local function mark_read_and_archive()
  local buf = vim.api.nvim_get_current_buf()
  local cursor_row = vim.api.nvim_win_get_cursor(0)[1]  -- 1-indexed
  local lines = vim.api.nvim_buf_get_lines(buf, 0, -1, false)

  -- Find the paragraph boundaries (delimited by blank lines or ---)
  local block_start = cursor_row - 1  -- 0-indexed
  local block_end   = cursor_row - 1  -- 0-indexed

  -- Walk upward to find block start
  while block_start > 0 do
    local line = lines[block_start]  -- 1-indexed into lines table
    if line == "" or line == "---" then
      block_start = block_start + 1
      break
    end
    block_start = block_start - 1
  end

  -- Walk downward to find block end
  while block_end < #lines - 1 do
    local next_line = lines[block_end + 2]  -- +2: 1-indexed + next
    if next_line == "" or next_line == "---" then
      break
    end
    block_end = block_end + 1
  end

  -- Extract the block (0-indexed slice)
  local block = vim.api.nvim_buf_get_lines(buf, block_start, block_end + 1, false)

  if #block == 0 then
    vim.notify("No block found at cursor", vim.log.levels.WARN)
    return
  end

  -- Build done.md entry
  local timestamp = get_timestamp()
  local machine   = (vim.fn.stdpath("config"):match("nvim%-dmice%-(%a+)$") or "local")
  local entry = {
    "",
    "## " .. timestamp .. " — " .. machine,
  }
  for _, l in ipairs(block) do
    table.insert(entry, l)
  end

  -- Append to done.md
  local done_path = get_done_path()
  local done_file = io.open(done_path, "a")
  if not done_file then
    vim.notify("Could not open " .. done_path, vim.log.levels.ERROR)
    return
  end
  done_file:write(table.concat(entry, "\n") .. "\n")
  done_file:close()

  -- Delete block from current buffer
  vim.api.nvim_buf_set_lines(buf, block_start, block_end + 1, false, {})

  vim.notify(
    string.format("Archived %d lines to done.md", #block),
    vim.log.levels.INFO
  )
end

-- ── Commands ────────────────────────────────────────────────────────────────

vim.api.nvim_create_user_command("DMiceMarkRead", function(opts)
  if opts.args == "done" then
    mark_read_and_archive()
  else
    mark_read_at_cursor()
  end
end, {
  nargs = "?",
  desc  = "Mark inbox item as read. Pass 'done' to also archive to done.md",
})

-- ── Keybindings ─────────────────────────────────────────────────────────────

vim.keymap.set("n", "<leader>dr", "<cmd>DMiceMarkRead<cr>",
  { desc = "Mark inbox item as read", silent = true })

vim.keymap.set("n", "<leader>dR", "<cmd>DMiceMarkRead done<cr>",
  { desc = "Mark read + archive to done.md", silent = true })

-- ── Auto-open inbox on launch (optional) ────────────────────────────────────
-- Uncomment if you want the inbox to open automatically when Neovim starts
-- with a DMice profile that has unread messages.
--
-- vim.api.nvim_create_autocmd("VimEnter", {
--   callback = function()
--     local inbox = get_inbox_path()
--     if vim.fn.filereadable(inbox) == 1 then
--       local content = vim.fn.readfile(inbox)
--       -- Only open if there's content below the header (line 4+)
--       if #content > 4 then
--         vim.cmd("edit " .. inbox)
--       end
--     end
--   end,
-- })
