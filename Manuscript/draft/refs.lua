local main = "acedisparity_draft"   -- main file name, no extension
local labels = {}

local words = { figure = "Figure", table = "Table", section = "Section",
                subsection = "Section", subsubsection = "Section",
                equation = "Equation", appendix = "Appendix" }
local keyword = { fig = "Figure", tab = "Table", sec = "Section", eq = "Equation" }

local function load_aux(name, prefix)
  local f = io.open(name .. ".aux")
  if not f then
    io.stderr:write("refs.lua: cannot open " .. name .. ".aux\n")
    return
  end
  for line in f:lines() do
    local k, v = line:match("^\\newlabel{([^}]*)}{{([^}]*)}")
    if k and not k:match("@cref$") then
      -- hyperref stores the type, e.g. {figure.1}
      local typ = line:match("{(%a+)%.[^{}]*}{}}%s*$")
      local word = typ and words[typ]
      if not word then
        word = keyword[k:match("^(%a+)[:%-_]") or ""]
      end
      local info = { num = v, word = word }
      labels[prefix .. k] = info
      if prefix ~= "" and labels[k] == nil then labels[k] = info end
    end
  end
  f:close()
end

load_aux(main, "")

local f = io.open(main .. ".tex")
if f then
  local src = f:read("a")
  f:close()
  for pre, name in src:gmatch("\\externaldocument%[([^%]]*)%]{([^}]*)}") do
    load_aux(name, pre)
  end
  for name in src:gmatch("\\externaldocument{([^}]*)}") do
    load_aux(name, "")
  end
end

local figure_count = 0
local table_count = 0

local function caption_starts_with_figure(caption)
  local blocks = caption.long or caption
  local text = pandoc.utils.stringify(blocks)
  return text:match("^Figure%s+") ~= nil
end

local function caption_starts_with_table(caption)
  local blocks = caption.long or caption
  local text = pandoc.utils.stringify(blocks)
  return text:match("^Table%s+") ~= nil
end

local function prepend_caption_prefix(caption, word, number)
  local blocks = caption.long or caption
  local first_block = blocks[1]

  if first_block and (first_block.t == "Para" or first_block.t == "Plain") then
    local prefix = {
      pandoc.Str(word),
      pandoc.Space(),
      pandoc.Str(number),
      pandoc.Str("."),
      pandoc.Space()
    }

    for i = #prefix, 1, -1 do
      table.insert(first_block.content, 1, prefix[i])
    end
  end
end

function Figure(el)
  figure_count = figure_count + 1

  local figure_number = tostring(figure_count)
  if el.identifier and labels[el.identifier] then
    figure_number = labels[el.identifier].num
  end

  if el.caption and not caption_starts_with_figure(el.caption) then
    prepend_caption_prefix(el.caption, "Figure", figure_number)
  end

  return el
end

function Table(el)
  table_count = table_count + 1

  local table_number = tostring(table_count)
  if el.identifier and labels[el.identifier] then
    table_number = labels[el.identifier].num
  end

  if el.caption and not caption_starts_with_table(el.caption) then
    prepend_caption_prefix(el.caption, "Table", table_number)
  end

  return el
end

function Inlines(inlines)
  for i, el in ipairs(inlines) do
    if el.t == "Link" and el.attributes["reference"] then
      local ref = el.attributes["reference"]
      local info = labels[ref]
      if info then
        local before = ""
        if i > 2 then
          before = pandoc.utils.stringify(pandoc.Inlines{inlines[i-2], inlines[i-1]})
        elseif i > 1 then
          before = pandoc.utils.stringify(inlines[i-1])
        end
        local has_word = before:lower():match("fig%.?[%s\u{a0}]*$")
          or before:lower():match("figure[%s\u{a0}]*$")
          or before:lower():match("table[%s\u{a0}]*$")
          or before:lower():match("sec%.?[%s\u{a0}]*$")
          or before:lower():match("section[%s\u{a0}]*$")
          or before:lower():match("eq%.?[%s\u{a0}]*$")
          or before:lower():match("equation[%s\u{a0}]*$")
        if info.word and not has_word then
          inlines[i] = pandoc.Str(info.word .. "\u{a0}" .. info.num)
        else
          inlines[i] = pandoc.Str(info.num)
        end
      else
        io.stderr:write("refs.lua: unresolved label: " .. ref .. "\n")
      end
    end
  end
  return inlines
end


-- pandoc acedisparity_draft.tex -o acedisparity.docx   --citeproc --bibliography=references.bib   --number-sections --lua-filter=refs.lua
