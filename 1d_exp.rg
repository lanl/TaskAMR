import "regent"
local C = regentlib.c

local BLOCK_X = 4
local LEVEL_1_NX_BLOCKS = 3
local MAX_LEVEL = 3

fspace CellValues
{
  phi : double
}

function make_level_region(n)
  local level_region = regentlib.newsymbol("level_" .. n .. "_region")
  local value = rquote var [level_region] = region(ispace(int1d, BLOCK_X * LEVEL_1_NX_BLOCKS), CellValues) end
  return level_region, value
end



function make_top_level_task()
  local level_regions = terralib.newlist()
  local values = terralib.newlist()
  for n = 1, MAX_LEVEL do
    local level_region, value = make_level_region(n)
    level_regions:insert(level_region)
    values:insert(value)
  end
  local loops = terralib.newlist()
  for n = 1, MAX_LEVEL do
    loops:insert(rquote
      for cell in [level_regions[n]] do
        C.printf(["level " .. n .. " cell %d\n"], cell)
      end
    end)
  end
  local task top_level()
    [values];
    [loops]
  end
  return top_level
end

local main = make_top_level_task()
main:printpretty()
regentlib.start(main)

