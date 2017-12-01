import "regent"
local C = regentlib.c

local BLOCK_X = 4
local LEVEL_1_NX_BLOCKS = 3
local MAX_LEVEL = 3

fspace CellValues
{
  phi : double
}

-- like stdlib's
function pow(x, y)
  local value = x
  for i = 2,y do
    value = value * x
  end
  return value
end

-- for creating regions for levels 1 to MAX_LEVEL
function make_level_region(n)
  local level_region = regentlib.newsymbol("level_" .. n .. "_region")
  local ratio_to_level1 = pow(2, n) / 2
  local num_cells = BLOCK_X * LEVEL_1_NX_BLOCKS * ratio_to_level1
  local value = rquote var [level_region] = region(ispace(int1d, num_cells), CellValues) end
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
      var total : int64 = 0
      for cell in [level_regions[n]] do
        total = total + 1
      end
      C.printf(["level " .. n .. " cells %d \n"], total)
    end)
  end
  local task top_level()
    [values];
    [loops];
    initialize_cells([level_regions[MAX_LEVEL]])
    print_cells([level_regions[MAX_LEVEL]])
  end
  return top_level
end

task initialize_cells(cell_region: region(ispace(int1d), CellValues))
where
  writes(cell_region.phi)
do
  var size : int64 = cell_region.ispace.bounds.hi - cell_region.ispace.bounds.lo + 1
  for cell in cell_region.ispace do
    if [int64](cell) < (size/2) then
      cell_region[cell].phi = 1.0
    else
      cell_region[cell].phi = 0.0
    end
  end
  C.printf("initialize_cells %d cells\n", size)
end

task print_cells(cell_region: region(ispace(int1d), CellValues))
where
  reads(cell_region.phi)
do
  C.printf("phi: ")
  for cell in cell_region do
    C.printf("%f ", cell_region[cell].phi)
  end
  C.printf("\n")
end

local main = make_top_level_task()
main:printpretty()
regentlib.start(main)

