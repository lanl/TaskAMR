import "regent"
local C = regentlib.c

local CELLS_PER_BLOCK_X = 4
local LEVEL_1_BLOCKS_X = 3
local MAX_REFINEMENT_LEVEL = 3

fspace CellValues
{
  phi : double
}

fspace FaceValues
{
  flux : double
}


-- like stdlib's
function pow(x, y)
  local value = x
  for i = 2,y do
    value = value * x
  end
  return value
end

-- for creating regions for levels 1 to MAX_REFINEMENT_LEVEL
function make_level_regions(n)
  local level_cell = regentlib.newsymbol("level_" .. n .. "_cell")
  local ratio_to_level1 = pow(2, n) / 2
  local num_cells = CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X * ratio_to_level1
  local cells = rquote var [level_cell] = region(ispace(int1d, num_cells), CellValues) end

  local level_face = regentlib.newsymbol("level_" .. n .. "_face")
  local ratio_to_level1 = pow(2, n) / 2
  local num_faces = num_cells + 1
  local faces = rquote var [level_face] = region(ispace(int1d, num_faces), FaceValues) end
  return level_cell, cells, level_face, faces
end

function make_top_level_task()
  local level_cells = terralib.newlist()
  local level_faces = terralib.newlist()
  local declare_level_regions = terralib.newlist()
  for n = 1, MAX_REFINEMENT_LEVEL do
    local level_cell, declare_cells, level_face, declare_faces = make_level_regions(n)
    level_cells:insert(level_cell)
    declare_level_regions:insert(declare_cells)
    level_faces:insert(level_face)
    declare_level_regions:insert(declare_faces)
  end
  local loops = terralib.newlist()
  for n = 1, MAX_REFINEMENT_LEVEL do
    loops:insert(rquote
      var total : int64 = 0
      for cell in [level_cells[n]] do
        total = total + 1
      end
      C.printf(["level " .. n .. " cells %d \n"], total)
    end)
  end
  local task top_level()
    [declare_level_regions];
    [loops];
    initialize_cells([level_cells[MAX_REFINEMENT_LEVEL]])
    print_cells([level_cells[MAX_REFINEMENT_LEVEL]])
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

