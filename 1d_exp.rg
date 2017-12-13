-- AMR meta programming for all models
import "regent"
local C = regentlib.c

require("model")

-- like stdlib's pow
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
    initializeCells([level_cells[MAX_REFINEMENT_LEVEL]])
    var time : double = 0.0
    while time < T_FINAL - DT do 
      calculateFlux(DX, DT, [level_cells[MAX_REFINEMENT_LEVEL]],
                    [level_faces[MAX_REFINEMENT_LEVEL]])
      applyFlux(DX, DT, [level_cells[MAX_REFINEMENT_LEVEL]],
                    [level_faces[MAX_REFINEMENT_LEVEL]])
      time += DT
      C.printf("time = %f\n",time)
    end
    writeCells([level_cells[MAX_REFINEMENT_LEVEL]])
  end
  return top_level
end

-- main top level task

local main = make_top_level_task()
regentlib.start(main)

