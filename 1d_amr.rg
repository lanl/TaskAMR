-- 1d fixed grid meta programming for all models
import "regent"
local C = regentlib.c

-- implement all required model APIs and link model.rg to your file
require("model")

-- keeps track of refinement state of the binary tree
-- might be more performant if a bitfield was used
fspace RefinementBits
{
  isRefined: bool,
  needsRefinement: bool,
  wantsCoarsening: bool,
  plusXMoreRefined: bool,
  minusXMoreRefined: bool,
  plusXMoreCoarse: bool,
  minusXMoreCoarse: bool
}

-- like stdlib's pow
function pow(x, y)
  local value = x
  for i = 2,y do
    value = value * x
  end
  return value
end

-- meta programming to create regions for levels 1 to MAX_REFINEMENT_LEVEL
function make_level_regions(n)
  local ratio_to_level1 = pow(2, n) / 2

  local cell_region = regentlib.newsymbol("level_" .. n .. "_cell_region")
  local num_cells = CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X * ratio_to_level1
  local cell_declaration =
    rquote var [cell_region] = region(ispace(int1d, num_cells), CellValues) end

  local face_region = regentlib.newsymbol("level_" .. n .. "_face_region")
  local num_faces = num_cells + 1
  local face_declaration =
    rquote var [face_region] = region(ispace(int1d, num_faces), FaceValues) end

  local refinement_region = regentlib.newsymbol("level_" .. n .. "_refinement_region")
  local num_refinements = LEVEL_1_BLOCKS_X * ratio_to_level1
  local refinement_declaration =
    rquote var [refinement_region] = region(ispace(int1d, num_refinements), RefinementBits) end


  return cell_region, cell_declaration, face_region, face_declaration, refinement_region,
    refinement_declaration
end

-- meta programming to create top_level_task
function make_top_level_task()
  -- arrays of region by level
  local cell_region_for_level = terralib.newlist()
  local face_region_for_level = terralib.newlist()
  local refinement_region_for_level = terralib.newlist()
  -- array of region declarations
  local declare_level_regions = terralib.newlist()

  for n = 1, MAX_REFINEMENT_LEVEL do
    local cell_region, declare_cells, face_region, declare_faces, refinement_region,
      declare_refinements = make_level_regions(n)
    cell_region_for_level:insert(cell_region)
    declare_level_regions:insert(declare_cells)
    face_region_for_level:insert(face_region)
    declare_level_regions:insert(declare_faces)
    refinement_region_for_level:insert(refinement_region)
    declare_level_regions:insert(declare_refinements)
  end

  -- meta programming to create loop over levels
  local loops = terralib.newlist()
  for n = 1, MAX_REFINEMENT_LEVEL do
    loops:insert(rquote
      var total : int64 = 0
      for cell in [cell_region_for_level[n]] do
        total = total + 1
      end
      C.printf(["level " .. n .. " cells %d \n"], total)
    end)
  end

  -- create top_level task with previous meta programming
  local task top_level()
    [declare_level_regions];
    [loops];

    initializeCells([cell_region_for_level[1]])

    var time : double = 0.0
    while time < T_FINAL - DT do 
      calculateFlux(DX, DT, [cell_region_for_level[MAX_REFINEMENT_LEVEL]],
                    [face_region_for_level[MAX_REFINEMENT_LEVEL]])
      applyFlux(DX, DT, [cell_region_for_level[MAX_REFINEMENT_LEVEL]],
                    [face_region_for_level[MAX_REFINEMENT_LEVEL]])
      time += DT
      C.printf("time = %f\n",time)
    end
    writeCells([cell_region_for_level[1]])
  end
  return top_level
end

-- top level task

local top_level_task = make_top_level_task()
regentlib.start(top_level_task)
