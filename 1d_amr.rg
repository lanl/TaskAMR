-- 1d AMR grid meta programming for all models
import "regent"
local C = regentlib.c

-- implement all required model APIs and link model.rg to your file
require("model")
require("model_amr")
require("1d_make_levels")
require("1d_make_amr")

-- meta programming to create top_level_task
function make_top_level_task()

  -- arrays of region by level
  local meta_region_for_level = terralib.newlist()
  local cell_region_for_level = terralib.newlist()
  local face_region_for_level = terralib.newlist()
  local meta_partition_for_level = terralib.newlist()
  local cell_partition_for_level = terralib.newlist()
  local face_partition_for_level = terralib.newlist()
  local bloated_partition_for_level = terralib.newlist()  -- refactor to "cell"
  local bloated_meta_partition_for_level = terralib.newlist()
  local parent_cell_partition_for_level = terralib.newlist()
  local parent_meta_partition_for_level = terralib.newlist()

  -- array of region and partition declarations
  local declarations = declare_level_regions(meta_region_for_level,
                                             cell_region_for_level,
                                             face_region_for_level,
                                             meta_partition_for_level,
                                             cell_partition_for_level,
                                             face_partition_for_level,
                                             bloated_partition_for_level,
                                             bloated_meta_partition_for_level,
                                             MAX_REFINEMENT_LEVEL,
                                             NUM_PARTITIONS)

  insert_parent_partitions(parent_cell_partition_for_level,
                           parent_meta_partition_for_level)

  local init_parent_partitions = initialize_parent_partitions(cell_partition_for_level,
                                                              cell_region_for_level,
                                                              parent_cell_partition_for_level,
                                                              meta_partition_for_level,
                                                              meta_region_for_level,
                                                              parent_meta_partition_for_level)

  -- meta programming to initialize num_cells per level
  local num_cells = regentlib.newsymbol(int64[MAX_REFINEMENT_LEVEL+1], "num_cells")
  local init_num_cells = make_init_num_cells(num_cells, MAX_REFINEMENT_LEVEL,
                                               cell_region_for_level)

  local init_activity = make_init_activity(meta_region_for_level)

  local write_cells = make_write_cells(num_cells,
                                       meta_partition_for_level,
                                       cell_partition_for_level)

  -- top_level task using previous meta programming
  local task top_level()
    [declarations];
    [init_num_cells];
    [init_parent_partitions];

    var dx : double[MAX_REFINEMENT_LEVEL + 1]
    for level = 1, MAX_REFINEMENT_LEVEL + 1 do
      dx[level] = LENGTH_X / [double]([num_cells][level])
      C.printf("Level %d cells %d dx %e\n", level, [num_cells][level], dx[level])
    end

    [init_activity];
    initializeCells([num_cells][1], [cell_region_for_level[1]])

    __demand(__parallel)
    for color in [cell_partition_for_level[1]].colors do
      calculateGradient(num_cells[1], dx[1], [bloated_partition_for_level[1]][color],
                        [face_partition_for_level[1]][color])
    end

    __demand(__parallel)
    for color in [cell_partition_for_level[1]].colors do
      printFaces([face_partition_for_level[1]][color])
    end

    var needs_regrid = 0
    __demand(__parallel)
    for color in [cell_partition_for_level[1]].colors do
      needs_regrid += flagRegrid([meta_partition_for_level[1]][color],
                                 [face_partition_for_level[1]][color])
    end
    C.printf("Needs regrid = %d\n", needs_regrid)

    if needs_regrid > 0 then

      __demand(__parallel)
      for color in [cell_partition_for_level[1]].colors do
        interpolateToChildren(num_cells[2], [meta_partition_for_level[1]][color],
                              [bloated_partition_for_level[1]][color],
                              [parent_cell_partition_for_level[2]][color])
      end

      __demand(__parallel)
      for color in [meta_partition_for_level[1]].colors do
        updateRefinementBits(num_cells[1]/CELLS_PER_BLOCK_X, [meta_partition_for_level[1]][color],
                             [bloated_meta_partition_for_level[1]][color],
                             [parent_meta_partition_for_level[2]][color])
      end
      fill([meta_region_for_level[1]].needsRefinement, false)

    end -- needs_regrid

    --var time : double = 0.0
    --while time < T_FINAL - DT do 
      --time += DT
      --C.printf("time = %f\n",time)
    --end
    [write_cells];
  end
  return top_level
end

-- top level task

local top_level_task = make_top_level_task()
regentlib.start(top_level_task)

