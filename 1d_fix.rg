-- 1d fixed grid meta programming for all models
import "regent"
local C = regentlib.c

-- implement all required model APIs and link model.rg to your file
require("model")
require("1d_make_level_regions")

-- meta programming to create top_level_task
function make_top_level_task()

  -- arrays of region by level
  local meta_region_for_level = terralib.newlist()
  local cell_region_for_level = terralib.newlist()
  local face_region_for_level = terralib.newlist()
  local meta_partition_for_level = terralib.newlist()
  local cell_partition_for_level = terralib.newlist()
  local face_partition_for_level = terralib.newlist()
  local bloated_partition_for_level = terralib.newlist()

  -- array of region declarations
  local declarations = terralib.newlist()

  for n = 1, MAX_REFINEMENT_LEVEL do
    local cell_region, declare_cells, face_region, declare_faces, cell_partition, declare_cpart,
      face_partition, declare_fpart, bloated_partition, declare_bpart, meta_region,
      declare_meta, meta_partition, declare_mpart = make_level_regions(n, NUM_PARTITIONS)
    meta_region_for_level:insert(meta_region)
    declarations:insert(declare_meta)
    meta_partition_for_level:insert(meta_partition)
    declarations:insert(declare_mpart)
    cell_region_for_level:insert(cell_region)
    declarations:insert(declare_cells)
    face_region_for_level:insert(face_region)
    declarations:insert(declare_faces)
    cell_partition_for_level:insert(cell_partition)
    declarations:insert(declare_cpart)
    face_partition_for_level:insert(face_partition)
    declarations:insert(declare_fpart)
    bloated_partition_for_level:insert(bloated_partition)
    declarations:insert(declare_bpart)
  end

  local num_cells = regentlib.newsymbol(int64[MAX_REFINEMENT_LEVEL+1], "num_cells")

  -- meta programming to create loop over levels
  local init_num_cells = terralib.newlist()
  init_num_cells:insert(rquote var [num_cells] end)
  for n = 1, MAX_REFINEMENT_LEVEL do
    init_num_cells:insert(rquote
      var limits = [cell_region_for_level[n]].ispace.bounds
      var total : int64 = limits.hi - limits.lo + 1
      [num_cells][n] = total
    end)
  end

  -- top_level task using previous meta programming
  local task top_level()
    [declarations];
    [init_num_cells];

    var dx : double[MAX_REFINEMENT_LEVEL + 1]
    for level = 1, MAX_REFINEMENT_LEVEL + 1 do
      dx[level] = LENGTH_X / [double]([num_cells][level])
      C.printf("Level %d cells %d dx %e\n", level, [num_cells][level], dx[level])
    end

    initializeCells([cell_region_for_level[MAX_REFINEMENT_LEVEL]])
    var time : double = 0.0
    while time < T_FINAL - DT do

      __demand(__parallel)
      for color in [cell_partition_for_level[MAX_REFINEMENT_LEVEL]].colors do
        calculateFlux(num_cells[MAX_REFINEMENT_LEVEL], dx[MAX_REFINEMENT_LEVEL], DT,
        [bloated_partition_for_level[MAX_REFINEMENT_LEVEL]][color],
        [face_partition_for_level[MAX_REFINEMENT_LEVEL]][color])
      end

      __demand(__parallel)
      for color in [cell_partition_for_level[MAX_REFINEMENT_LEVEL]].colors do
          applyFlux(dx[MAX_REFINEMENT_LEVEL], DT,
          [cell_partition_for_level[MAX_REFINEMENT_LEVEL]][color],
          [face_partition_for_level[MAX_REFINEMENT_LEVEL]][color])
      end

      time += DT
      C.printf("time = %f\n",time)
    end
    writeCells([num_cells][MAX_REFINEMENT_LEVEL], [cell_region_for_level[MAX_REFINEMENT_LEVEL]])
  end
  return top_level
end

-- top level task

local top_level_task = make_top_level_task()
regentlib.start(top_level_task)

