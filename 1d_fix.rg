--Copyright (c) 2018, Los Alamos National Security, LLC
--All rights reserved.

--Copyright 2018. Los Alamos National Security, LLC. This software was produced under
--U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL),
--which is operated by Los Alamos National Security, LLC for the U.S. Department of
--Energy. The U.S. Government has rights to use, reproduce, and distribute this
--software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
--WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
--If software is modified to produce derivative works, such modified software should be
--clearly marked, so as not to confuse it with the version available from LANL.

-- 1d fixed grid meta programming for all models
import "regent"
local C = regentlib.c

-- implement all required model APIs and link model.rg to your file
require("model")
require("1d_make_levels")

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
  local bloated_meta_partition_for_level = terralib.newlist()

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

  -- meta programming to initialize num_cells per level
  local num_cells = regentlib.newsymbol(int64[MAX_REFINEMENT_LEVEL+1], "num_cells")
  local dx = regentlib.newsymbol(double[MAX_REFINEMENT_LEVEL+1], "dx")
  local needs_regrid = regentlib.newsymbol(int64[MAX_REFINEMENT_LEVEL+1], "needs_regrid")
  local init_num_cells = make_init_num_cells(num_cells,
                                             dx,
                                             needs_regrid,
                                             MAX_REFINEMENT_LEVEL,
                                             cell_region_for_level)

  -- top_level task using previous meta programming
  local task top_level()
    [declarations];
    [init_num_cells];

    for level = 1, MAX_REFINEMENT_LEVEL + 1 do
      [dx][level] = LENGTH_X / [double]([num_cells][level])
      C.printf("Level %d cells %d dx %e\n", level, [num_cells][level], [dx][level])
    end

    fill([meta_region_for_level[MAX_REFINEMENT_LEVEL]].isActive, true)

    __demand(__parallel)
    for color in [cell_partition_for_level[MAX_REFINEMENT_LEVEL]].colors do
      initializeCells([num_cells][MAX_REFINEMENT_LEVEL],
                      [cell_partition_for_level[MAX_REFINEMENT_LEVEL]][color])
    end

    var time : double = 0.0
    while time < T_FINAL - DT do

      __demand(__parallel)
      for color in [cell_partition_for_level[MAX_REFINEMENT_LEVEL]].colors do
        calculateFlux(num_cells[MAX_REFINEMENT_LEVEL], [dx][MAX_REFINEMENT_LEVEL], DT,
                      [meta_partition_for_level[MAX_REFINEMENT_LEVEL]][color],
                      [bloated_partition_for_level[MAX_REFINEMENT_LEVEL]][color],
                      [face_partition_for_level[MAX_REFINEMENT_LEVEL]][color])
      end

      __demand(__parallel)
      for color in [cell_partition_for_level[MAX_REFINEMENT_LEVEL]].colors do
          applyFlux([dx][MAX_REFINEMENT_LEVEL], DT,
                    [meta_partition_for_level[MAX_REFINEMENT_LEVEL]][color],
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

