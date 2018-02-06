-- 1d fixed grid meta programming for all models
import "regent"
local C = regentlib.c

MAX_REFINEMENT_LEVEL = 3
NUM_PARTITIONS = 7
CELLS_PER_BLOCK_X = 8
LEVEL_1_BLOCKS_X = 11

fspace CellValues
{
}

fspace FaceValues
{
}

require("linear_advection")
require("linear_advection_amr")
require("1d_make_levels")
require("1d_make_amr")

__demand(__inline)
task ASSERT_BOOL_EQUAL(actual : bool,
                        expected : bool,
                        descriptor : &int8)

  if actual == expected then
    C.printf("%s: \x1b[32mPASS\x1b[0m %d = %d\n", descriptor, actual, expected)
  else
    C.printf("%s: \x1b[31mFAIL\x1b[0m %d != %d\n", descriptor, actual, expected)
    C.exit(1)
  end

end

__demand(__inline)
task ASSERT_INT64_EQUAL(actual : int64,
                        expected : int64,
                        descriptor : &int8)

  if actual == expected then
    C.printf("%s: \x1b[32mPASS\x1b[0m %d = %d\n", descriptor, actual, expected)
  else
    C.printf("%s: \x1b[31mFAIL\x1b[0m %d != %d\n", descriptor, actual, expected)
    C.exit(1)
  end

end

-- meta programming to create top_level_task
function make_test_declarations()

  -- arrays of region by level
  local meta_region_for_level = terralib.newlist()
  local cell_region_for_level = terralib.newlist()
  local face_region_for_level = terralib.newlist()
  local meta_partition_for_level = terralib.newlist()
  local cell_partition_for_level = terralib.newlist()
  local face_partition_for_level = terralib.newlist()
  local bloated_partition_for_level = terralib.newlist()
  local bloated_meta_partition_for_level = terralib.newlist()
  local parent_cell_partition_for_level = terralib.newlist()
  local parent_meta_partition_for_level = terralib.newlist()
  local bloated_parent_meta_partition_for_level = terralib.newlist()
  local bloated_cell_partition_by_parent_for_level = terralib.newlist()

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
  local level_needs_regrid = regentlib.newsymbol(int64[MAX_REFINEMENT_LEVEL+1], "level_needs_regrid")
  local needs_regrid = regentlib.newsymbol(int64, "needs_regrid")
  local init_num_cells = make_init_num_cells(num_cells,
                                             dx,
                                             level_needs_regrid,
                                             MAX_REFINEMENT_LEVEL,
                                             cell_region_for_level)

  local init_activity = make_init_activity(meta_region_for_level)

  insert_parent_partitions(parent_cell_partition_for_level,
                           parent_meta_partition_for_level,
                           bloated_parent_meta_partition_for_level,
                           bloated_cell_partition_by_parent_for_level)

  local init_parent_partitions = initialize_parent_partitions(cell_partition_for_level,
                                                              cell_region_for_level,
                                                              parent_cell_partition_for_level,
                                                              meta_partition_for_level,
                                                              meta_region_for_level,
                                                              parent_meta_partition_for_level,
                                                              bloated_parent_meta_partition_for_level,
                                                              bloated_cell_partition_by_parent_for_level,
                                                              num_cells)
 
  local do_regrid = make_do_regrid(num_cells,
                                   meta_region_for_level,
                                   cell_region_for_level,
                                   meta_partition_for_level,
                                   cell_partition_for_level,
                                   parent_cell_partition_for_level,
                                   parent_meta_partition_for_level,
                                   bloated_partition_for_level,
                                   bloated_cell_partition_by_parent_for_level,
                                   bloated_parent_meta_partition_for_level,
                                   bloated_meta_partition_for_level
                                   )

  -- top_level task using previous meta programming
  local task test_declarations()
    [declarations];
    [init_num_cells];
    [init_parent_partitions];
    [init_activity];

    var result : bool;
    var test_name : &int8;
    test_name = [&int8](C.malloc(81));


    -- TEST do_regrid coarsening

    -- looping instead of Lua-recompiling [do_regrid] saves 12 minutes of time
    for test = 1, 7 do

      if test == 1 then
        -- test Level 1 does not coarsen
        fill([meta_region_for_level[1]].wantsCoarsening, true);
        C.sprintf(test_name, "do_regrid::Level 1 does not coarsen");
      elseif test == 2 then
        -- test only coarsen if sibling does
        [meta_region_for_level[1]][0].isActive = false;
        [meta_region_for_level[1]][0].isRefined = true;
        [meta_region_for_level[2]][0].isActive = true;
        [meta_region_for_level[2]][1].isActive = true;
        [meta_region_for_level[2]][1].wantsCoarsening = true;
        C.sprintf(test_name, "do_regrid::Only coarsen if sibling also wants");
      elseif test == 3 then
        -- test can't coarsen if neighbor refines
        fill([meta_region_for_level[1]].isRefined, true);
        fill([meta_region_for_level[1]].isActive, false);
        fill([meta_region_for_level[2]].isActive, true);
        [meta_region_for_level[2]][0].wantsCoarsening = true;
        [meta_region_for_level[2]][1].wantsCoarsening = true;
        [meta_region_for_level[2]][2].needsRefinement = true;
        C.sprintf(test_name, "do_regrid::Cannot coarsen if right neighbor refines");
      elseif test == 4 then
        fill([meta_region_for_level[1]].isRefined, true);
        fill([meta_region_for_level[1]].isActive, false);
        fill([meta_region_for_level[2]].isActive, true);
        [meta_region_for_level[2]][2].wantsCoarsening = true;
        [meta_region_for_level[2]][3].wantsCoarsening = true;
        [meta_region_for_level[2]][1].needsRefinement = true;
        C.sprintf(test_name, "do_regrid::Cannot coarsen if left neighbor refines");
      elseif test == 5 then
        -- test can't coarsen if neighbor is refined
        fill([meta_region_for_level[1]].isRefined, true);
        fill([meta_region_for_level[1]].isActive, false);
        [meta_region_for_level[1]][0].isActive = true;
        [meta_region_for_level[1]][0].isRefined = false;
        fill([meta_region_for_level[2]].isActive, true);
        fill([meta_region_for_level[2]].isRefined, false);
        [meta_region_for_level[2]][0].isActive = false;
        [meta_region_for_level[2]][1].isActive = false;
        [meta_region_for_level[2]][4].isActive = false;
        [meta_region_for_level[2]][4].isRefined = true;
        [meta_region_for_level[3]][8].isActive = true;
        [meta_region_for_level[3]][9].isActive = true;
        [meta_region_for_level[2]][1].isActive = false;
        [meta_region_for_level[2]][3].wantsCoarsening = true;
        [meta_region_for_level[2]][2].wantsCoarsening = true;
        C.sprintf(test_name, "do_regrid::Cannot coarsen if right neighbor is refined");
      elseif test == 6 then
        fill([meta_region_for_level[1]].isRefined, true);
        fill([meta_region_for_level[1]].isActive, false);
        fill([meta_region_for_level[2]].isActive, true);
        fill([meta_region_for_level[2]].isRefined, false);
        [meta_region_for_level[2]][1].isActive = false;
        [meta_region_for_level[2]][1].isRefined = true;
        [meta_region_for_level[3]][2].isActive = false;
        [meta_region_for_level[3]][3].isActive = false;
        [meta_region_for_level[2]][3].wantsCoarsening = true;
        [meta_region_for_level[2]][2].wantsCoarsening = true;
        C.sprintf(test_name, "do_regrid::Cannot coarsen if left neighbor is refined");
      end

      [do_regrid];

      if test == 1 then
        result = [meta_region_for_level[1]][0].isActive;
      elseif test == 2 then
        var not_sibling : bool = (not [meta_region_for_level[1]][0].isActive)
                                 and [meta_region_for_level[1]][0].isRefined
                                 and [meta_region_for_level[2]][0].isActive
                                 and [meta_region_for_level[2]][1].isActive;
        [meta_region_for_level[2]][0].wantsCoarsening = true;
        [do_regrid];
        result = not_sibling
                and [meta_region_for_level[1]][0].isActive
                and (not [meta_region_for_level[1]][0].isRefined)
                and (not [meta_region_for_level[2]][0].isActive)
                and (not [meta_region_for_level[2]][1].isActive)
      elseif test == 3 then
        result = [meta_region_for_level[1]][0].isRefined
                 and (not [meta_region_for_level[1]][0].isActive)
                 and [meta_region_for_level[2]][0].isActive
                 and [meta_region_for_level[2]][1].isActive
      elseif test == 4 then
        result = [meta_region_for_level[1]][1].isRefined
                 and (not [meta_region_for_level[1]][1].isActive)
                 and [meta_region_for_level[2]][2].isActive
                 and [meta_region_for_level[2]][3].isActive
      elseif test == 5 then
        result = [meta_region_for_level[1]][1].isRefined
                 and (not [meta_region_for_level[1]][1].isActive)
                 and [meta_region_for_level[2]][2].isActive
                 and [meta_region_for_level[2]][3].isActive
      elseif test == 6 then
        result = [meta_region_for_level[1]][1].isRefined
                 and (not [meta_region_for_level[1]][1].isActive)
                 and [meta_region_for_level[2]][2].isActive
                 and [meta_region_for_level[2]][3].isActive
      end

      ASSERT_BOOL_EQUAL(result, true, test_name);
    end -- test

    -- test cell count

    var ncells_level1 : int64 = [cell_region_for_level[1]].ispace.bounds.hi
                              - [cell_region_for_level[1]].ispace.bounds.lo + 1
    C.sprintf(test_name, "declarations::Level 1 cell count")
    ASSERT_INT64_EQUAL(ncells_level1, CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X, test_name)

    
    -- all tests passed
    C.printf("Unit tests: \x1b[32mPASS\x1b[0m\n")
    C.free([&opaque](test_name))
  end
  return test_declarations
end

-- top level task

local top_level_task = make_test_declarations()
regentlib.start(top_level_task)

