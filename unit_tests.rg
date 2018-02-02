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


-- implement all required model APIs and link model.rg to your file
require("1d_make_levels")

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

  -- top_level task using previous meta programming
  local task test_declarations()
    [declarations];

    var test_name : &int8
    test_name = [&int8](C.malloc(81))

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

