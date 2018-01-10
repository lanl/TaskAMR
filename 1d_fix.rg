-- 1d fixed grid meta programming for all models
import "regent"
local C = regentlib.c

-- implement all required model APIs and link model.rg to your file
require("model")

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

  local cell_region = regentlib.newsymbol("level_" .. n .. "_cell_region")
  local ratio_to_level1 = pow(2, n) / 2
  local num_cells = CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X * ratio_to_level1
  local cell_declaration =
    rquote var [cell_region] = region(ispace(int1d, num_cells), CellValues) end

  local cell_partition = regentlib.newsymbol("level_" .. n .. "_cell_partition")
  local cpart_declaration =
    rquote var [cell_partition] = partition(equal, [cell_region], ispace(int1d, NUM_PARTITIONS)) end

  local bloated_partition = regentlib.newsymbol("level_" .. n .. "_bloated_partition")
  local bpart_declaration = rquote
    var coloring = C.legion_domain_point_coloring_create()
    for color in cell_partition.colors do
      var limits = cell_partition[color].bounds
      var first_cell : int64 = limits.lo - 1
      var last_cell : int64 = limits.hi + 1
      if first_cell < 0 then
        first_cell = 0
      end
      if last_cell >= num_cells then
        last_cell = num_cells - 1
      end
      C.legion_domain_point_coloring_color_domain(coloring, [int1d](color), rect1d {first_cell,
        last_cell})
    end
    var [bloated_partition] = partition(aliased, [cell_region], coloring, cell_partition.colors)
  end

  local face_region = regentlib.newsymbol("level_" .. n .. "_face_region")
  local ratio_to_level1 = pow(2, n) / 2
  local num_faces = num_cells + NUM_PARTITIONS
  local face_declaration =
    rquote var [face_region] = region(ispace(int1d, num_faces), FaceValues) end

  local face_partition = regentlib.newsymbol("level_" .. n .. "_face_partition")
  local fpart_declaration = rquote
    var coloring = C.legion_domain_point_coloring_create()
    var first_face : int64 = 0
    for color in cell_partition.colors do
      var limits = cell_partition[color].bounds
      var num_cells : int64 = limits.hi - limits.lo + 1
      var last_face : int64 = first_face + num_cells
      C.legion_domain_point_coloring_color_domain(coloring, [int1d](color), rect1d {first_face,
        last_face})
      first_face = first_face + num_cells + 1
    end
    var [face_partition] = partition(disjoint, [face_region], coloring, cell_partition.colors)
  end

  return cell_region, cell_declaration, face_region, face_declaration, cell_partition,
    cpart_declaration, face_partition, fpart_declaration, bloated_partition, bpart_declaration
end

-- meta programming to create top_level_task
function make_top_level_task()

  -- arrays of region by level
  local cell_region_for_level = terralib.newlist()
  local face_region_for_level = terralib.newlist()
  local cell_partition_for_level = terralib.newlist()
  local face_partition_for_level = terralib.newlist()
  local bloated_partition_for_level = terralib.newlist()

  -- array of region declarations
  local declarations = terralib.newlist()

  for n = 1, MAX_REFINEMENT_LEVEL do
    local cell_region, declare_cells, face_region, declare_faces, cell_partition, declare_cpart,
      face_partition, declare_fpart, bloated_partition, declare_bpart = make_level_regions(n)
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

  -- create top_level task with previous meta programming
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

