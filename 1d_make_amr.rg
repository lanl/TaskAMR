-- Lua meta programming for AMR
import "regent"
local C = regentlib.c


-- meta programming to create partition by parent
function make_parent_partitions(n)

  local parent_cell_partition = regentlib.newsymbol("level_" .. n .. "_parent_cell_partition")

  local parent_meta_partition = regentlib.newsymbol("level_" .. n .. "_parent_meta_partition")

  return parent_cell_partition, parent_meta_partition
end -- make_parent_partitions


function insert_parent_partitions(parent_cell_partitions, parent_meta_partitions)

  for n = 1, MAX_REFINEMENT_LEVEL do
    local cell_partition, meta_partition = make_parent_partitions(n)
    parent_cell_partitions:insert(cell_partition)
    parent_meta_partitions:insert(meta_partition)
  end

end -- insert_parent_partitions


function initialize_parent_partitions(cell_partition_for_level,
                                      cell_region_for_level,
                                      parent_cell_partition_for_level,
                                      meta_partition_for_level,
                                      meta_region_for_level,
                                      parent_meta_partition_for_level)

  local init_parent_partitions =  terralib.newlist()
  for n = 2, MAX_REFINEMENT_LEVEL do
    init_parent_partitions:insert(rquote
      var coloring = C.legion_domain_point_coloring_create()

      for color in [cell_partition_for_level[n-1]].colors do
        var limits = [cell_partition_for_level[n-1]][color].bounds
        var first_cell : int64 = 2 * limits.lo
        var last_cell : int64 = 2 * (limits.hi + 1) - 1
        C.legion_domain_point_coloring_color_domain(coloring, [int1d](color), rect1d{first_cell,
          last_cell})
      end -- color

      var [parent_cell_partition_for_level[n]] = partition(disjoint, [cell_region_for_level[n]],
        coloring, [cell_partition_for_level[n-1]].colors)

      var meta_coloring = C.legion_domain_point_coloring_create()

      for color in [meta_partition_for_level[n-1]].colors do
        var limits = [meta_partition_for_level[n-1]][color].bounds
        var first_meta : int64 = 2 * limits.lo
        var last_meta : int64 = 2 * (limits.hi + 1) - 1
        C.legion_domain_point_coloring_color_domain(meta_coloring, [int1d](color), rect1d{first_meta,
          last_meta})
      end -- color

      var [parent_meta_partition_for_level[n]] = partition(disjoint, [meta_region_for_level[n]], meta_coloring,
        [meta_partition_for_level[n-1]].colors)

    end)
  end -- n

  return init_parent_partitions
end -- initialize_parent_partitions


function make_init_activity(meta_region_for_level)

  local init_activity = terralib.newlist()

  init_activity:insert(rquote fill([meta_region_for_level[1]].isActive, true) end)
  for n = 2, MAX_REFINEMENT_LEVEL do
    init_activity:insert(rquote
      fill([meta_region_for_level[n]].isActive, false)
    end)
  end

  return init_activity
end -- make_init_activity


function make_write_cells(num_cells,
                          meta_partition_for_level,
                          cell_partition_for_level)

  local write_cells = terralib.newlist()

  for n = 1, MAX_REFINEMENT_LEVEL do
    write_cells:insert(rquote
      __demand(__parallel)
      for color in [meta_partition_for_level[1]].colors do
        writeAMRCells([num_cells][n], [meta_partition_for_level[n]][color],
                      [cell_partition_for_level[n]][color])
      end
    end)
  end

  return write_cells
end -- make_write_cells


function make_init_grid_and_values(num_cells,
                                   dx,
                                   cell_region_for_level,
                                   cell_partition_for_level,
                                   bloated_partition_for_level,
                                   face_partition_for_level,
                                   meta_partition_for_level,
                                   parent_cell_partition_for_level,
                                   bloated_meta_partition_for_level,
                                   parent_meta_partition_for_level,
                                   meta_region_for_level)

  local init_grid_and_values = terralib.newlist()

  local level = 1

  init_grid_and_values:insert(rquote

    initializeCells([num_cells][level], [cell_region_for_level[level]])

    __demand(__parallel)
    for color in [cell_partition_for_level[level]].colors do
      calculateGradient(num_cells[level], [dx][level], [bloated_partition_for_level[level]][color],
                        [face_partition_for_level[level]][color])
    end

    __demand(__parallel)
    for color in [cell_partition_for_level[level]].colors do
      printFaces([face_partition_for_level[level]][color])
    end

    var needs_regrid = 0
    __demand(__parallel)
    for color in [cell_partition_for_level[level]].colors do
      needs_regrid += flagRegrid([meta_partition_for_level[level]][color],
                                 [face_partition_for_level[level]][color])
    end
    C.printf("Needs regrid = %d\n", needs_regrid)

    if needs_regrid > 0 then

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        interpolateToChildren(num_cells[level+1], [meta_partition_for_level[level]][color],
                              [bloated_partition_for_level[level]][color],
                              [parent_cell_partition_for_level[level+1]][color])
      end

      __demand(__parallel)
      for color in [meta_partition_for_level[level]].colors do
        updateRefinementBits(num_cells[level]/CELLS_PER_BLOCK_X, [meta_partition_for_level[level]][color],
                             [bloated_meta_partition_for_level[level]][color],
                             [parent_meta_partition_for_level[level+1]][color])
      end
      fill([meta_region_for_level[level]].needsRefinement, false)

    end -- needs_regrid

  end)

  return init_grid_and_values
end -- make_init_grid_and_values

