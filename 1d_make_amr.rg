-- Lua meta programming for AMR
import "regent"
local C = regentlib.c


-- meta programming to create partition by parent
function make_parent_partitions(n)

  local parent_cell_partition = regentlib.newsymbol("level_" .. n .. "_parent_cell_partition")

  local parent_meta_partition = regentlib.newsymbol("level_" .. n .. "_parent_meta_partition")

  local bloated_parent_meta_partition = regentlib.newsymbol("level_" .. n .. "_bloated_parent_meta_partition")

  local bloated_cell_partition_by_parent = regentlib.newsymbol("level_" .. n .. "_bloated_cell_partition_by_parent")

  return parent_cell_partition, parent_meta_partition, bloated_parent_meta_partition,
    bloated_cell_partition_by_parent
end -- make_parent_partitions


function insert_parent_partitions(parent_cell_partitions,
                                  parent_meta_partitions,
                                  bloated_parent_meta_partitions,
                                  bloated_cell_partitions_by_parent)

  for n = 1, MAX_REFINEMENT_LEVEL do
    local cell_partition, meta_partition, bloated_meta_partition, bloated_cell_partition
      = make_parent_partitions(n)
    parent_cell_partitions:insert(cell_partition)
    parent_meta_partitions:insert(meta_partition)
    bloated_parent_meta_partitions:insert(bloated_meta_partition)
    bloated_cell_partitions_by_parent:insert(bloated_cell_partition)
  end

end -- insert_parent_partitions


function initialize_parent_partitions(cell_partition_for_level,
                                      cell_region_for_level,
                                      parent_cell_partition_for_level,
                                      meta_partition_for_level,
                                      meta_region_for_level,
                                      parent_meta_partition_for_level,
                                      bloated_parent_meta_partition_for_level,
                                      bloated_cell_partition_by_parent_for_level,
                                      num_cells)

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

      var bloated_coloring = C.legion_domain_point_coloring_create()

      for color in [parent_meta_partition_for_level[n]].colors do
        var limits = [parent_meta_partition_for_level[n]][color].bounds
        var first_bloated : int64 = limits.lo - 1
        var last_bloated : int64 = limits.hi + 1
        if first_bloated < 0 then
          first_bloated = 0
        end
        if last_bloated >= num_cells[n] / CELLS_PER_BLOCK_X then
          last_bloated = (num_cells[n] / CELLS_PER_BLOCK_X) - 1
        end
        C.legion_domain_point_coloring_color_domain(bloated_coloring, [int1d](color), rect1d{first_bloated,
                                                                                             last_bloated})
      end -- color

      var [bloated_parent_meta_partition_for_level[n]] = partition(aliased,
                                                                   [meta_region_for_level[n]],
                                                                   bloated_coloring,
                                                                   [meta_partition_for_level[n-1]].colors)

      bloated_coloring = C.legion_domain_point_coloring_create()

      for color in [parent_cell_partition_for_level[n]].colors do
        var limits = [parent_cell_partition_for_level[n]][color].bounds
        var first_bloated : int64 = limits.lo - 1
        var last_bloated : int64 = limits.hi + 1
        if first_bloated < 0 then
          first_bloated = 0
        end
        if last_bloated >= num_cells[n] then
          last_bloated = num_cells[n] - 1
        end
        C.legion_domain_point_coloring_color_domain(bloated_coloring, [int1d](color),
                                                    rect1d{first_bloated, last_bloated})
      end -- color

      var [bloated_cell_partition_by_parent_for_level[n]] = partition(aliased,
                                                                      [cell_region_for_level[n]],
                                                                      bloated_coloring,
                                                             [cell_partition_for_level[n-1]].colors)


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
      for color in [meta_partition_for_level[n]].colors do
        writeAMRCells([num_cells][n], [meta_partition_for_level[n]][color],
                      [cell_partition_for_level[n]][color])
      end
    end)
  end

  return write_cells
end -- make_write_cells


function make_print_grid(meta_partition_for_level,
                         cell_partition_for_level)

  local print_grid = terralib.newlist()

  for n = 1, MAX_REFINEMENT_LEVEL do
    print_grid:insert(rquote
      __demand(__parallel)
      for color in [meta_partition_for_level[n]].colors do
        printAMRCells([n], [meta_partition_for_level[n]][color],
                      [cell_partition_for_level[n]][color])
      end
    end)
  end

  return print_grid
end -- make_print_grid


function make_init_regrid_and_values(num_cells,
                                     dx,
                                     cell_partition_for_level,
                                     bloated_partition_for_level,
                                     face_partition_for_level,
                                     meta_partition_for_level)

  local init_regrid_and_values = terralib.newlist()

  for level = 1, MAX_REFINEMENT_LEVEL do

    init_regrid_and_values:insert(rquote

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        initializeCells([num_cells][level], [cell_partition_for_level[level]][color])
      end

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        calculateGradient(num_cells[level], [dx][level], [bloated_partition_for_level[level]][color],
                          [face_partition_for_level[level]][color])
      end

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        flagRegrid([meta_partition_for_level[level]][color],
                                   [face_partition_for_level[level]][color])
      end

    end)

  end -- level

  return init_regrid_and_values
end -- make_init_regrid_and_values


function make_init_grid_refinement(num_cells,
                                   cell_partition_for_level,
                                   meta_partition_for_level,
                                   bloated_meta_partition_for_level,
                                   parent_meta_partition_for_level,
                                   parent_cell_partition_for_level,
                                   bloated_parent_meta_partition_for_level,
                                   meta_region_for_level)

  local init_grid_refinement = terralib.newlist()

  for max_level = 1, MAX_REFINEMENT_LEVEL - 1 do

    for level = 1, max_level do

      init_grid_refinement:insert(rquote

        __demand(__parallel)
        for color in [cell_partition_for_level[level]].colors do
          updateRefinement(num_cells[level]/CELLS_PER_BLOCK_X,
                           [meta_partition_for_level[level]][color],
                           [cell_partition_for_level[level]][color],
                           [bloated_meta_partition_for_level[level]][color],
                           [parent_meta_partition_for_level[level+1]][color],
                           [parent_cell_partition_for_level[level+1]][color],
                           [bloated_parent_meta_partition_for_level[level+1]][color])
        end

        fill([meta_region_for_level[level]].needsRefinement, false)

        __demand(__parallel)
        for color in [cell_partition_for_level[level]].colors do
          smoothGrid([meta_partition_for_level[level]][color])
        end

      end)

    end -- level

  end -- max_level

  return init_grid_refinement
end -- make_init_grid_refinement


function make_time_step(num_cells,
                        dx,
                        cell_region_for_level,
                        face_partition_for_level,
                        cell_partition_for_level,
                        meta_partition_for_level,
                        bloated_partition_for_level,
                        bloated_cell_partition_by_parent_for_level,
                        parent_cell_partition_for_level)

  local time_step = terralib.newlist()

  -- interpolateToChildren works from phi_copy not phi or __demand(__parallel) fails
  time_step:insert(rquote
    copy([cell_region_for_level[MAX_REFINEMENT_LEVEL]].phi,
         [cell_region_for_level[MAX_REFINEMENT_LEVEL]].phi_copy)
  end)

  for level = 1, MAX_REFINEMENT_LEVEL - 1 do

    time_step:insert(rquote

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        copyToChildren([meta_partition_for_level[level]][color],
                       [cell_partition_for_level[level]][color],
                       [parent_cell_partition_for_level[level+1]][color])
      end

    end)

  end -- level

  for level = 1, MAX_REFINEMENT_LEVEL - 1 do

    time_step:insert(rquote

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        interpolateToChildren(num_cells[level+1],
                              [meta_partition_for_level[level]][color],
                              [bloated_partition_for_level[level]][color],
                              [bloated_cell_partition_by_parent_for_level[level+1]][color],
                              [parent_cell_partition_for_level[level+1]][color])
      end

    end)

  end -- level

  for level = 1, MAX_REFINEMENT_LEVEL - 1 do

    time_step:insert(rquote

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        calculateAMRFlux(num_cells[level],
                         dx[level],
                         DT,
                         [meta_partition_for_level[level]][color],
                         [bloated_partition_for_level[level]][color],
                         [bloated_cell_partition_by_parent_for_level[level+1]][color],
                         [face_partition_for_level[level]][color])
      end

    end)

  end -- level
  
  time_step:insert(rquote

    __demand(__parallel)
    for color in [cell_partition_for_level[MAX_REFINEMENT_LEVEL]].colors do
      calculateFlux(num_cells[MAX_REFINEMENT_LEVEL],
                    dx[MAX_REFINEMENT_LEVEL],
                    DT,
                    [meta_partition_for_level[MAX_REFINEMENT_LEVEL]][color],
                    [bloated_partition_for_level[MAX_REFINEMENT_LEVEL]][color],
                    [face_partition_for_level[MAX_REFINEMENT_LEVEL]][color])
    end

  end)

  for level = 1, MAX_REFINEMENT_LEVEL do

    time_step:insert(rquote

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        applyFlux(dx[level],
                  DT,
                  [meta_partition_for_level[level]][color],
                  [cell_partition_for_level[level]][color],
                  [face_partition_for_level[level]][color])
      end

    end)

  end -- level

  return time_step
end -- make_time_step


function make_flag_regrid(num_cells,
                     dx,
                     needs_regrid,
                     do_regrid,
                     face_partition_for_level,
                     meta_partition_for_level,
                     bloated_partition_for_level,
                     bloated_cell_partition_by_parent_for_level)

  local flag_regrid = terralib.newlist()


  for level = 1, MAX_REFINEMENT_LEVEL - 1 do

    flag_regrid:insert(rquote

      __demand(__parallel)
      for color in [meta_partition_for_level[level]].colors do
        calculateAMRGradient(num_cells[level],
                         dx[level],
                         [meta_partition_for_level[level]][color],
                         [bloated_partition_for_level[level]][color],
                         [bloated_cell_partition_by_parent_for_level[level+1]][color],
                         [face_partition_for_level[level]][color])
      end

    end)

  end -- level
  
  flag_regrid:insert(rquote

    __demand(__parallel)
    for color in [meta_partition_for_level[MAX_REFINEMENT_LEVEL]].colors do
      calculateGradient(num_cells[MAX_REFINEMENT_LEVEL],
                        dx[MAX_REFINEMENT_LEVEL],
                        [bloated_partition_for_level[MAX_REFINEMENT_LEVEL]][color],
                        [face_partition_for_level[MAX_REFINEMENT_LEVEL]][color])
    end

  end)

  for level = 1, MAX_REFINEMENT_LEVEL do

    flag_regrid:insert(rquote

      needs_regrid[level] = 0
      __demand(__parallel)
      for color in [meta_partition_for_level[level]].colors do
        needs_regrid[level] += flagRegrid([meta_partition_for_level[level]][color],
                                          [face_partition_for_level[level]][color])
      end

    end)

  end -- level
  
  flag_regrid:insert(rquote
    var [do_regrid] = 0
   end)

  for level = 1, MAX_REFINEMENT_LEVEL do

    flag_regrid:insert(rquote

      [do_regrid] += needs_regrid[level]

    end)

  end -- level

  return flag_regrid
end -- make_flag_regrid


function make_do_regrid(num_cells,
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

  local do_regrid = terralib.newlist()

  -- interpolateToChildren works from phi_copy not phi or __demand(__parallel) fails
  do_regrid:insert(rquote
    copy([cell_region_for_level[MAX_REFINEMENT_LEVEL]].phi,
         [cell_region_for_level[MAX_REFINEMENT_LEVEL]].phi_copy)
  end)

  for level = 1, MAX_REFINEMENT_LEVEL - 1 do

    do_regrid:insert(rquote

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        copyToChildren([meta_partition_for_level[level]][color],
                       [cell_partition_for_level[level]][color],
                       [parent_cell_partition_for_level[level+1]][color])
      end

    end)

  end -- level

  for level = 1, MAX_REFINEMENT_LEVEL - 1 do

    do_regrid:insert(rquote

      __demand(__parallel)
      for color in [cell_partition_for_level[level]].colors do
        interpolateToChildren(num_cells[level+1],
                              [meta_partition_for_level[level]][color],
                              [bloated_partition_for_level[level]][color],
                              [bloated_cell_partition_by_parent_for_level[level+1]][color],
                              [parent_cell_partition_for_level[level+1]][color])
      end

    end)

  end -- level

  for max_level = 1, MAX_REFINEMENT_LEVEL - 1 do

    for level = 1, max_level do

      do_regrid:insert(rquote

        __demand(__parallel)
        for color in [cell_partition_for_level[level]].colors do
          updateRefinement(num_cells[level]/CELLS_PER_BLOCK_X,
                           [meta_partition_for_level[level]][color],
                           [cell_partition_for_level[level]][color],
                           [bloated_meta_partition_for_level[level]][color],
                           [parent_meta_partition_for_level[level+1]][color],
                           [parent_cell_partition_for_level[level+1]][color],
                           [bloated_parent_meta_partition_for_level[level+1]][color])
        end

        fill([meta_region_for_level[level]].needsRefinement, false)

        __demand(__parallel)
        for color in [cell_partition_for_level[level]].colors do
          smoothGrid([meta_partition_for_level[level]][color])
        end

      end)

    end -- level

  end -- max_level
  
  return do_regrid
end -- make_do_regrid


