import "regent"
local C = regentlib.c

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
function make_level_regions(n, num_partitions)

  local ratio_to_level1 = pow(2, n) / 2

  local cell_region = regentlib.newsymbol("level_" .. n .. "_cell_region")
  local num_cells = CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X * ratio_to_level1
  local cell_declaration =
    rquote var [cell_region] = region(ispace(int1d, num_cells), CellValues) end

  -- meta data is refinement status
  local meta_region = regentlib.newsymbol("level_" .. n .. "_meta_region")
  local num_metas = LEVEL_1_BLOCKS_X * ratio_to_level1
  local meta_declaration =
    rquote var [meta_region] = region(ispace(int1d, num_metas), RefinementBits) end

  local meta_partition = regentlib.newsymbol("level_" .. n .. "_meta_partition")
  local mpart_declaration =
    rquote var [meta_partition] = partition(equal, [meta_region], ispace(int1d, num_partitions)) end

  local cell_partition = regentlib.newsymbol("level_" .. n .. "_cell_partition")
  local cpart_declaration = rquote
    var coloring = C.legion_domain_point_coloring_create()
    for color in meta_partition.colors do
      var limits = meta_partition[color].bounds
      var first_cell : int64 = limits.lo * CELLS_PER_BLOCK_X
      var last_cell : int64 = (limits.hi + 1) * CELLS_PER_BLOCK_X - 1
      C.legion_domain_point_coloring_color_domain(coloring, [int1d](color), rect1d {first_cell,
        last_cell})
    end
    var [cell_partition] = partition(disjoint, [cell_region], coloring, meta_partition.colors)
  end

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
  local num_faces = num_cells + num_partitions
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
    cpart_declaration, face_partition, fpart_declaration, bloated_partition, bpart_declaration,
    meta_region, meta_declaration, meta_partition, mpart_declaration
end


