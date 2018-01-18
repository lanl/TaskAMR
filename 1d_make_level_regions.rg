import "regent"
local C = regentlib.c

require("refinement_bits")

-- like stdlib's pow
function pow(x, y)
  local value = x
  for i = 2,y do
    value = value * x
  end
  return value
end

local function declare_bloated_partition(bloated_partition, old_partition, old_region, num_pts)
  local declaration = rquote
    var coloring = C.legion_domain_point_coloring_create()
    for color in old_partition.colors do
      var limits = old_partition[color].bounds
      var first : int64 = limits.lo - 1
      var last : int64 = limits.hi + 1
      if first < 0 then
        first = 0
      end
      if last >= num_pts then
        last = num_pts - 1
      end
      C.legion_domain_point_coloring_color_domain(coloring, [int1d](color), rect1d {first, last})
    end
    var [bloated_partition] = partition(aliased, [old_region], coloring, old_partition.colors)
  end
  return declaration
end

-- meta programming to create regions and partitions for levels 1 to MAX_REFINEMENT_LEVEL
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
    rquote var [meta_partition] = partition(equal, [meta_region], ispace(int1d, num_partitions))
      end

  local bloated_meta_partition = regentlib.newsymbol("level_" .. n .. "_bloated_meta_partition")
  local bmeta_declaration = declare_bloated_partition(bloated_meta_partition, meta_partition,
    meta_region, num_metas)

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
  local bpart_declaration = declare_bloated_partition(bloated_partition, cell_partition,
    cell_region, num_cells)

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
    meta_region, meta_declaration, meta_partition, mpart_declaration, bloated_meta_partition,
    bmeta_declaration

end

function declare_level_regions(meta_region_for_level,
                               cell_region_for_level,
                               face_region_for_level,
                               meta_partition_for_level,
                               cell_partition_for_level,
                               face_partition_for_level,
                               bloated_partition_for_level,
                               bloated_meta_partition_for_level,
                               MAX_REFINEMENT_LEVEL,
                               NUM_PARTITIONS)
  -- array of region and partition declarations
  local declarations = terralib.newlist()

  for n = 1, MAX_REFINEMENT_LEVEL do
    local cell_region, declare_cells, face_region, declare_faces, cell_partition, declare_cpart,
      face_partition, declare_fpart, bloated_partition, declare_bpart, meta_region,
      declare_meta, meta_partition, declare_mpart, bmeta_partition, declare_bmpart
      = make_level_regions(n, NUM_PARTITIONS)
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
    bloated_meta_partition_for_level:insert(bmeta_partition)
    declarations:insert(declare_bmpart)
  end

  return declarations
end

function init_num_cells_levels(num_cells,
                               MAX_REFINEMENT_LEVEL,
                               cell_region_for_level)

  local init_num_cells = terralib.newlist()
  init_num_cells:insert(rquote var [num_cells] end)
  for n = 1, MAX_REFINEMENT_LEVEL do
    init_num_cells:insert(rquote
      var limits = [cell_region_for_level[n]].ispace.bounds
      var  total : int64 = limits.hi - limits.lo + 1
      [num_cells][n] = total
    end)
  end

  return init_num_cells
end

