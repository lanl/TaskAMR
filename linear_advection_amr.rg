-- linear advection with Lax-Friedrichs
import "regent"
local C = regentlib.c
local MATH = terralib.includec("math.h")

require("global_const")
require("refinement_bits")
require("linear_constants")

-- model specific tasks must implement following API:

task calculateGradient(num_cells : int64,
                   dx : double,
                   cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
where
  reads(cells.phi),
  writes(faces.grad)
do
  var vel : double = U
  var left_boundary_face : int64 = faces.ispace.bounds.lo
  var right_boundary_face : int64 = faces.ispace.bounds.hi

  var left_boundary_cell : int64 = cells.ispace.bounds.lo
  var right_boundary_cell : int64 = cells.ispace.bounds.hi

  var start_face : int64 = left_boundary_face
  var stop_face : int64 = right_boundary_face + 1

  if left_boundary_cell == 0 then
    start_face  = start_face + 1
  end
  if right_boundary_cell == (num_cells - 1) then
    stop_face  = stop_face - 1
  end

  -- loop on inner faces
  var cell_index : int64 = left_boundary_cell
  for face = start_face, stop_face do
    var left : double = cells[cell_index].phi
    cell_index = cell_index + 1
    var right : double = cells[cell_index].phi
    var grad : double = (right - left) / dx
    faces[face].grad = grad
  end

  -- boundary conditions: hold end cells constant in time
  if left_boundary_cell == 0 then
    faces[0].grad = 0.0
  end
  if right_boundary_cell == (num_cells - 1) then
    faces[right_boundary_face].grad = 0.0
  end
end -- calculateGradient


task flagRegrid(blocks: region(ispace(int1d), RefinementBits),
                faces: region(ispace(int1d), FaceValues))
                
where
  reads(faces.grad),
  reads writes(blocks.needsRefinement)
do
  var needs_regrid : int64 = 0
  var first_face : int64 = faces.ispace.bounds.lo

  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1

  for block = start_block, stop_block do
    var start_face : int64 = first_face + CELLS_PER_BLOCK_X * (block - start_block)
    var stop_face : int64 = start_face + CELLS_PER_BLOCK_X + 1
    for face = start_face, stop_face do
      if MATH.fabs(faces[face].grad) > MAX_GRAD then
        blocks[block].needsRefinement = true
        needs_regrid = 1
        C.printf("block %d <= %d < %d REFINE\n", start_block, block, stop_block)
      end
    end -- for face
  end -- for block

  return needs_regrid
end -- flagRegrid


__demand(__inline)
task parent_to_left(child : int64)
  return (child - 1) / 2
end


__demand(__inline)
task parent_to_right(child : int64)
  return (child + 1) / 2
end


__demand(__inline)
task linear_interpolate(x : double,
                        x0 : double,
                        x1 : double,
                        y0 : double,
                        y1 : double)
  return y0 + (x - x0) * (y1 - y0) / (x1 - x0)
end


task interpolateToChildren(num_children: int64,
                           blocks: region(ispace(int1d), RefinementBits),
                           parents: region(ispace(int1d), CellValues),
                           children: region(ispace(int1d), CellValues))
where
  reads (blocks.needsRefinement),
  reads (parents.phi),
  reads writes (children.phi)
do
  var first_child : int64 = children.ispace.bounds.lo
  var first_block : int64 = blocks.ispace.bounds.lo

  for block in blocks do
    if blocks[block].needsRefinement then
      var start_child : int64 = first_child +  2 * CELLS_PER_BLOCK_X * (block - first_block)
      var stop_child : int64 = start_child + 2 * CELLS_PER_BLOCK_X

      if start_child == 0 then
        children[0].phi = parents[0].phi
        start_child += 1
      end

      if stop_child == num_children then
        stop_child -= 1
        children[stop_child].phi = parents[parent_to_left(stop_child)].phi
      end

      for child = start_child, stop_child do
        var left_parent : int64 = parent_to_left(child)
        var right_parent : int64 = parent_to_right(child)
        var left_parent_x : double = 2.0 * [double](left_parent) + 0.5
        var right_parent_x : double = 2.0 * [double](right_parent) + 0.5
        children[child].phi = linear_interpolate([double](child), left_parent_x, right_parent_x,
                                                 parents[left_parent].phi, parents[right_parent].phi)
      end

    end -- needsRefinement
  end -- block

end -- interpolateToChildren


__demand(__inline)
task left_child(block : int64)
  return 2 * block
end


__demand(__inline)
task right_child(block : int64)
  return 2 * block + 1
end


__demand(__inline)
task left(block : int64)
  return block - 1
end


__demand(__inline)
task right(block : int64)
  return block + 1
end


task smoothGrid(blocks: region(ispace(int1d), RefinementBits))
where
  reads (blocks.cascadeRefinement),
  writes (blocks.{needsRefinement,
                  cascadeRefinement})
do
  for block in blocks do
    if blocks[block].cascadeRefinement then
      blocks[block].needsRefinement = true
      blocks[block].cascadeRefinement = false
      C.printf("block %d needsRefine\n", block)
    end
  end
end -- smoothGrid


task updateRefinementBits(num_blocks: int64,
                          blocks: region(ispace(int1d), RefinementBits),
                          ghosts: region(ispace(int1d), RefinementBits),
                          children: region(ispace(int1d), RefinementBits),
                          ghost_children: region(ispace(int1d), RefinementBits))
where
  reads (ghosts.needsRefinement),
  reads (ghost_children.{needsRefinement,
                         isRefined}),
  reads writes (blocks.{isActive,
                  isRefined,
                  cascadeRefinement,
                  plusXMoreRefined,
                  minusXMoreRefined,
                  plusXMoreCoarse,
                  minusXMoreCoarse}),
  writes (children.isActive)
do
  var needs_regrid : int64 = 0

  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1

  for block = start_block, stop_block do
    C.printf("block: %d <= %d < %d\n", start_block, block, stop_block)
    var left_refinement_delta : int64

    if block == 0 then
      blocks[block].minusXMoreCoarse = false
      blocks[block].minusXMoreRefined = false
    else
      if ghost_children[right_child(left(block))].needsRefinement then
        left_refinement_delta = 2
      elseif ghost_children[right_child(left(block))].isRefined then
        left_refinement_delta = 2
      elseif blocks[block].minusXMoreRefined then
        left_refinement_delta = 1
      elseif ghosts[left(block)].needsRefinement then
        left_refinement_delta = 1
      elseif blocks[block].minusXMoreCoarse then
        left_refinement_delta = -1
      else
        left_refinement_delta = 0
      end
    end -- left boundary

    var right_refinement_delta : int64

    if block == (num_blocks -1 ) then
      blocks[block].plusXMoreCoarse = false
      blocks[block].plusXMoreRefined = false
    else
      if ghost_children[left_child(right(block))].needsRefinement then
        right_refinement_delta = 2
      elseif ghost_children[left_child(right(block))].isRefined then
        right_refinement_delta = 2
      elseif blocks[block].plusXMoreRefined then
        right_refinement_delta = 1
      elseif ghosts[right(block)].needsRefinement then
        right_refinement_delta = 1
      elseif blocks[block].plusXMoreCoarse then
        right_refinement_delta = -1
      else
        right_refinement_delta = 0
      end
    end -- right boundary

    var my_refinement_delta : int64

    if ghosts[block].needsRefinement then

      blocks[block].isRefined = true
      blocks[block].isActive = false
      children[left_child(block)].isActive = true
      children[right_child(block)].isActive = true
      my_refinement_delta = 1
      C.printf("refine\n");

    elseif blocks[block].isActive then

      blocks[block].isRefined = false
      my_refinement_delta = 0

    end -- needsRefinement or isActive

    if ghosts[block].needsRefinement or blocks[block].isActive then
      if not (block == 0) then
        if left_refinement_delta > my_refinement_delta then
          blocks[block].minusXMoreCoarse = false
          blocks[block].minusXMoreRefined = true
        elseif left_refinement_delta == my_refinement_delta then
          blocks[block].minusXMoreCoarse = false
          blocks[block].minusXMoreRefined = false
        else
          blocks[block].minusXMoreCoarse = true
          blocks[block].minusXMoreRefined = false
        end
        C.printf("left delta %d\n", left_refinement_delta - my_refinement_delta)
        if (left_refinement_delta - my_refinement_delta) > 1 then
          blocks[block].cascadeRefinement = true
          needs_regrid = 1
          C.printf("cascade\n")
        end
      end -- not block 0

      if not (block == (num_blocks-1)) then
        if right_refinement_delta > my_refinement_delta then
          blocks[block].plusXMoreCoarse = false
          blocks[block].plusXMoreRefined = true
        elseif right_refinement_delta == my_refinement_delta then
          blocks[block].plusXMoreCoarse = false
          blocks[block].plusXMoreRefined = false
        else
          blocks[block].plusXMoreCoarse = true
          blocks[block].plusXMoreRefined = false
        end
        C.printf("right delta %d\n", right_refinement_delta - my_refinement_delta)
        if (right_refinement_delta - my_refinement_delta) > 1 then
          blocks[block].cascadeRefinement = true
          needs_regrid = 1
          C.printf("cascade\n")
        end
      end -- not last block
    end -- needsRefinement or isActive

  end -- block

  return needs_regrid
end -- updateRefinementBits


task writeAMRCells(ncells : int64,
                   blocks: region(ispace(int1d), RefinementBits),
                   cells: region(ispace(int1d), CellValues))
where
  reads(cells.phi),
  reads(blocks.isActive)
do
  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1
  var buf : &int8
  buf = [&int8](C.malloc(60))
  C.sprintf(buf, "linear_amr.%d.%d.txt", ncells, start_block)
  var fp = C.fopen(buf,"w")
  for block = start_block, stop_block do
    if blocks[block].isActive then
      var start_cell : int64 = block * CELLS_PER_BLOCK_X
      var stop_cell : int64 = (block + 1) * CELLS_PER_BLOCK_X
      for cell = start_cell, stop_cell do
        C.fprintf(fp, "%f %f\n", LENGTH_X * (cell + 0.5) / [double](ncells), cells[cell].phi)
      end
    end -- is Active
  end -- block
  C.fclose(fp)
  C.free([&opaque](buf))
end -- writeAMRCells


task copyToChildren(blocks: region(ispace(int1d), RefinementBits),
                    cells: region(ispace(int1d), CellValues),
                    children: region(ispace(int1d), CellValues))
where
  reads(cells.phi),
  reads(blocks.isActive),
  writes(children.phi)
do
end -- copyToChildren
