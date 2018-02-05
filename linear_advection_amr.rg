-- linear advection with Lax-Friedrichs
import "regent"
local C = regentlib.c
local MATH = terralib.includec("math.h")

require("global_const")
require("refinement_bits")
require("linear_constants")

-- model specific tasks must implement following API:


__demand(__inline)
task parent_to_left(child : int64)
  return (child - 1) / 2
end


__demand(__inline)
task parent_to_right(child : int64)
  return (child + 1) / 2
end


__demand(__inline)
task parent_x(parent : int64)
  return 2.0 * [double](parent) + 0.5
end


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


__demand(__inline)
task linear_interpolate(x : double,
                        x0 : double,
                        x1 : double,
                        y0 : double,
                        y1 : double)
  return y0 + (x - x0) * (y1 - y0) / (x1 - x0)
end


__demand(__inline)
task start_ghost_cell(block : int64)
  var start_cell : int64 = block * CELLS_PER_BLOCK_X - 1
  if block == 0 then
    start_cell += 1 -- no ghost cell on left boundary
  end
  return start_cell
end


__demand(__inline)
task stop_ghost_cell(block : int64)
  return (block + 1) * CELLS_PER_BLOCK_X + 1
end


__demand(__inline)
task first_face(block : int64,
                blocks_ispace : ispace(int1d),
                faces_ispace : ispace(int1d))
  var start_block : int64 = blocks_ispace.bounds.lo
  var start_face : int64 = faces_ispace.bounds.lo
  return start_face + (block - start_block) * CELLS_PER_BLOCK_X
end


__demand(__inline)
task last_face(block : int64,
               blocks_ispace : ispace(int1d),
               faces_ispace : ispace(int1d))
  var start_block : int64 = blocks_ispace.bounds.lo
  var start_face : int64 = faces_ispace.bounds.lo
  return start_face + (block - start_block + 1) * CELLS_PER_BLOCK_X
end


task calculateGradient(num_cells : int64,
                   dx : double,
                   bloated_cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
where
  reads(bloated_cells.phi),
  writes(faces.grad)
do
  var vel : double = U
  var left_boundary_face : int64 = faces.ispace.bounds.lo
  var right_boundary_face : int64 = faces.ispace.bounds.hi

  var left_boundary_cell : int64 = bloated_cells.ispace.bounds.lo
  var right_boundary_cell : int64 = bloated_cells.ispace.bounds.hi

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
    var left : double = bloated_cells[cell_index].phi
    cell_index = cell_index + 1
    var right : double = bloated_cells[cell_index].phi
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
  writes(blocks.needsRefinement)
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
      end
    end -- for face
  end -- for block

  return needs_regrid
end -- flagRegrid


task interpolateToChildren(num_children: int64,
                           blocks: region(ispace(int1d), RefinementBits),
                           ghost_parents: region(ispace(int1d), CellValues),
                           ghost_children: region(ispace(int1d), CellValues),
                           children: region(ispace(int1d), CellValues))
where
  reads (blocks.{isActive,
                 minusXMoreRefined,
                 minusXMoreCoarse,
                 plusXMoreRefined,
                 plusXMoreCoarse}),
  reads (ghost_parents.phi_copy),
  reads (ghost_children.phi_copy),
  writes (children.phi)
do
  var first_child : int64 = children.ispace.bounds.lo
  var first_block : int64 = blocks.ispace.bounds.lo

  for block in blocks do
    if blocks[block].isActive then
      var start_child : int64 = first_child +  2 * CELLS_PER_BLOCK_X * (block - first_block)
      var stop_child : int64 = start_child + 2 * CELLS_PER_BLOCK_X

      if start_child == 0 then
        children[0].phi = ghost_parents[0].phi_copy
        start_child += 1
      end

      if stop_child == num_children then
        stop_child -= 1
        children[stop_child].phi = ghost_parents[parent_to_left(stop_child)].phi_copy
      end

      var left_phi : double
      var left_x : double
      var right_phi : double
      var right_x : double

      if (blocks[block].minusXMoreRefined or blocks[block].minusXMoreCoarse) then

        if blocks[block].minusXMoreRefined then
          left_x = [double](left(start_child)) 
        elseif blocks[block].minusXMoreCoarse then
          left_x = parent_x(parent_to_left(start_child))
        end

        left_phi = ghost_children[left(start_child)].phi_copy
        right_x = parent_x(parent_to_right(start_child))
        right_phi = ghost_parents[parent_to_right(start_child)].phi_copy

        children[start_child].phi = linear_interpolate([double](start_child), left_x, right_x,
                                                       left_phi, right_phi)
        start_child += 1

      end

      if (blocks[block].plusXMoreRefined or blocks[block].plusXMoreCoarse) then

        stop_child -= 1

        if blocks[block].plusXMoreRefined then
          right_x = [double](right(stop_child)) 
        elseif blocks[block].plusXMoreCoarse then
          right_x = parent_x(parent_to_right(stop_child))
        end

        right_phi = ghost_children[right(stop_child)].phi_copy
        left_x = parent_x(parent_to_left(stop_child))
        left_phi = ghost_parents[parent_to_left(stop_child)].phi_copy

        children[stop_child].phi = linear_interpolate([double](stop_child), left_x, right_x,
                                                      left_phi, right_phi)

      end

      for child = start_child, stop_child do
        left_x = parent_x(parent_to_left(child))
        left_phi = ghost_parents[parent_to_left(child)].phi_copy
        right_x = parent_x(parent_to_right(child))
        right_phi = ghost_parents[parent_to_right(child)].phi_copy
        children[child].phi = linear_interpolate([double](child), left_x, right_x,
                                                 left_phi, right_phi)
      end

    end -- isActive
  end -- block

end -- interpolateToChildren


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
    end
  end
end -- smoothGrid


task updateRefinement(num_blocks: int64,
                      blocks: region(ispace(int1d), RefinementBits),
                      cells: region(ispace(int1d), CellValues),
                      ghosts: region(ispace(int1d), RefinementBits),
                      children: region(ispace(int1d), RefinementBits),
                      child_cells: region(ispace(int1d), CellValues),
                      ghost_children: region(ispace(int1d), RefinementBits))
where
  reads (ghosts.needsRefinement,
         children.wantsCoarsening,
         child_cells.phi,
         ghost_children.{needsRefinement,
                         isRefined}),
  reads writes (blocks.{isActive,
                  isRefined,
                  cascadeRefinement,
                  plusXMoreRefined,
                  minusXMoreRefined,
                  plusXMoreCoarse,
                  minusXMoreCoarse}),
  writes (children.isActive,
          cells.phi)
do
  var needs_regrid : int64 = 0

  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1

  for block = start_block, stop_block do
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

    if blocks[block].isActive then
      if ghosts[block].needsRefinement then

        blocks[block].isRefined = true
        blocks[block].isActive = false
        children[left_child(block)].isActive = true
        children[right_child(block)].isActive = true
        my_refinement_delta = 1

      else

        blocks[block].isRefined = false
        my_refinement_delta = 0

      end -- needsRefinement
    else
      if children[left_child(block)].wantsCoarsening
         and children[right_child(block)].wantsCoarsening then

        cells[block].phi = 0.5 * (child_cells[left_child(block)].phi
                                   + child_cells[right_child(block)].phi)
        blocks[block].isRefined = false
        blocks[block].isActive = true
        children[left_child(block)].isActive = false
        children[right_child(block)].isActive = false
        my_refinement_delta = 0

      end -- children want coarsening
    end -- isActive

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
        if (left_refinement_delta - my_refinement_delta) > 1 then
          blocks[block].cascadeRefinement = true
          needs_regrid = 1
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
        if (right_refinement_delta - my_refinement_delta) > 1 then
          blocks[block].cascadeRefinement = true
          needs_regrid = 1
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
  writes(cells.phi_copy),
  writes(children.phi_copy)
do
  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1
  for block = start_block, stop_block do
    if blocks[block].isActive then
      var start_cell : int64 = block * CELLS_PER_BLOCK_X
      var stop_cell : int64 = (block + 1) * CELLS_PER_BLOCK_X
      for cell = start_cell, stop_cell do
        cells[cell].phi_copy = cells[cell].phi
        children[left_child(cell)].phi_copy = cells[cell].phi
        children[right_child(cell)].phi_copy = cells[cell].phi
      end
    end -- is Active
  end -- block
end -- copyToChildren


-- duplicates code from calculateAMRGrad
task calculateAMRFlux(num_cells : int64,
                   dx : double,
                   dt : double,
                   blocks: region(ispace(int1d), RefinementBits),
                   bloated_cells: region(ispace(int1d), CellValues),
                   bloated_children: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
where
  reads(bloated_cells.phi,
        bloated_children.phi,
        blocks.{isActive,
                minusXMoreRefined,
                plusXMoreRefined},
        faces.flux),
  writes(faces.flux)
do
  var vel : double = U

  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1

  for block = start_block, stop_block do
    if blocks[block].isActive then
      var start_cell : int64 = start_ghost_cell(block)
      var stop_cell : int64 = stop_ghost_cell(block)
      var start_face : int64 = first_face(block, blocks.ispace, faces.ispace)
      var stop_face : int64 = last_face(block, blocks.ispace, faces.ispace) + 1

      if start_cell == 0 then
        start_face += 1 -- handle left boundary special
      end
      if stop_cell > num_cells then
        stop_face -= 1 -- handle right boundary special
      end
 
      var cell_index : int64 = start_cell

      if blocks[block].minusXMoreRefined then
        var left : double = bloated_children[right_child(cell_index)].phi
        cell_index = cell_index + 1
        var right : double = bloated_children[left_child(cell_index)].phi
        var flux : double = 0.5 * vel * (left + right) +0.25 * dx * (left - right)/dt
        faces[start_face].flux = flux
        start_face += 1
      end
 
      if blocks[block].plusXMoreRefined then
        stop_face -= 1
        var left : double = bloated_children[right_child(stop_cell - 2)].phi
        var right : double = bloated_children[left_child(stop_cell - 1)].phi
        var flux : double = 0.5 * vel * (left + right) +0.25 * dx * (left - right)/dt
        faces[stop_face].flux = flux
      end
 
      for face = start_face, stop_face do
        var left : double = bloated_cells[cell_index].phi
        cell_index = cell_index + 1
        var right : double = bloated_cells[cell_index].phi
        var flux : double = 0.5 * vel * (left + right) +0.5 * dx * (left - right)/dt
        faces[face].flux = flux
      end -- face

      -- boundary conditions: hold end cells constant in time
      if start_cell == 0 then
        faces[0].flux = faces[1].flux
      end
      if stop_cell > num_cells then
        faces[stop_face].flux = faces[stop_face-1].flux
      end

    end -- isActive
  end -- block

end --calculateAMRFlux


-- duplicates code from calculateAMRFlux
task calculateAMRGradient(num_cells : int64,
                   dx : double,
                   blocks: region(ispace(int1d), RefinementBits),
                   bloated_cells: region(ispace(int1d), CellValues),
                   bloated_children: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
where
  reads(bloated_cells.phi,
        bloated_children.phi,
        blocks.{isActive,
                minusXMoreRefined,
                plusXMoreRefined},
        faces.grad),
  writes(faces.grad)
do
  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1

  for block = start_block, stop_block do
    if blocks[block].isActive then
      var start_cell : int64 = start_ghost_cell(block)
      var stop_cell : int64 = stop_ghost_cell(block)
      var start_face : int64 = first_face(block, blocks.ispace, faces.ispace)
      var stop_face : int64 = last_face(block, blocks.ispace, faces.ispace) + 1

      if start_cell == 0 then
        start_face += 1 -- handle left boundary special
      end
      if stop_cell > num_cells then
        stop_face -= 1 -- handle right boundary special
      end
 
      var cell_index : int64 = start_cell

      if blocks[block].minusXMoreRefined then
        var left : double = bloated_children[right_child(cell_index)].phi
        cell_index = cell_index + 1
        var right : double = bloated_children[left_child(cell_index)].phi
        var grad : double = 2.0 * (right - left) / dx
        faces[start_face].grad = grad
        start_face += 1
      end
 
      if blocks[block].plusXMoreRefined then
        stop_face -= 1
        var left : double = bloated_children[right_child(stop_cell - 2)].phi
        var right : double = bloated_children[left_child(stop_cell - 1)].phi
        var grad : double = 2.0 * (right - left) / dx
        faces[stop_face].grad = grad
      end
 
      for face = start_face, stop_face do
        var left : double = bloated_cells[cell_index].phi
        cell_index = cell_index + 1
        var right : double = bloated_cells[cell_index].phi
        var grad : double = (right - left) / dx
        faces[face].grad = grad
      end -- face

      -- boundary conditions: hold end cells constant in time
      if start_cell == 0 then
        faces[0].grad = 0.0
      end
      if stop_cell > num_cells then
        faces[stop_face].grad = 0.0
      end

    end -- isActive
  end -- block

end --calculateAMRGradient



