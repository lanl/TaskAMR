--Copyright (c) 2018, Triad National Security, LLC
--All rights reserved.

--This program was produced under U.S. Government contract 89233218CNA000001 for
--Los Alamos National Laboratory (LANL), which is operated by Triad National
--Security, LLC for the U.S. Department of Energy/National Nuclear Security
--Administration.

--THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS
--IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
--IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
--DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE
--LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
--CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
--SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
--INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
--CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
--ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
--POSSIBILITY OF SUCH DAMAGE.

--If software is modified to produce derivative works, such modified software should be
--clearly marked, so as not to confuse it with the version available from LANL.
-- linear advection with Lax-Friedrichs
import "regent"
local C = regentlib.c
local MATH = terralib.includec("math.h")

require("global_const")
require("refinement_bits")
require("linear_constants")


task initializeCells(num_cells : int64,
                     cell_region: region(ispace(int1d), CellValues))
where
  writes(cell_region.phi)
do
  for cell in cell_region.ispace do
    if [int64](cell) < (num_cells/2) then
      cell_region[cell].phi = 1.0
    else
      cell_region[cell].phi = 0.0
    end
  end
end -- initializeCells



task applyFlux(dx : double,
               dt : double,
               blocks: region(ispace(int1d), RefinementBits),
               cells: region(ispace(int1d), CellValues),
               faces: region(ispace(int1d), FaceValues))
where
  reads(blocks.isActive,
        cells.phi,
        faces.flux),
  writes(cells.phi)
do

  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1
  var first_face : int64 = faces.ispace.bounds.lo

  for block = start_block, stop_block do
    if blocks[block].isActive then
      var start_cell : int64 = block * CELLS_PER_BLOCK_X  -- should be modularized
      var stop_cell : int64 = (block + 1) * CELLS_PER_BLOCK_X
      var face_index : int64 = first_face + (block - start_block) * CELLS_PER_BLOCK_X

      for cell = start_cell, stop_cell do
        cells[cell].phi = cells[cell].phi - dt * (faces[face_index+1].flux
                                                  - faces[face_index].flux) / dx
        face_index = face_index + 1
      end -- cell

    end -- isActive
  end -- block
end -- applyFlux


task calculateFlux(num_cells : int64,
                   dx : double,
                   dt : double,
                   blocks: region(ispace(int1d), RefinementBits),
                   bloated_cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
where
  reads(bloated_cells.phi,
        blocks.isActive,
        faces.flux),
  writes(faces.flux)
do
  var vel : double = U

  var start_block : int64 = blocks.ispace.bounds.lo
  var stop_block : int64 = blocks.ispace.bounds.hi + 1
  var first_face : int64 = faces.ispace.bounds.lo

  for block = start_block, stop_block do
    if blocks[block].isActive then
      var start_cell : int64 = block * CELLS_PER_BLOCK_X - 1  -- should be modularized
      var stop_cell : int64 = (block + 1) * CELLS_PER_BLOCK_X + 1
      var start_face : int64 = first_face + (block - start_block) * CELLS_PER_BLOCK_X
      var stop_face : int64 = start_face + CELLS_PER_BLOCK_X + 1

      if start_cell == - 1 then
        start_face += 1
        start_cell += 1
      end
      if stop_cell > num_cells then
        stop_face -= 1
      end
  
      var cell_index : int64 = start_cell
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

end --calculateFlux


task writeCells(nx : int64,
                cells: region(ispace(int1d), CellValues))
where
  reads(cells.phi)
do
  var first_cell : int64 = cells.ispace.bounds.lo
  var last_cell : int64 = cells.ispace.bounds.hi
  var buf : &int8
  buf = [&int8](C.malloc(40))
  C.sprintf(buf, "linear.%d.txt", nx)
  var fp = C.fopen(buf,"w")
  for cell in cells do
    C.fprintf(fp, "%f\n", cells[cell].phi)
  end
  C.fclose(fp)
  C.free([&opaque](buf))
end -- writeCells


task printCells(cell_region: region(ispace(int1d), CellValues))
where
  reads(cell_region.phi)
do
  C.printf("phi: ")
  for cell in cell_region do
    C.printf("%f ", cell_region[cell].phi)
  end
  C.printf("\n")
end -- printCells


task printFaces(faces: region(ispace(int1d), FaceValues))
where
  reads(faces.grad)
do
  C.printf("grad: ")
  for face in faces do
    C.printf("%d:%f ", face, faces[face].grad)
  end
  C.printf("\n")
end -- printFaces

