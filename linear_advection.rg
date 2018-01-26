-- linear advection with Lax-Friedrichs
import "regent"
local C = regentlib.c
local MATH = terralib.includec("math.h")

require("global_const")
require("refinement_bits")
require("linear_constants")

-- model specific tasks must implement following API:
-- task initializeCells(num_cells : int64,
--                      cell_region: region(ispace(int1d), CellValues))
-- task calculateFlux(num_cells : int64,
--                    dx : double,
--                    dt : double,
--                    cells: region(ispace(int1d), CellValues),
--                    faces: region(ispace(int1d), FaceValues))
-- task applyFlux(dx : double,
--                dt : double,
--                cells: region(ispace(int1d), CellValues),
--                faces: region(ispace(int1d), FaceValues))
-- task writeCells(num_cells : int64,
--                 cells: region(ispace(int1d), CellValues))


task initializeCells(num_cells : int64,
                     cell_region: region(ispace(int1d), CellValues))
where
  writes(cell_region.phi)
do
  C.printf("initializeCells %d cells\n", num_cells)
  for cell in cell_region.ispace do
    C.printf("%d:", cell)
    if [int64](cell) < (num_cells/2) then
      cell_region[cell].phi = 1.0
      C.printf("1 ")
    else
      cell_region[cell].phi = 0.0
      C.printf("0 ")
    end
  end
  C.printf("\n")
end -- initializeCells


task applyFlux(dx : double,
               dt : double,
               cells: region(ispace(int1d), CellValues),
               faces: region(ispace(int1d), FaceValues))
where
  reads(cells.phi,
        faces.flux),
  writes(cells.phi)
do
  var face_index : int64 = faces.ispace.bounds.lo
  for cell in cells do
    cells[cell].phi = cells[cell].phi - dt * (faces[face_index+1].flux - faces[face_index].flux) / dx
    face_index = face_index + 1
  end
end -- applyFlux


task calculateFlux(num_cells : int64,
                   dx : double,
                   dt : double,
                   cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
where
  reads(cells.phi,
        faces.flux),
  writes(faces.flux)
do
  C.printf("calculateFlux %d cells\n", num_cells)
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
    var flux : double = 0.5 * vel * (left + right) +0.5 * dx * (left - right)/dt
    faces[face].flux = flux
    C.printf("face %d, flux %f\n",face,flux)
  end

  -- boundary conditions: hold end cells constant in time
  if left_boundary_cell == 0 then
    faces[0].flux = faces[1].flux
  end
  if right_boundary_cell == (num_cells - 1) then
    faces[right_boundary_face].flux = faces[right_boundary_face-1].flux
  end
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

