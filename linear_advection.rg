import "regent"
local C = regentlib.c

-- required AMR global constants
CELLS_PER_BLOCK_X = 2 -- must be multilpe of 2
LEVEL_1_BLOCKS_X = 5
MAX_REFINEMENT_LEVEL = 4
NUM_PARTITIONS = 7

-- constants that should not exist
local NX = 80
local MAX_NX = 640

-- model specific global constants
local CFL = 0.5
T_FINAL = 0.25
local U = 1.0
DX = 1.0 / NX
local MIN_DX = 1.0 / MAX_NX
DT = CFL * MIN_DX / U

-- model specific fields must be in fspace's CellValues and FaceValues

fspace CellValues
{
  phi : double
}

fspace FaceValues
{
  flux : double
}

-- model specific tasks must implement following API:
-- task initializeCells(cell_region: region(ispace(int1d), CellValues))
-- task calculateFlux(dx : double,
--                    dt : double,
--                    cells: region(ispace(int1d), CellValues),
--                    faces: region(ispace(int1d), FaceValues))
-- task applyFlux(dx : double,
--                dt : double,
--                cells: region(ispace(int1d), CellValues),
--                faces: region(ispace(int1d), FaceValues))
-- task writeCells(cells: region(ispace(int1d), CellValues))

task initializeCells(cell_region: region(ispace(int1d), CellValues))
where
  writes(cell_region.phi)
do
  var size : int64 = cell_region.ispace.bounds.hi - cell_region.ispace.bounds.lo + 1
  for cell in cell_region.ispace do
    if [int64](cell) < (size/2) then
      cell_region[cell].phi = 1.0
    else
      cell_region[cell].phi = 0.0
    end
  end
  C.printf("initializeCells %d cells\n", size)
end

-- hard coded linear advection with Lax-Friedrichs for now, metaprog soon
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
end

-- hard coded linear advection with Lax-Friedrichs for now, metaprog soon
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
  end

  -- boundary conditions: hold end cells constant in time
  if left_boundary_cell == 0 then
    faces[0].flux = faces[1].flux
  end
  if right_boundary_cell == (num_cells - 1) then
    faces[right_boundary_face].flux = faces[right_boundary_face-1].flux
  end
end

task writeCells(cells: region(ispace(int1d), CellValues))
where
  reads(cells.phi)
do
  var first_cell : int64 = cells.ispace.bounds.lo
  var last_cell : int64 = cells.ispace.bounds.hi
  --var nx : int64 = last_cell - first_cell + 1
  var fp = C.fopen(["output." .. NX .. ".txt"],"w")
  for cell in cells do
    C.fprintf(fp, "%f\n", cells[cell].phi)
  end
  C.fclose(fp)
end


task printCells(cell_region: region(ispace(int1d), CellValues))
where
  reads(cell_region.phi)
do
  C.printf("phi: ")
  for cell in cell_region do
    C.printf("%f ", cell_region[cell].phi)
  end
  C.printf("\n")
end

