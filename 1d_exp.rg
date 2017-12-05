import "regent"
local C = regentlib.c

-- global constants

local CELLS_PER_BLOCK_X = 4
local LEVEL_1_BLOCKS_X = 3
local MAX_REFINEMENT_LEVEL = 3
local NX = 48
local CFL = 0.5
local U = 1.0
local DX = 1.0 / (NX - 1)
local DT = CFL * DX / U

fspace CellValues
{
  phi : double
}

fspace FaceValues
{
  flux : double
}


-- like stdlib's
function pow(x, y)
  local value = x
  for i = 2,y do
    value = value * x
  end
  return value
end

-- for creating regions for levels 1 to MAX_REFINEMENT_LEVEL
function make_level_regions(n)
  local level_cell = regentlib.newsymbol("level_" .. n .. "_cell")
  local ratio_to_level1 = pow(2, n) / 2
  local num_cells = CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X * ratio_to_level1
  local cells = rquote var [level_cell] = region(ispace(int1d, num_cells), CellValues) end

  local level_face = regentlib.newsymbol("level_" .. n .. "_face")
  local ratio_to_level1 = pow(2, n) / 2
  local num_faces = num_cells + 1
  local faces = rquote var [level_face] = region(ispace(int1d, num_faces), FaceValues) end
  return level_cell, cells, level_face, faces
end

function make_top_level_task()
  local level_cells = terralib.newlist()
  local level_faces = terralib.newlist()
  local declare_level_regions = terralib.newlist()
  for n = 1, MAX_REFINEMENT_LEVEL do
    local level_cell, declare_cells, level_face, declare_faces = make_level_regions(n)
    level_cells:insert(level_cell)
    declare_level_regions:insert(declare_cells)
    level_faces:insert(level_face)
    declare_level_regions:insert(declare_faces)
  end
  local loops = terralib.newlist()
  for n = 1, MAX_REFINEMENT_LEVEL do
    loops:insert(rquote
      var total : int64 = 0
      for cell in [level_cells[n]] do
        total = total + 1
      end
      C.printf(["level " .. n .. " cells %d \n"], total)
    end)
  end
  local task top_level()
    [declare_level_regions];
    [loops];
    initializeCells([level_cells[MAX_REFINEMENT_LEVEL]])
    printCells([level_cells[MAX_REFINEMENT_LEVEL]])
    calculateFlux(U, DX, DT, [level_cells[MAX_REFINEMENT_LEVEL]],
                  [level_faces[MAX_REFINEMENT_LEVEL]])
  end
  return top_level
end

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
task calculateFlux(vel : double,
                   dx : double,
                   dt : double,
                   cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
where
  reads(cells.phi,
        faces.flux),
  writes(faces.flux)
do
  var left_boundary_face : int64 = faces.ispace.bounds.lo
  var right_boundary_face : int64 = faces.ispace.bounds.hi
  C.printf("Faces from %d to %d\n", left_boundary_face, right_boundary_face)
  -- loop on inner faces
  for face = left_boundary_face+1, right_boundary_face do
    var left : double = cells[face - 1].phi
    var right : double = cells[face].phi
    var flux : double = 0.5 * vel * (left + right) +0.5 * dx * (left - right)/dt
    faces[face].flux = flux
  end

  -- boundary conditions: hold end cells constant in time
  faces[0].flux = faces[1].flux
  faces[right_boundary_face].flux = faces[right_boundary_face-1].flux

  -- debug
  for face in faces do
    C.printf("face %d flux %f\n", face, faces[face].flux)
  end
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

local main = make_top_level_task()
main:printpretty()
regentlib.start(main)

