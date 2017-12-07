-- 1D Euler equations with Lax-Friedrichs method

import "regent"
local C = regentlib.c

-- required AMR global constants
CELLS_PER_BLOCK_X = 5
LEVEL_1_BLOCKS_X = 5
MAX_REFINEMENT_LEVEL = 1

-- constants that should not exist
local NX = 25
local MAX_NX = 3200

-- model specific global constants
DX = 1.0 / NX
local MIN_DX = 1.0 / MAX_NX
T_FINAL = 0.142681382
DT = 0.2 * MIN_DX  -- dt < dx / (2^0.5 * (u+c))

-- model specific local constants
local GAMMA = 1.4
local BETA = (GAMMA + 1.0) / (GAMMA - 1.0)
-- left and right initial states for Sod shock
local P_L = 1.0
local RHO_L = 1.0
local V_L = 0.0
local P_R = 0.1
local RHO_R = 0.125
local V_R = 0.0

-- model specific fields

fspace CellValues
{
  density : double,
  velocity : double,
  momentum : double,
  pressure : double,
  energy : double,
}

fspace FaceValues
{
  density_flux : double,
  momentum_flux : double,
  energy_flux : double,
}

-- model specific tasks

task initializeCells(cell_region: region(ispace(int1d), CellValues))
where
  writes(cell_region.{density,
                      velocity,
                      momentum,
                      pressure,
                      energy})
do
  var size : int64 = cell_region.ispace.bounds.hi - cell_region.ispace.bounds.lo + 1
  for cell in cell_region.ispace do
    var P : double
    if [int64](cell) < (size/2) then
      P = P_L
      cell_region[cell].density = RHO_L
      cell_region[cell].velocity = V_L
      cell_region[cell].momentum = RHO_L * V_L
    else
      P = P_R
      cell_region[cell].density = RHO_R
      cell_region[cell].velocity = V_R
      cell_region[cell].momentum = RHO_R * V_R
    end
    cell_region[cell].pressure = P
    -- this should be meta-programmed
    cell_region[cell].energy = P / (GAMMA - 1.0)
  end
  C.printf("initializeCells %d cells\n", size)
end

-- hard coded euler with Lax-Friedrichs for now, metaprog soon
task applyFlux(dx : double,
               dt : double,
               cells: region(ispace(int1d), CellValues),
               faces: region(ispace(int1d), FaceValues))
where
  reads(cells.{density,
               momentum,
               energy},
        faces.{density_flux,
               momentum_flux,
               energy_flux}),
  writes(cells.{density,
               momentum,
               energy,
               velocity,
               pressure})
do
  for cell in cells do
    cells[cell].density = cells[cell].density
             - dt * (faces[cell+1].density_flux - faces[cell].density_flux) / dx
    cells[cell].momentum = cells[cell].momentum
             - dt * (faces[cell+1].momentum_flux - faces[cell].momentum_flux) / dx
    cells[cell].energy = cells[cell].energy
             - dt * (faces[cell+1].energy_flux - faces[cell].energy_flux) / dx
    var velocity : double =  cells[cell].momentum /  cells[cell].density
    cells[cell].velocity =  velocity
    -- this should be meta programmed
    cells[cell].pressure = (cells[cell].energy - 0.5 * cells[cell].momentum
                          * velocity) * (GAMMA - 1.0)
  end
end

-- hard coded euler Lax-Friedrichs for now, metaprog maybe
task calculateFlux(dx : double,
                   dt : double,
                   cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
where
  reads(cells.{density,
               velocity,
               momentum,
               pressure,
               energy},
        faces.{density_flux,
              momentum_flux,
              energy_flux}),
  writes(faces.{density_flux,
                momentum_flux,
                energy_flux})
do
  var left_boundary_face : int64 = faces.ispace.bounds.lo
  var right_boundary_face : int64 = faces.ispace.bounds.hi
  -- loop on inner faces
  for face = left_boundary_face+1, right_boundary_face do

    var rho_l : double = cells[face - 1].density
    var rho_r : double = cells[face].density

    var v_l : double = cells[face - 1].velocity
    var v_r : double = cells[face].velocity

    var rhov_l : double = cells[face - 1].momentum
    var rhov_r : double = cells[face].momentum

    var P_l : double = cells[face - 1].pressure
    var P_r : double = cells[face].pressure

    var E_l : double = cells[face - 1].energy
    var E_r : double = cells[face].energy

    var density_flux : double = 0.5 * (rhov_l + rhov_r)
                                + 0.5 * dx * (rho_l - rho_r)/dt
    faces[face].density_flux = density_flux

    var momentum_flux : double = 0.5 * (rhov_l * v_l + P_l + rhov_r * v_r + P_r)
                                + 0.5 * dx * (rhov_l - rhov_r)/dt
    faces[face].momentum_flux = momentum_flux

    var energy_flux : double = 0.5 * (v_l * (P_l + E_l) + v_r * (P_r + E_r)) 
                                + 0.5 * dx * (E_l - E_r)/dt
    faces[face].energy_flux = energy_flux
  end

  -- boundary conditions: hold end cells constant in time
  faces[0].density_flux = faces[1].density_flux
  faces[0].momentum_flux = faces[1].momentum_flux
  faces[0].energy_flux = faces[1].energy_flux
  faces[right_boundary_face].density_flux = faces[right_boundary_face-1].density_flux
  faces[right_boundary_face].momentum_flux = faces[right_boundary_face-1].momentum_flux
  faces[right_boundary_face].energy_flux = faces[right_boundary_face-1].energy_flux
end

task writeCells(cells: region(ispace(int1d), CellValues))
where
  reads(cells.{density,
               momentum,
               energy})
do
  var first_cell : int64 = cells.ispace.bounds.lo
  var last_cell : int64 = cells.ispace.bounds.hi
  --var nx : int64 = last_cell - first_cell + 1
  var fp = C.fopen(["euler." .. NX .. ".txt"],"w")
  for cell in cells do
    C.fprintf(fp, "%f %f %f\n", cells[cell].density, cells[cell].momentum,
              cells[cell].energy)
  end
  C.fclose(fp)
end

