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
-- 1D Euler equations with Lax-Friedrichs method

import "regent"
local C = regentlib.c

require("global_const")
require("refinement_bits")

-- model specific local constants
local MAX_NX = 3200
local MIN_DX = 1.0 / MAX_NX

-- required global constants
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

task initializeCells(num_cells : int64,
                     cell_region: region(ispace(int1d), CellValues))
where
  writes(cell_region.{density,
                      velocity,
                      momentum,
                      pressure,
                      energy})
do
  for cell in cell_region.ispace do
    var P : double
    if [int64](cell) < (num_cells/2) then
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
  C.printf("initializeCells %d cells\n", num_cells)
end

-- hard coded euler with Lax-Friedrichs for now, metaprog soon
task applyFlux(dx : double,
               dt : double,
               blocks: region(ispace(int1d), RefinementBits),
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
  var face_index : int64 = faces.ispace.bounds.lo
  for cell in cells do
    cells[cell].density = cells[cell].density
             - dt * (faces[face_index+1].density_flux - faces[face_index].density_flux) / dx
    cells[cell].momentum = cells[cell].momentum
             - dt * (faces[face_index+1].momentum_flux - faces[face_index].momentum_flux) / dx
    cells[cell].energy = cells[cell].energy
             - dt * (faces[face_index+1].energy_flux - faces[face_index].energy_flux) / dx
    var velocity : double =  cells[cell].momentum /  cells[cell].density
    cells[cell].velocity =  velocity
    -- this should be meta programmed
    cells[cell].pressure = (cells[cell].energy - 0.5 * cells[cell].momentum
                          * velocity) * (GAMMA - 1.0)
    face_index = face_index + 1
  end
end

-- hard coded euler Lax-Friedrichs for now, metaprog maybe
task calculateFlux(num_cells : int64,
                   dx : double,
                   dt : double,
                   blocks: region(ispace(int1d), RefinementBits),
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

  var left_boundary_cell : int64 = cells.ispace.bounds.lo
  var right_boundary_cell : int64 = cells.ispace.bounds.hi

  var start_face : int64 = left_boundary_face
  var stop_face : int64 = right_boundary_face + 1

  if left_boundary_cell == 0 then
    start_face = start_face + 1
  end
  if right_boundary_cell == (num_cells - 1) then
    stop_face = stop_face - 1
  end

  -- loop on inner faces
  var cell_index : int64 = left_boundary_cell
  for face = start_face, stop_face do

    var rho_l : double = cells[cell_index].density
    var rho_r : double = cells[cell_index + 1].density

    var v_l : double = cells[cell_index].velocity
    var v_r : double = cells[cell_index + 1].velocity

    var rhov_l : double = cells[cell_index].momentum
    var rhov_r : double = cells[cell_index + 1].momentum

    var P_l : double = cells[cell_index].pressure
    var P_r : double = cells[cell_index + 1].pressure

    var E_l : double = cells[cell_index].energy
    var E_r : double = cells[cell_index + 1].energy

    var density_flux : double = 0.5 * (rhov_l + rhov_r)
                                + 0.5 * dx * (rho_l - rho_r)/dt
    faces[face].density_flux = density_flux

    var momentum_flux : double = 0.5 * (rhov_l * v_l + P_l + rhov_r * v_r + P_r)
                                + 0.5 * dx * (rhov_l - rhov_r)/dt
    faces[face].momentum_flux = momentum_flux

    var energy_flux : double = 0.5 * (v_l * (P_l + E_l) + v_r * (P_r + E_r)) 
                                + 0.5 * dx * (E_l - E_r)/dt
    faces[face].energy_flux = energy_flux

    cell_index = cell_index + 1
  end

  -- boundary conditions: hold end cells constant in time
  if left_boundary_cell == 0 then
    faces[0].density_flux = faces[1].density_flux
    faces[0].momentum_flux = faces[1].momentum_flux
    faces[0].energy_flux = faces[1].energy_flux
  end
  if right_boundary_cell == (num_cells - 1) then
    faces[right_boundary_face].density_flux = faces[right_boundary_face-1].density_flux
    faces[right_boundary_face].momentum_flux = faces[right_boundary_face-1].momentum_flux
    faces[right_boundary_face].energy_flux = faces[right_boundary_face-1].energy_flux
  end
end

task writeCells(nx : int64,
                cells: region(ispace(int1d), CellValues))
where
  reads(cells.{density,
               momentum,
               energy})
do
  var first_cell : int64 = cells.ispace.bounds.lo
  var last_cell : int64 = cells.ispace.bounds.hi
  var buf : &int8
  buf = [&int8](C.malloc(40))
  C.sprintf(buf, "euler.%d.txt", nx)
  var fp = C.fopen(buf ,"w")
  for cell in cells do
    C.fprintf(fp, "%f %f %f\n", cells[cell].density, cells[cell].momentum,
              cells[cell].energy)
  end
  C.fclose(fp)
  C.free([&opaque](buf))
end

