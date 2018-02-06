-- linear advection with Lax-Friedrichs
import "regent"

require("global_const")

-- model specific local constants
local CFL = 0.5
local MAX_NX = 640
local MIN_DX = 1.0 / MAX_NX

-- model specific global constants
U = 1.0
MAX_GRAD = 1.0
MIN_GRAD = 0.1

-- required global constants
DT = CFL * MIN_DX / U

-- model specific fields must be in fspace's CellValues and FaceValues

fspace CellValues
{
  phi : double,
  phi_copy : double
}

fspace FaceValues
{
  flux : double,
  grad : double
}

