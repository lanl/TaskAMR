--Copyright (c) 2018, Los Alamos National Security, LLC
--All rights reserved.

--Copyright 2018. Los Alamos National Security, LLC. This software was produced under
--U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL),
--which is operated by Los Alamos National Security, LLC for the U.S. Department of
--Energy. The U.S. Government has rights to use, reproduce, and distribute this
--software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
--WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
--If software is modified to produce derivative works, such modified software should be
--clearly marked, so as not to confuse it with the version available from LANL.

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
MIN_GRAD = 1.0e-4

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

