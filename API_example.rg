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

import "regent"
local C = regentlib.c

-- required AMR global constants

CELLS_PER_BLOCK_X = -- must be multiple of 2
LEVEL_1_BLOCKS_X = -- number of blocks at coarsest level
MAX_REFINEMENT_LEVEL = -- every level above 1 doubles max resolution
NUM_PARTITIONS = -- how many parallel pieces to break problem into (suggestion)
LENGTH_X = -- DX = LENGTH_X / NX
DT = -- currently time steps are fixed
T_FINAL == -- simulation ends at T_FINAL <= time < T_FINAL + DT

-- model specific fields must be in fspace's CellValues and FaceValues

fspace CellValues
{
}

fspace FaceValues
{
}

-- model specific tasks must implement following API:
-- task initializeCells(cell_region: region(ispace(int1d), CellValues))
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

