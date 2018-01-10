Regent *executables* include:
  1d_fix.rg: 1 dimensional fixed grid (at MAX_REFINEMENT_LEVEL resolution)

To use, implement or link "model.rg" to fit the following API:

## required global constants

```
CELLS_PER_BLOCK_X = -- must be multiple of 2
LEVEL_1_BLOCKS_X = -- number of blocks at coarsest level
MAX_REFINEMENT_LEVEL = -- every level above 1 doubles max resolution
NUM_PARTITIONS = -- how many parallel pieces to break problem into (suggestion)
LENGTH_X = -- DX = LENGTH_X / NX
DT = -- currently time steps are fixed
T_FINAL == -- simulation ends at T_FINAL <= time < T_FINAL + DT
```
## model specific fields must be in fspace's CellValues and FaceValues

```
fspace CellValues
{
}

fspace FaceValues
{
}
```

## model specific tasks must implement following API:
```
task initializeCells(cell_region: region(ispace(int1d), CellValues))
task calculateFlux(num_cells : int64,
                   dx : double,
                   dt : double,
                   cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))
task applyFlux(dx : double,
               dt : double,
               cells: region(ispace(int1d), CellValues),
               faces: region(ispace(int1d), FaceValues))
task writeCells(num_cells : int64,
                cells: region(ispace(int1d), CellValues))
```
Then, regent.py 1d_fix.rg

