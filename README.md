# TaskAMR
Adaptive Mesh Refinement for Finite Volumes using Legion

## License

Please read the [License](./LICENSE).

## Installation instructions

1. At minimum, you need a C compiler, CLang (that includes LLVM), and MPI.
  * I used gcc/4.9.3, clang/3.7.0, and openmpi/1.10.0-gcc_4.9.3
2. `git clone git@gitlab.com:StanfordLegion/legion.git`
3. `cd legion/language/`
4. `./install.py --debug`
  * For performance runs, repeat without the `--debug`.
5. Test your installation with `./regent.py ./examples/circuit.rg`.

## 1D fixed-grid linear advection

1. Setup `model.rg`
```
ln -sf linear_advection.rg model.rg
```
2. Run the model
```
mpirun -n 4 <PATH_TO>/regent.py ./1d_fix.rg 
```
3. Plot the final time result
```
./analyze_linear.py linear.80.txt
```
 * Note that output files from multiple resolutions can be given in sequence to measure the error convergence.

### Model configuration

#### Global constants
`global_const.rg` requires the settings:

```
CELLS_PER_BLOCK_X = -- must be multiple of 2 for AMR
LEVEL_1_BLOCKS_X = -- number of blocks at coarsest level
MAX_REFINEMENT_LEVEL = -- every level above 1 doubles max resolution
NUM_PARTITIONS = -- how many parallel pieces to break problem into (suggestion)
LENGTH_X = -- DX = LENGTH_X / NX
T_FINAL == -- simulation ends at T_FINAL <= time < T_FINAL + DT
```
These settings are shared with the AMR version.  For fix-grid calculations, the resolution is fixed at
`CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X * 2 ** (MAX_REFINEMENT_LEVEL - 1)`.

#### Linear model constants
`linear_constants.rg` requires the settings:

```
U = -- velocity for linear advection
DT = -- fixed time step
fspace CellValues = -- Do not change these
fspace FaceValues = -- Do not change these
```

#### Initial conditions

In `linear_advection.rg`, the task `initializeCells()` can be altered to change the initial conditions.



## 1D fixed-grid Euler equations

1. Setup `model.rg`
```
ln -sf euler.rg model.rg
```
2. Change simulation run time in `global_const.rg`:
```
T_FINAL = 0.142681382
```
3. Run the model
```
mpirun -n 4 <PATH_TO>/regent.py ./1d_fix.rg 
```
4. Change simulation resolution in `global_const.rg`:
```
CELLS_PER_BLOCK_X = 32
```
5. Run the model
```
mpirun -n 4 <PATH_TO>/regent.py ./1d_fix.rg 
```
6. Plot the final time results
```
./analyze_euler.py euler.80.txt euler.1280.txt
```

### Model configuration

#### Global constants
`global_const.rg` settings are the same as for linear advection.

#### Euler model constants
`euler.rg` requires the settings:

```
GAMMA = -- ratio of specific heats for ideal gas equation of state
DT = -- fixed time step
fspace CellValues = -- Do not change these
fspace FaceValues = -- Do not change these
```

#### Initial conditions

In `euler.rg`, the task `initializeCells()` can be altered to change the initial conditions.


## 1D AMR-grid linear advection

1. Setup `model.rg` and `model_amr.rg`
```
ln -sf linear_advection.rg model.rg
ln -sf linear_advection_amr.rg model_amr.rg
```
2. Run the model
```
mpirun -n 4 <PATH_TO>/regent.py ./1d_amr.rg 
```
3. Plot the final time result
```
./analyze_amr_linear.py linear_amr.*.txt
```

### Model configuration

#### Global constants
`global_const.rg` settings are the same as for fixed-grid linear advection
with the exception that `LENGTH_X / (CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X * 2 ** (MAX_REFINEMENT_LEVEL - 1))` is now the minumum grid size instead
of the fixed grid size.

#### Linear model constants
`linear_constants.rg` settings are the same as for fixed-grid linear advection.

#### Initial conditions

Initial conditions settings are the same as for fixed-grid linear advection.

## Tests

### Convergence tests

To run the linear fixed-grid convergence test:
```
./test_linear.py
```

To run the euler fixed-grid convergence test:
```
./test_euler.py
```

To run the linear AMR-grid convergence test:
```
./test_linear_amr.py
```

### Unit tests

To run the unit tests for AMR grid refinement and coarsening:
```
mpirun -n 4 <PATH_TO>/regent.py ./unit_tests.rg 
```

## Create a new physics model for 1D fixed-grid

#### Global constants
`global_const.rg` is required:

```
CELLS_PER_BLOCK_X = -- must be multiple of 2
LEVEL_1_BLOCKS_X = -- number of blocks at coarsest level
MAX_REFINEMENT_LEVEL = -- every level above 1 doubles max resolution
NUM_PARTITIONS = -- how many parallel pieces to break problem into (suggestion)
DT = -- fixed time step
LENGTH_X = -- DX = LENGTH_X / NX
T_FINAL == -- simulation ends at T_FINAL <= time < T_FINAL + DT
```
These settings are shared with the AMR version.  For fix-grid calculations, the resolution is fixed at
`CELLS_PER_BLOCK_X * LEVEL_1_BLOCKS_X * 2 ** (MAX_REFINEMENT_LEVEL - 1)`.

### Model specific fields must be in fspace's CellValues and FaceValues

```
fspace CellValues
{
}

fspace FaceValues
{
}
```

### Model specific tasks

`model.rg` must implement following API:
```
task initializeCells(num_cells : int64,
                     cell_region: region(ispace(int1d), CellValues))

task calculateFlux(num_cells : int64,
                   dx : double,
                   dt : double,
                   blocks: region(ispace(int1d), RefinementBits),
                   bloated_cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))


task applyFlux(dx : double,
               dt : double,
               blocks: region(ispace(int1d), RefinementBits),
               cells: region(ispace(int1d), CellValues),
               faces: region(ispace(int1d), FaceValues))

task writeCells(nx : int64,
                cells: region(ispace(int1d), CellValues))
```


## Create a new physics model for 1D AMR-grid

1. Follow all steps to create a new physics model for 1D fixed-grid except that you
will not need `writeCells()`.
2. `model_amr.rg` must implement following API:
```
task calculateGradient(num_cells : int64,
                   dx : double,
                   bloated_cells: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))

task flagRegrid(blocks: region(ispace(int1d), RefinementBits),
                faces: region(ispace(int1d), FaceValues))

task interpolateToChildren(num_children: int64,
                           blocks: region(ispace(int1d), RefinementBits),
                           ghost_parents: region(ispace(int1d), CellValues),
                           ghost_children: region(ispace(int1d), CellValues),
                           children: region(ispace(int1d), CellValues))

task smoothGrid(blocks: region(ispace(int1d), RefinementBits))

task updateRefinement(num_blocks: int64,
                      blocks: region(ispace(int1d), RefinementBits),
                      cells: region(ispace(int1d), CellValues),
                      ghosts: region(ispace(int1d), RefinementBits),
                      children: region(ispace(int1d), RefinementBits),
                      child_cells: region(ispace(int1d), CellValues),
                      ghost_children: region(ispace(int1d), RefinementBits))

task writeAMRCells(ncells : int64,
                   blocks: region(ispace(int1d), RefinementBits),
                   cells: region(ispace(int1d), CellValues))

task printAMRCells(level : int64,
                   blocks: region(ispace(int1d), RefinementBits),
                   cells: region(ispace(int1d), CellValues))

task copyToChildren(blocks: region(ispace(int1d), RefinementBits),
                    cells: region(ispace(int1d), CellValues),
                    children: region(ispace(int1d), CellValues))

task calculateAMRFlux(num_cells : int64,
                   dx : double,
                   dt : double,
                   blocks: region(ispace(int1d), RefinementBits),
                   bloated_cells: region(ispace(int1d), CellValues),
                   bloated_children: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))

task calculateAMRGradient(num_cells : int64,
                   dx : double,
                   blocks: region(ispace(int1d), RefinementBits),
                   bloated_cells: region(ispace(int1d), CellValues),
                   bloated_children: region(ispace(int1d), CellValues),
                   faces: region(ispace(int1d), FaceValues))


```


