# CloverLeaf_Nested Parallelism

This is the nested parallelism version of CloverLeaf version 1.3. 

## Release Notes

### Version 1.3

CloverLeaf 1.3 contains a number of optimisations over previous releases.
These include a number of loop fusion optimisations and the use of scalar variables over work arrays.
Overall this improves cache efficiency.

This version also contains some support for explicit tiling.
This is activated through the two input deck parameters:

* `tiles_per_chunk` To specify how many tiles per MPI ranks there are.
* `tiles_per_problem` To specify how many global tiles there are, this is rounded down to be an even number per MPI rank.


### Nested Parallelism

Now with the inclusion of 'tiles' we can have an additional level of parallelism.
Firstly we have parallelism at the tile decomposition level, and again at the loop level.

Running this version without nested parallelism will default to the top level of parallelism, the tiles.
This means there must be at least as many tiles as there are OpenMP threads.
This can be enabled through the 'tiles_per_chunk' input parameter.

#### Running Nested - Intel

We advise the following set up to make use of nested parallelism

```
    KMP_HOT_TEAMS_MODE=1 \
    KMP_HOT_TEAMS_MAX_LEVEL=2 \
    OMP_NESTED=1 \
    OMP_NUM_THREADS=4,2 \
    mpirun -np 2 ./clover_leaf
```
Giving a total thread count of 16 (2 MPI each with 4 OMP each with 2 OMP).


## Performance

Expected performance is give below.

If you do not see this performance, or you see variability, then is it recommended that you check MPI task placement and OpenMP thread affinities, because it is essential these are pinned and placed optimally to obtain best performance.

Note that performance can depend on compiler (brand and release), memory speed, system settings (e.g. turbo, huge pages), system load etc. 

### Performance Table

| Test Problem  | Time                         | Time                        | Time                        |
| ------------- |:----------------------------:|:---------------------------:|:---------------------------:|
| Hardware      |  E5-2670 0 @ 2.60GHz Core    | E5-2670 0 @ 2.60GHz Node    | E5-2698 v3 @ 2.30GHz Node   |
| Options       |  make COMPILER=INTEL         | make COMPILER=INTEL         | make COMPILER=CRAY          |
| Options       |  mpirun -np 1                | mpirun -np 16               | aprun -n4 -N4 -d8           |
| 2             | 20.0                         | 2.5                         | 0.9                         |
| 3             | 960.0                        | 100.0                       |                             |
| 4             | 460.0                        | 40.0                        | 23.44                       |
| 5             | 13000.0                      | 1700.0                      |                             |

### Weak Scaling - Test 4

| Node Count | Time         |
| ---------- |:------------:|
| 1          |   40.0       |
| 2          |              |
| 4          |              |
| 8          |              |
| 16         |              |


### Strong Scaling - Test 5

| Node Count | Time          | Speed Up |
| ---------- |:-------------:|:--------:|
| 1          |   1700        |  1.0     |
| 2          |               |          |
| 4          |               |          |
| 8          |               |          |
| 16         |               |          |


