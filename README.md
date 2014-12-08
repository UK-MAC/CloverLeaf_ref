# CloverLeaf_ref

This is the reference version of CloverLeaf. Expected performance is give below.

If you do not see this performance, or you see variability, then is it recommended that you check MPI task placement and OpenMP thread affinities, because it is essential these are pinned and placed optimally to obtain best performance.

### Performance Table

| Test Problem  | Time                         | Time                        |
| ------------- |:----------------------------:|:---------------------------:|
| Hardware      |  E5-2670 0 @ 2.60GHz Core    | E5-2670 0 @ 2.60GHz Node    |
| Options       |  make COMPILER=INTEL         | mpirun -np 16               |
| 2             | 20.0                         | 2.5                         |
| 3             | 960.0                        | 100.0                       |
| 4             | 460.0                        | 40.0                        |
| 5             | 13000.0                      | 1700.0                      |

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


