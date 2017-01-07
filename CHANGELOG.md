# Change Log

## Unreleased
- Added residue eta calculation (average eta per residue).
- New struct-based API.
- Added -nthreads option to set the number of threads to use at runtime.
- Fixed memory leaks in svm training.
- Atom IDs are output next to eta values.
- Shows progress while constructing svm problem structs (can take a long time for large trajectories).

## 2.0.1 - 2016-02-11
- Added libdl to link stage which is necessary in some environments.

## 2.0.0 - 2015-10-03
- Added OpenMP parallelization, with significant improvement in performance on multicore/multiprocessor systems.
- Build outputs now go in build directory.

## 1.0.1 - 2015-09-22
- Added MAE of method (3.26%) to description.

## 1.0.0 - 2015-09-15
- First stable release.
