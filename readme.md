# Smith Waterman - local sequence alignment

This project was a university assignment on the course of parallel processing. Given the serial [Smith Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) algorithm, a parallel implementation was created preserving the same input-output constraints. The parallelism was achieved using the [OpenMp](http://openmp.org/wp/) API for multi-threaded applications and vectorization for finner tunning. 

## Algorithm Design

This implementation uses the wavefront method, for anti-diagonal calculation of the score matrix. The last, is a `NxM` sized array of floats, where `N,M` are the sizes of the input sequences. The calculation of the score matrix happens in blocks, of variable size. The blocking size is determined  by the number of threads available, as well as sequence sizes and other constants calculated through benchmarking. Each block is computed by a single thread, using OpenMP Tasks and vectorization instructions when applicable. Score matrix calculation, max sequence value and backtracking operations, have been altered for best performance. For more information on the desciscions of this design and overall walkthrought on the development check the `report.pdf`, currently in greek.

## Usage

The `data` directory includes some test input sequences for the application. Simply run the `compare.sh` script on the main directory of this repo, for the #tests `0-2`, to benchmark the two implementations. 
``` shell
$ ./compare.sh 1
```
Troubleshooting: setting the `__SIMD__` flag to zero, will resolve vectorization related errors, on systems which don't support [Intel AVX Intrinsics](https://software.intel.com/sites/landingpage/IntrinsicsGuide/#techs=AVX). The build requirements are listed below:

* g++ compiler with OpenMp v4.0 or greater
* AVX enabled CPU

## Disclaimer 

This is not a "great" implementation to be used for real purposes. There are really faster and more  flexible implementations out there. This is just a university project, heavily optimized to run on a specific machine. 

