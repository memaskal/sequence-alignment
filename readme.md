# Smith Waterman - local sequence alignment

A simple parallel implementation of the [Smith Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)'s algorithm, created as a university assignment. Parallelism was achieved using the [OpenMp](http://openmp.org/wp/) API for multi-threaded applications and AVX vectorization for finer tuning. 

## Algorithm Design

This implementation uses the wave-front method, for anti-diagonal calculation of the score matrix. The last, is a `NxM` float matrix, where `N,M` are the lengths of the two input sequences. The score matrix calculation, happens in blocks of variable size. The blocking size, is determined by the number of threads available, as well as the sequences' size and other fine-tuned constants defined through benchmarking on our machine.

The computational phase follows a consumer producer pattern. The master thread, produces OpenMP Tasks for each block in the score matrix, which are later consumed by working threads using AVX vectorization instructions, when applicable, to speed-up the calculations. For more information on the decisions of this design and overall walk-through on the development check the project's [report](./report.pdf).

## Results

Each of these implementation were tested with no optimization `-O0` and all optimizations `-O3` flags, using one to four threads (systems maximum). The table bellow contains the speed-up measurement against the simple serial implementation:

| # Threads    |   1 	|   2	|   4 |
| ---	 |  :-:	    |:-:	| :-:  |	
| OpenMP-O0  	    | 1.15 | 2.12 | 2.63 |
| OpenMP-O3 	    | 1.58 | 3.00 | 3.67 | 
| OpenMP-O0 + AVX  	| 2.93 | 5.41 | 6.74 |
| OpenMP-O3 + AVX  	| 3.48 | 6.18 | 8.14 |

## Testing

Input test sequences are provided, in the [data](./data) directory of this project. To compare the execution of the original and the parallel implementation, simply run the `compare.sh` script from the root directory, of this project. This script, takes as argument a number #`0-2`, which represent small, medium and a large test cases.  

``` shell
# Test on a medium length sequence
$ ./compare.sh 1
```

## Compiling 

To compile each implementation, simply run the `make` command. The build requirements are listed below:
* g++ compiler with OpenMp v4.0 or greater
* AVX enabled CPU

> Troubleshooting: setting the `__SIMD__` flag to zero, will resolve vectorization related errors, on systems which don't support [Intel AVX Intrinsics](https://software.intel.com/sites/landingpage/IntrinsicsGuide/#techs=AVX).

## Disclaimer 

This is not a "great" implementation to be used for real purposes. There are really faster and more  flexible implementations out there. This is just a university project, heavily optimized to run on a specific machine. 

