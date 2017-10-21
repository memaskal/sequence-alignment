#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <unistd.h> 
#include <sys/time.h>
#include <limits>

#include <omp.h>
#include <x86intrin.h>

#define __SIMD__ 1

// Fine tuning constants
#if __SIMD__ 
	const int I = 128, p = 2;
#else
	const int I = 64, p = 4;
#endif

// Alignment constants
#define _AVX_ALLIGN_ 	32
#define _CACHE_ALLIGN_ 	64

#define _FLP_LINE_  (_CACHE_ALLIGN_ / sizeof(float))
#define _FFLP_REG_  (_AVX_ALLIGN_ / sizeof(float))

#define _MSCORE_ 1.
#define MIN(a,b) ((a)<(b)?(a):(b))

using namespace std;

// Blocking size
int Bw, Bh, N_a, N_b;
const int minBw = 4 * _FLP_LINE_,
		  minBh = 4;

// Some global vars
string seq_a, seq_b;
float mu, delta, **H;
float *seq_b_floats;

// Point type declaration
typedef struct {
	int row, col;
} Point;

// Dynamic alloc functions
void freeMemory(int);
void allocateMemory(int, int);

// Read from file function
string read_sequence_from_file(char *);

// Computational - helper functions
Point movePoint(Point);
void computeBlock(Point);
inline float find_max_score(float []);
inline float similarity_score(char, char);

int main (int argc, char** argv) {

    // Input arguments check
    if( argc != 5 ) {
        printf("Give me the propen number of input arguments:\n");
        printf("1 : mu\n");
        printf("2 : delta\n");
        printf("3 : filename sequence A\n");
        printf("4 : filename sequence B\n");
        exit(1);
    }
    
    mu    = atof(argv[1]);
    delta = atof(argv[2]);

    char *nameof_seq_a = argv[3];
    char *nameof_seq_b = argv[4];
      
	struct timeval StartTime, EndTime;									
    
    // Read the sequences into strings
    seq_a = read_sequence_from_file(nameof_seq_a);
    N_a   = seq_a.length();

	seq_b = read_sequence_from_file(nameof_seq_b);
    N_b   = seq_b.length();

    printf("First sequence has length  : %6d\n", N_a);
    printf("Second sequence has length : %6d\n\n", N_b);

	// Allocate dynamic memory - exit on failure
	allocateMemory(N_a + 1, N_b + 1);
			
	Point upper = { 0, 0 }, 
		  lower = upper, tmp;	
	
	int newRows, newCols, numOfDiags;
	
	// Declare reduction such it will work in serial execution.
	// Maximum values found first would be set as the best score for the algorithm
	struct Comparator { float val; Point p; } max;	
	#pragma omp declare reduction(max : struct Comparator : omp_out = (omp_in.val > omp_out.val || \
	(omp_in.val == omp_out.val && (omp_in.p.row < omp_out.p.row || omp_in.p.col < omp_out.p.col))) ? omp_in : omp_out)
	// Set max staring val as -Inf
	max.val = -numeric_limits<float>::infinity();
			
	// Start counting
    gettimeofday(&StartTime, NULL);

	#if __SIMD__
		const char *char_seqb = seq_b.c_str();
	#endif

	// Create threads
	#pragma omp parallel
    {
		#if __SIMD__
			// Initialize seq_b_flaots to contain seq_b
			#pragma omp for schedule(static, 16)
			for (int i = 0; i < N_b; ++i) {
				seq_b_floats[i] = char_seqb[i];
			}
		#endif

		#pragma omp single
		{		
			// Get the number of working threads
			char T = omp_get_num_threads();
			
			Bh = (int)ceil(((double)N_a)/(T * I));
			Bw = (int)(((double)N_b)/(T * p));
			Bw += _FLP_LINE_ - Bw % _FLP_LINE_;
			
			// If the blocks are really small set a constant size
			if (Bh < minBh) Bh = minBh;
			if (Bw < minBw) Bw = minBw;
			
			// Calculate blocks total
			newRows = (int)ceil(((double)N_a) / Bh);
			newCols = (int)ceil(((double)N_b) / Bw);
			numOfDiags = newRows + newCols - 1;
			
			// Create all blocks
			for (int i = 0; i < numOfDiags; ++i) {
				tmp = lower;
				while (tmp.row >= upper.row && tmp.col <= upper.col) {
					// Create a new task per block
					#pragma omp task firstprivate(tmp)
					computeBlock(tmp);
					--tmp.row; 
					++tmp.col;
				}
				#pragma omp taskwait
				// Execute the increments only when tasks have finished
				if (upper.col < newCols - 1) ++upper.col; else ++upper.row;
				if (lower.row < newRows - 1) ++lower.row; else ++lower.col;
			}
		}
		// Threads wait here, at the implied barrier, while they complete tasks
		
		// Search H for the maximal score
		#pragma omp for reduction(max: max)	
		for (int i = 1; i <= N_a; ++i) {
			for (int j = 1; j <= N_b; ++j) {
				if( H[i][j] > max.val ) {
					max.val = H[i][j];
					max.p.row = i;
					max.p.col = j;
				}
			}
		}
	}
	
    Point curr = { max.p.row, max.p.col },
		  next = movePoint(curr);

    int tick = 0;
    
    char consensus_a[ N_a + N_b + 2 ], 									
		 consensus_b[ N_a + N_b + 2 ];
    
	// Backtracking from H_max
	while (((curr.row!=next.row) || (curr.col!=next.col)) && (next.row!=0) && (next.col!=0)) {
				
		if ( next.row == curr.row ) {
			consensus_a[tick] = '-';                  					// deletion in A
		} else {
			consensus_a[tick] = seq_a[curr.row - 1];   					// match/mismatch in A
		}
		
		if ( next.col == curr.col ) { 
			consensus_b[tick] = '-';                  					// deletion in B
        } else {
			consensus_b[tick] = seq_b[curr.col - 1];   					// match/mismatch in B
		}
	
		++tick;
		curr = next;
		next = movePoint(next);    
    }

    gettimeofday(&EndTime, NULL);


    printf("\n***********************************************\n");
    printf("The alignment of the sequences\n\n");
    for(int i = 0; i < N_a; ++i) {
        printf("%c", seq_a[i]);
    }
    printf("  and\n");
    for(int i = 0; i < N_b; ++i) {
        printf("%c", seq_b[i]);
    }
    printf("\n\nis for the parameters  mu = %.2lf and delta = %.2lf given by\n\n", mu, delta);
    for(int i = tick - 1; i >= 0; --i) {
		printf("%c", consensus_a[i]);
	}
    
    printf("\n");
    for(int j = tick - 1; j >= 0; --j) {
		printf("%c", consensus_b[j]);
	}
    printf("\n");

    if (EndTime.tv_usec < StartTime.tv_usec) {
        int nsec = (StartTime.tv_usec - EndTime.tv_usec) / 1000000 + 1;
        StartTime.tv_usec -= 1000000 * nsec;
        StartTime.tv_sec += nsec;
    }
    if (EndTime.tv_usec - StartTime.tv_usec > 1000000) {
        int nsec = (EndTime.tv_usec - StartTime.tv_usec) / 1000000;
        StartTime.tv_usec += 1000000 * nsec;
        StartTime.tv_sec -= nsec;
    }

    printf("\n\nParallel calculation time: %ld.%.6ld seconds\n", EndTime.tv_sec  - StartTime.tv_sec, EndTime.tv_usec - StartTime.tv_usec);

	// Release dynamic memory
	freeMemory(N_a + 1);
	return 0;
}

/******************************************************************************/
/* Computational - helper functions                                        	 */
/******************************************************************************/

Point movePoint(Point old) {
	
	int row = old.row,
		col = old.col;
	
	float Val = H[row][col];
	
	if ( Val == H[row-1][col-1] + similarity_score(seq_a[row-1], seq_b[col-1])) {
		--old.row;
		--old.col;
	} 
	else if ( Val == H[row-1][col] - delta) {
		--old.row;
	}
	else if ( Val == H[row][col-1] - delta ) {
		--old.col;
	}
	return old;
}


/**
 * Computes the given block int the array
**/
void computeBlock(Point block) {
	
	int i, j;
	const int row = block.row * Bh,
			  col = block.col * Bw; 
	
	int Rows = MIN(N_a, row + Bh),
	    Cols = MIN(N_b, col + Bw);

	float tp, tmp[3];	

	#if __SIMD__
	float constants[] = { _MSCORE_, -mu, delta, 0.f };	
	float __attribute__((aligned(_AVX_ALLIGN_))) max_mem[_FFLP_REG_];

	const int unrolled  = Cols / _FFLP_REG_ * _FFLP_REG_,
			  remainder = Cols - Cols % _FFLP_REG_;

	__m256 	seqa, seqb, mask, diag, up, max,
	mscore = _mm256_broadcast_ss( constants ), 				// matching score
	nscore = _mm256_broadcast_ss(constants+1),				// not-matching score
	dvect  = _mm256_broadcast_ss(constants+2),				// delta vector
	zvect  = _mm256_broadcast_ss(constants+3);				// zero vector
	#else 
	const int remainder = col;
	#endif

	for (i = row + 1; i <= Rows; ++i) 
	{	
		tp = H[i][col];
		#if __SIMD__		
		// Load seq_a to register
		seqa = _mm256_set1_ps(seq_a[i-1]);
		for (j = col + 1; j <= unrolled; j += _FFLP_REG_)
		{					
			// Load seq_b to register
			seqb = _mm256_load_ps(&seq_b_floats[j-1]);

			// Perform the similarity ifs
			mask  = _mm256_cmp_ps(seqa, seqb, _CMP_EQ_OQ);
			diag  = _mm256_or_ps(_mm256_and_ps(mask, mscore),
					    _mm256_andnot_ps(mask, nscore));	

			// Calculate diag, up scores
			diag  = _mm256_add_ps(_mm256_loadu_ps(&H[i-1][j-1]), diag);
			up    = _mm256_sub_ps(_mm256_loadu_ps(&H[i-1][ j ]), dvect);

			// Find max score beetween diag, up, 0.f
			max	  = _mm256_max_ps(diag, _mm256_max_ps(up, zvect));

			// Store max back to memory
			_mm256_store_ps(max_mem, max);

			// Depedent data max here
			H[i][ j ] = tp = (tp - delta > max_mem[0]) ? tp - delta : max_mem[0];
			H[i][j+1] = tp = (tp - delta > max_mem[1]) ? tp - delta : max_mem[1];
			H[i][j+2] = tp = (tp - delta > max_mem[2]) ? tp - delta : max_mem[2];
			H[i][j+3] = tp = (tp - delta > max_mem[3]) ? tp - delta : max_mem[3];
			H[i][j+4] = tp = (tp - delta > max_mem[4]) ? tp - delta : max_mem[4];
			H[i][j+5] = tp = (tp - delta > max_mem[5]) ? tp - delta : max_mem[5];
			H[i][j+6] = tp = (tp - delta > max_mem[6]) ? tp - delta : max_mem[6];
			H[i][j+7] = tp = (tp - delta > max_mem[7]) ? tp - delta : max_mem[7];
		}		
		#endif
		// Loop for the reminder
		for (j = remainder + 1; j <= Cols; ++j)
		{		
			tmp[0] = H[i-1][j-1] + similarity_score(seq_a[i-1], seq_b[j-1]);
			tmp[1] = H[i-1][ j ] - delta;
			tmp[2] = tp - delta;

			// Compute action that produces maximum score
			H[i][j] = tp = find_max_score(tmp);
		}
	}
}


/**
 * Returns the similarity score 
 * 
**/
inline float similarity_score(char a, char b) {
	
    return ((a == b) ? _MSCORE_ : -mu);
}


/**
 * Calculates max between 3 input scores and zero
 * Returns the index of max in values[]
**/
inline float find_max_score(float values[]) {
	 
	float max = (values[0] > values[1]) ? values[0] : values[1];	
	if (values[2] > max) max = values[2];
	return ((max > 0.) ? max : 0.);
}


/******************************************************************************/
/* auxiliary functions used by main                                           */
/******************************************************************************/

/**
 * Allocates memory for the Dynamic Arrays
 * int rows Rows of arrays
 * int cols Columns
**/
void allocateMemory(int rows, int cols) {
	
	int i;

	printf("Allocating memory for matrix H\n");
	H = (float **)aligned_alloc(_CACHE_ALLIGN_, rows * sizeof(float *));
	if (H == NULL) {
		fprintf(stderr, "Could not allocate memory for matrix H\n");
		exit(1);
	}
	for (i = 0; i < rows; ++i) {
		H[i] = (float *)aligned_alloc(_CACHE_ALLIGN_, cols * sizeof(float));
		if (H[i] == NULL) {
		    fprintf(stderr, "Could not allocate memory for matrix H[%6d]\n", i);
		    exit(1);
		}
	}
	printf("Memory for matrix H allocated\n\n");

	printf("Initializing matrix H\n");
	for (i = 0; i < rows; ++i) {
		H[i][0] = 0.;
	}
	for (i = 0; i < rows; ++i) {
		H[0][i] = 0.;
	}        
	printf("Matrix H initialized\n\n");

	#if __SIMD__
		printf("Allocating memory for vector seq_b_floats\n");
		seq_b_floats = (float *)aligned_alloc(_AVX_ALLIGN_, cols * sizeof(float));
		if (seq_b_floats == NULL) {
			fprintf(stderr, "Could not allocate memory for vector seq_b_floats\n");
			exit(1);
		}
		printf("Memory for matrix H allocated\n\n");
	#endif
}


/** 
 * Deallocates dynamic memory
**/
void freeMemory(int rows) {
	
	if ( H ) {
		for ( int i = 0; i < rows; ++i ) {
			if ( H[i] ) free(H[i]);
		}
		free(H);
	}
	
	#if __SIMD__
		if ( seq_b_floats ) {
			free(seq_b_floats);
		}
	#endif	
}


/**
 * Reads sting sequence from file and retturns
 * it to a new string
**/
string read_sequence_from_file(char *filename) {
	
	ifstream file(filename);

	if (file) {
		printf("Opened file \"%s\"\n", filename);
	} else {
		fprintf(stderr, "Error: Can't open the file %s\n", filename);
		exit(1);
	}

	printf("Reading file \"%s\"\n", filename);

	string seq;
	char line[20000];

	while(file.good()) 
	{	
		file.getline(line, 20000);
		if( line[0] == 0 || line[0]== '#' )
			continue;
		for( int i = 0; line[i] != 0; ++i )
		{
			int c = toupper(line[i]);
			if( c != 'A' && c != 'G' && c != 'C' && c != 'T' )
			continue;

			seq.push_back(char(c));
		}
	}
	printf("File \"%s\" read\n\n", filename);
	return seq;
}

