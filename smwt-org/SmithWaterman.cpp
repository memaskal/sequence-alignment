#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>
using namespace std;

float  similarity_score(char a, char b);
float  find_array_max(float array[], int length);
void   checkfile(int open, char filename[]);
string read_sequence(ifstream& f);

int ind;
unsigned short **I_i, **I_j;
float mu, delta, **H;

struct timeval	StartTime, EndTime;

int main(int argc, char** argv) {

    float temp[4];

    // read info from arguments
    if(argc!=5) {
        cout<<"Give me the propen number of input arguments:"<<endl<<"1 : mu"<<endl;
        cout<<"2 : delta"<<endl<<"3 : filename sequence A"<<endl<<"4 : filename sequence B"<<endl;
        exit(1);
    }

    mu    = atof(argv[1]);
    delta = atof(argv[2]);

    char *nameof_seq_a = argv[3];
    char *nameof_seq_b = argv[4];
    string seq_a,seq_b;

    // read the sequences into two vectors:
    ifstream stream_seq_a;
    stream_seq_a.open(nameof_seq_a);
    checkfile(! stream_seq_a,nameof_seq_a);
    cout << "Reading file \"" << nameof_seq_a << "\"\n";
    seq_a = read_sequence(stream_seq_a);
    cout << "File \"" << nameof_seq_a << "\" read\n\n";

    ifstream stream_seq_b;
    stream_seq_b.open(nameof_seq_b);
    checkfile(! stream_seq_b,nameof_seq_b);
    cout << "Reading file \"" << nameof_seq_b << "\"\n";
    seq_b = read_sequence(stream_seq_b);
    cout << "File \"" << nameof_seq_b << "\" read\n\n";

    int N_a = seq_a.length();
    int N_b = seq_b.length();

    cout << "First sequence has length  : " << setw(6) << N_a <<endl;
    cout << "Second sequence has length : " << setw(6) << N_b << endl << endl;


    cout << "Allocating memory for matrix H\n";
    H = (float **)malloc((N_a + 1) * sizeof(float *));
    if (H == NULL) {
        cout << "Could not allocate memory for matrix H\n";
        exit(1);
    }
    for (int i = 0; i < (N_a + 1); i++) {
        H[i] = (float *)malloc((N_b + 1) * sizeof(float));
        if (H[i] == NULL) {
            cout << "Could not allocate memory for matrix H[" << setw(6) << i << "]\n";
            exit(1);
        }
    }
    cout << "Memory for matrix H allocated\n\n";

    cout << "Initializing matrix H\n";
    for(int i=0; i<=N_a; i++) {
        for(int j=0; j<=N_b; j++) {
            H[i][j]=0.;
        }
    }
    cout << "Matrix H initialized\n\n";

    cout << "Allocating memory for matrix I_i\n";
    I_i = (unsigned short **)malloc((N_a + 1) * sizeof(unsigned short *));
    if (I_i == NULL) {
        cout << "Could not allocate memory for matrix I_i\n";
        exit(1);
    }
    for (int i = 0; i < (N_a + 1); i++) {
        I_i[i] = (unsigned short *)malloc((N_b + 1) * sizeof(unsigned short));
        if (I_i[i] == NULL) {
            cout << "Could not allocate memory for matrix I_i[" << setw(6) << i << "]\n";
            exit(1);
        }
    }
    cout << "Memory for matrix I_i allocated\n\n";

    cout << "Allocating memory for matrix I_j\n";
    I_j = (unsigned short **)malloc((N_a + 1) * sizeof(unsigned short *));
    if (I_j == NULL) {
        cout << "Could not allocate memory for matrix I_j\n";
        exit(1);
    }
    for (int i = 0; i < (N_a + 1); i++) {
        I_j[i] = (unsigned short *)malloc((N_b + 1) * sizeof(unsigned short));
        if (I_j[i] == NULL) {
            cout << "Could not allocate memory for matrix I_j[" << setw(6) << i << "]\n";
            exit(1);
        }
    }
    cout << "Memory for matrix I_j allocated\n\n";

    // here comes the actual algorithm

    gettimeofday(&StartTime, NULL);

    for(unsigned short i=1; i<=N_a; i++) {
        for(unsigned short j=1; j<=N_b; j++) {
            temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1]);
            temp[1] = H[i-1][j]-delta;
            temp[2] = H[i][j-1]-delta;
            temp[3] = 0.;
            H[i][j] = find_array_max(temp, 4);
            switch(ind) {
            case 0:                                  // score in (i,j) stems from a match/mismatch
                I_i[i][j] = i-1;
                I_j[i][j] = j-1;
                break;
            case 1:                                  // score in (i,j) stems from a deletion in sequence A
                I_i[i][j] = i-1;
                I_j[i][j] = j;
                break;
            case 2:                                  // score in (i,j) stems from a deletion in sequence B
                I_i[i][j] = i;
                I_j[i][j] = j-1;
                break;
            case 3:                                  // (i,j) is the beginning of a subsequence
                I_i[i][j] = i;
                I_j[i][j] = j;
                break;
            }
        }
    }

#if 0
    // Print the matrix H to the console
    cout<<"**********************************************"<<endl;
    cout<<"The scoring matrix is given by  "<<endl<<endl;
    for(int i=1; i<=N_a; i++) {
        for(int j=1; j<=N_b; j++) {
            cout<< H[i][j]<<"\t";
        }
        cout<<endl;
    }
#endif

    // search H for the maximal score
    float H_max = 0.;
    int i_max=0,j_max=0;
    for(int i=1; i<=N_a; i++) {
        for(int j=1; j<=N_b; j++) {
            if(H[i][j]>H_max) {
                H_max = H[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }

    // Backtracking from H_max
    int current_i = i_max, current_j = j_max;
    int next_i = I_i[current_i][current_j];
    int next_j = I_j[current_i][current_j];
    int tick = 0;
    char consensus_a[N_a+N_b+2], consensus_b[N_a+N_b+2];

    while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)) {

        if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
        else                   consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A

        if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
        else                   consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B

        current_i = next_i;
        current_j = next_j;
        next_i = I_i[current_i][current_j];
        next_j = I_j[current_i][current_j];
        tick++;
    }

    gettimeofday(&EndTime, NULL);

    cout<<endl<<"***********************************************"<<endl;
    cout<<"The alignment of the sequences"<<endl<<endl;
    for(int i=0; i<N_a; i++) {
        cout<<seq_a[i];
    };
    cout<<"  and"<<endl;
    for(int i=0; i<N_b; i++) {
        cout<<seq_b[i];
    };
    cout<<endl<<endl;
    cout<<"is for the parameters  mu = "<<mu<<" and delta = "<<delta<<" given by"<<endl<<endl;
    for(int i=tick-1; i>=0; i--) cout<<consensus_a[i];
    cout<<endl;
    for(int j=tick-1; j>=0; j--) cout<<consensus_b[j];
    cout<<endl;

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

}

/******************************************************************************/
/* auxiliary functions used by main                                           */
/******************************************************************************/

void checkfile(int open, char filename[]) {

    if (open) {
        cout << "Error: Can't open the file "<<filename<<endl;
        exit(1);
    }
    else cout<<"Opened file \"" << filename << "\"\n";
}

/******************************************************************************/

float similarity_score(char a,char b) {

    float result;
    if(a==b) {
        result=1.;
    }
    else {
        result=-mu;
    }
    return result;
}

/******************************************************************************/

float find_array_max(float array[],int length) {

    float max = array[0];            // start with max = first element
    ind = 0;

    for(int i = 1; i<length; i++) {
        if(array[i] > max) {
            max = array[i];
            ind = i;
        }
    }
    return max;                    // return highest value in array
}

/******************************************************************************/

string read_sequence(ifstream& f)
{
    string seq;
    char line[20000];
    while( f.good() )
    {
        f.getline(line,20000);
        if( line[0] == 0 || line[0]=='#' )
            continue;
        for(int i = 0; line[i] != 0; ++i)
        {
            int c = toupper(line[i]);
            if( c != 'A' && c != 'G' && c != 'C' && c != 'T' )
                continue;

            seq.push_back(char(c));
        }
    }
    return seq;
}

/******************************************************************************/
