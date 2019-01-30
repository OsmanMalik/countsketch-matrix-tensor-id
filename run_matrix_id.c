/* RUN_MATRIX_ID
 * Read matrix from file and apply matrix ID to it 
 *
 * This program reads a matrix stored in a file and then applies 
 * matrix ID from RSVDPACK [Vo16] to it. It then stores the
 * resulting decomposition in the form of a vector J and a matrix
 * T. More specifically, we have 
 * 		A \approx A(:, J(1:k))*P,
 * where k is the target rank and P = [Ik T']*Perm, and where Ik is the
 * k by k identity matrix, Perm is the n by n permutation matrix
 * Perm = In(:, J), and In is the n by n identity matrix; see
 * Section 2.4 in [Vo16] for further details.
 *
 * The program takes one optional input:
 * 		./run_matrix_id k
 * where k is the target rank.
 *
 * Note that this code is not intended to be called directly from Matlab.
 * Instead, we recommend using the Matlab wrapper 
 * run_matrix_id_externally.m.
 *
 * REFERENCES:
 *
 * [Vo16]	S. Voronin, and P. G. Martinsson. RSVDPACK: An implementation 
 *			of randomized algorithms for computing the singular value,
 * 			interpolative, and CUR decompositions of matrices on multi-core
 *		 	and GPU architectures. arXiv:1502.05366v3 [math.NA], 2016.
 * */

/* Author:   Osman Asif Malik
 * Email:    osman.malik@colorado.edu
 * Date:     January 29, 2019
 * */
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rank_revealing_algorithms_intel_mkl.h"

static double TimeSpecToSeconds(struct timespec *ts) {
	return (double)ts->tv_sec + (double)ts->tv_nsec / 1000000000.0;	
}

int main(int argc, char **argv) {

/* Declare variables */
int default_k, frank, k, m, n, i, j;
double TOL, elapsedSeconds;
mat *M, *T;
vec *I;
struct timespec start_time, end_time;
char *M_file, *output_frank, *output_I, *output_T, *output_time;
FILE *fp;

/* Settings */
M_file = "data/A_mat.bin";
output_frank = "output/run_matrix_id_frank.txt";
output_I = "output/run_matrix_id_I.txt";
output_T = "output/run_matrix_id_T.txt";
output_time = "output/run_matrix_id_time.txt";
default_k = 1;

/* Handle optional input arguments */
if(argc < 2) {
	k = default_k;
	printf("No rank provided. Using default of %d\n", k);
} else {
	k = atoi(argv[1]);
	printf("Using provided rank %d\n", k);
}
TOL = 0.0;

/* Load matrix from file */
printf("Loading matrix from %s...\n", M_file);
M = matrix_load_from_binary_file(M_file);
m = M->nrows;
n = M->ncols;
printf("Finished loading matrix of size %d by %d!\n", m, n);

/* Run matrix ID  */
clock_gettime(CLOCK_MONOTONIC, &start_time);
id_decomp_fixed_rank_or_prec(M, k, TOL, &frank, &I, &T);
clock_gettime(CLOCK_MONOTONIC, &end_time);
elapsedSeconds = TimeSpecToSeconds(&end_time) - TimeSpecToSeconds(&start_time);
printf("Finished running ID decomposition in %.4f s. Resulting rank is %d\n", elapsedSeconds, frank);

/* Save frank to file */
printf("Output matrix frank is %d. Saving frank... ", frank);
fp = fopen(output_frank, "w+");
fprintf(fp, "%d\n", frank);
fclose(fp);
printf("Done!\n");

/* Save I to file */
printf("Saving all %d rows of I... ", I->nrows);
fp = fopen(output_I, "w+");
for(i = 0; i < I->nrows; ++i) {
	fprintf(fp, "%.0f\n", I->d[i]);
}
fclose(fp);
printf("Done!\n");

/* Saving T to file */
printf("Saving T to file... ");
fp = fopen(output_T, "w+");
printf("T has %d rows and %d columns. Now printing all of T in column major order.\n", T->nrows, T->ncols);
for(i = 0; i < (T->nrows)*(T->ncols); ++i) {
	fprintf(fp, "%.20e\n", T->d[i]);
}
fclose(fp);
printf("Done!\n");

/* Saving time to file */
printf("Saving time to file... ");
fp = fopen(output_time, "w+");
fprintf(fp, "%.20e\n", elapsedSeconds);
fclose(fp);
printf("Done!\n");

/* Check decomposition */
use_id_decomp_for_approximation(M, T, I, k);

return 0;

}
