#!/bin/bash

RSVDPACK_INC="/projects/osma9213/own_build/LowRankMatrixDecompositionCodes/multi_core_mkl_code/"


icc -O2 -mkl -qopenmp -I "$RSVDPACK_INC" run_matrix_id.c "$RSVDPACK_INC"rank_revealing_algorithms_intel_mkl.c "$RSVDPACK_INC"matrix_vector_functions_intel_mkl.c -o run_matrix_id

icc -O2 -mkl -qopenmp -I "$RSVDPACK_INC" run_gaussian_matrix_id.c "$RSVDPACK_INC"rank_revealing_algorithms_intel_mkl.c "$RSVDPACK_INC"matrix_vector_functions_intel_mkl.c -o run_gaussian_matrix_id
