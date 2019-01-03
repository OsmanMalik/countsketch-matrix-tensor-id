% Compile without blas
%mex countSketch.c

% Compile with blas
mex countSketch.c -DUSE_BLAS