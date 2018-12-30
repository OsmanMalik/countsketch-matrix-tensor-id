function SA = TensorSketch(A, h, s, J)
%TENSORSKETCH Computes the TensorSketch of a matrix Khatri-Rao product
%
%   SA = TENSORSKETCH(A, h, s, J) returns the TensorSketch of the
%   Khatri-Rao product of the matrices in A computed using the hash 
%   functions in h and s, using a target sketch dimension J. A should be a 
%   (row or column) cell containing the matrices, h should be a (row or 
%   column) cell containing hash functions, and s should be a (row or 
%   column) cell containing sign hash functions.

%% Include relevant files

addpath(genpath('help_functions'));

%% Computations

N = length(A);
R = size(A{1}, 2);
Acs = cell(size(A)); % To store CountSketch of each matrix in A. FFT and transpose are also applied.
P = ones(J, R);

for n = 1:N
    Acs{n} = fft(countSketch(A{n}.', int64(h{n}), J, s{n}, 1).');
    P = P.*Acs{n};
end

SA = ifft(P);

end