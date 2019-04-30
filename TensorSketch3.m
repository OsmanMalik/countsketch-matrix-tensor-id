function SA = TensorSketch3(A, J)
%TENSORSKETCH Computes the TensorSketch of a matrix Khatri-Rao product
%
%   SA = TENSORSKETCH(A, h, s, J) returns the TensorSketch of the
%   Khatri-Rao product of the matrices in A computed using the hash 
%   functions in h and s, using a target sketch dimension J. A should be a 
%   (row or column) cell containing the matrices, h should be a (row or 
%   column) cell containing hash functions, and s should be a (row or 
%   column) cell containing sign hash functions. The matrices in A can be
%   either dense or sparse: The appropriate countSketch function will be
%   used in each case.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

%% Include relevant files

%addpath(genpath('help_functions'));

%% Computations

N = length(A);
R = size(A{1}, 2);
Acs = cell(size(A)); % To store CountSketch of each matrix in A. FFT and transpose are also applied.
P = ones(J, R);

% Define permutation map and multinomial count and sign function
perm_map = cell(N, 1);
mult_count = cell(N, 1);
mult_count_cumsum = cell(N, 1);
s = cell(N, 1);
for n = 1:N
    perm_map{n} = randperm(size(A{n}, 1));
    mult_count{n} = mnrnd(size(A{n}, 1), ones(1, J) / J);
    mult_count_cumsum{n} = cumsum(mult_count{n});
    s{n} = randi(2, size(A{n}, 1), 1)*2-3;
end

for n = 1:N
    Acs{n} = fft(countSketch_sparse_2(A{n}.', int64(perm_map{n}), int64(mult_count{n}), int64(mult_count_cumsum{n}), J, s{n}).');
    P = P.*Acs{n};
end

SA = ifft(P);

end