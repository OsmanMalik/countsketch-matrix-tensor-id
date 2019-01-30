function [P, J, time] = run_matrix_id_externally(k, n)
% RUN_MATRIX_ID_EXTERNALLY Call external matrix ID C code
%
%   [P, J, time] = RUN_MATRIX_ID_EXTERNALLY(k, n) calls the external
%   executable run_matrix_id, which is compiled from C code. That
%   executable utilizes the functionality provided in RSVDPACK [Vo16] for
%   computing the matrix ID of a matrix saved on disk. The inputs k and
%   n are the target rank and number of columns of the decomposed matrix,
%   respectively. The resulting decomposition is of the form
%       A \approx A(:, J)*P.
%   The output variable time contains the time it took for the C code to
%   produce the decomposition.
%
% REFERENCES:
%
%   [Vo16]  S. Voronin, and P. G. Martinsson. RSVDPACK: An implementation
%           of randomized algorithms for computing the singular value, 
%           interpolative, and CUR decompositions of matrices on multi-core
%           and GPU architectures. arXiv:1502.05366v3 [math.NA], 2016.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

% Execute C code for matrix ID
system(['./run_matrix_id ', num2str(k)]);

% Load outputs from C file
frank = csvread('output/run_matrix_id_frank.txt');
I = csvread('output/run_matrix_id_I.txt');
T = csvread('output/run_matrix_id_T.txt');
time = csvread('output/run_matrix_id_time.txt');

% Check outputs
if frank ~= k
    warning(['frank not equal to rank: frank = ', num2str(frank), ', rank = ', num2str(k)]);
end

% Compute function outputs
J = I(1:k)+1;
Tmat = reshape(T, k, n-k); 
Ji(I+1) = 1:length(I);
P = [eye(k) Tmat];
P = P(:, Ji);

end

