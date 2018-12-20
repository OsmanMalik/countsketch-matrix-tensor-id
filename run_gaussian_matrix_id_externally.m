function [P, J, time] = run_gaussian_matrix_id_externally(k, l, n)
%RUN_GAUSSIAN_MATRIX_ID_EXTERNALLY Call external Gaussian matrix ID C code
%
%   [P, J, time] = RUN_GAUSSIAN_MATRIX_ID_EXTERNALLY(k, l, n) 

% Execute C code for Gaussian matrix ID
system(['./run_gaussian_matrix_id ', num2str(k), ' ', num2str(l)]);

% Load outputs from C file
I = csvread('output/run_gaussian_matrix_id_I.txt');
T = csvread('output/run_gaussian_matrix_id_T.txt');
time = csvread('output/run_gaussian_matrix_id_time.txt');

% Compute function outputs
J = I(1:k)+1;
Tmat = reshape(T, k, n-k); 
Ji(I+1) = 1:length(I);
P = [eye(k) Tmat];
P = P(:, Ji);

end