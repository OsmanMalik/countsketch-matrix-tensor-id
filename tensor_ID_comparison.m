%TENSOR_ID_COMPARISON A comparison between different tensor ID algorithms
%                     [WORK IN PROGRESS]
%
%   This function requires Tensor Toolbox version 2.6 [Ba15].
%
%   The purpose of this script is to compare different versions of tensor
%   ID algorithms:
%       1.  Tensor ID using Gram matrix [Bi15].
%       2.  Gaussian tensor ID [Bi15].
%       3.  CountSketch tensor ID (proposal).
%
%   All of the methods utilize column pivoted QR instead of the strongly
%   rank-revealing QR factorization of [Gu96].
%
% NOTES:
%   1.  This script is a work in progress, and is a bit messy at the
%       moment. I expect to clean it up at a later date.
%
% REFERENCES:
%
%   [Ba15]  B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%           Version 2.6, Available online, February 2015. 
%           URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.
%   
%   [Bi15]  D. J. Biagioni, D. Beylkin, G. Beylkin. Randomized 
%           interpolative decomposition of separated representations. J. 
%           Comput. Phys. 281, pp. 116-134, 2015.
%
%   [Gu96]  M. Gu, and S. C. Eisenstat. Efficient algorithms for computing
%           a strong rank-revealing QR factorization. SIAM J. Sci. Comput.
%           17(1), pp. 848-869, 1996.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

%% Settings

N = 3;
I = 1e+3*ones(1, N);
%I = [532000 17000 3];
%R = 1000;
%k = R/2;
k = 100;
k_mult = 2;
l = k + 20;
mn = 10;
s_norm_tol = 1e-14;
s_norm_verbosity = 1;
s_norm_init = 'mean';
compute_s_norm = true;
lambda = 10.^(-((0:k_mult*k-1)/(k_mult*k))*mn);
%lambda = [1:k 10^(-mn)*ones(1, (k_mult-1)*k)];
%lambda = lambda(randperm(length(lambda)));
%lambda = 10.^(-((0:k_mult*k-1)/(k_mult*k))*mn);
%lambda = fliplr(lambda);
%lambda = 1e-2*randn(k_mult*k,1);
%lambda = 10*ones(k_mult*k,1);
density = 0.01;
compute_gram_id = true;
repeat = k-20;

%% Print program settings

fprintf('Running program with N = %d, I = %d, k = %d\n', N, I(1), k);

%% Generate test tensor

fprintf('Generating tensor... ')
X = generate_tensor(N, I, k_mult*k, 'lambda_type', 'custom', 'lambda', lambda, 'density', density, 'repeat', repeat);
%X = generate_dense_tensor(N, I, R, 'lambda_type', 'exp', 'repeat', 2);
fprintf('Done!\n')

%% Compute gram tensor ID

if compute_gram_id
    fprintf('Computing gram tensor ID... ')
    gram_tic = tic;
    Xk_gram = gram_tensor_ID(X, k);
    gram_toc = toc(gram_tic);
    fprintf('Done!\n')
else
    gram_toc = nan;
end

%% Compute Gaussian tensor ID

fprintf('Computing Gaussian tensor ID... ')
gaussian_tic = tic;
Xk_gaussian = gaussian_tensor_ID(X, k, l, 'qr');
gaussian_toc = toc(gaussian_tic);
fprintf('Done!\n')

%% Compute CountSketched tensor ID

fprintf('Computing CountSketch tensor Id... ')
CS_tic = tic;
Xk_CS = CS_tensor_ID(X, k, l, 'qr');
CS_toc = toc(CS_tic);
fprintf('Done!\n')

%% Print comparison of results

fprintf('Computing errors...\n')
if compute_s_norm
    if compute_gram_id
        gram_s_error = s_norm(X-Xk_gram, s_norm_tol, 'verbosity', s_norm_verbosity, 'init', s_norm_init);
    else
        gram_s_error = nan;
    end
    gaussian_s_error = s_norm(X-Xk_gaussian, s_norm_tol, 'verbosity', s_norm_verbosity, 'init', s_norm_init);
    CS_s_error = s_norm(X-Xk_CS, s_norm_tol, 'verbosity', s_norm_verbosity, 'init', s_norm_init);
else
    gram_s_error = 0;
    gaussian_s_error = 0;
    CS_s_error = 0;
end

fprintf('s-norm for error: %.8e. Time for gram tensor ID: %.2f s.\n', gram_s_error, gram_toc);
fprintf('s-norm for error: %.8e. Time for Gaussian tensor ID: %.2f s.\n', gaussian_s_error, gaussian_toc);
fprintf('s-norm for error: %.8e. Time for CountSketched tensor ID: %.2f s.\n', CS_s_error, CS_toc);
fprintf('Speed gain: %.3f, %.3f\n', gram_toc/CS_toc, gaussian_toc/CS_toc);
