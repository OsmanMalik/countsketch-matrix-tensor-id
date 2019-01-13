%TENSOR_ID_COMPARISON A comparison between different tensor ID algorithms.
%
%   The purpose of this script is to compare different versions of tensor
%   ID algorithms.

%% Settings

N = 3;
%I = 30*1000*ones(1, N);
I = [532000 17000 3];
R = 1/2*2000;
k = 1/2*1000;
l = k + 10;
small_val = 1e-4;
s_norm_tol = 1e-14;
s_norm_verbosity = 1;
s_norm_init = 'mean';
compute_s_norm = true;
lambda = [ones(1, R/2) small_val*ones(1, R/2)];
density = 0.01;

%% Print program settings

fprintf('Running program with N = %d, I = %d, R = %d, k = %d\n', N, I(1), R, k);

%% Generate test tensor

X = generate_tensor(N, I, R, 'lambda_type', 'custom', 'lambda', lambda, 'density', density);
%X = generate_dense_tensor(N, I, R, 'lambda_type', 'exp', 'repeat', 2);

%% Compute gram tensor ID

%{
gram_tic = tic;
Xk_gram = gram_tensor_ID(X, k);
gram_toc = toc(gram_tic);
%}
gram_toc = 0;

%% Compute Gaussian tensor ID

gaussian_tic = tic;
Xk_gaussian = gaussian_tensor_ID(X, k, l, 'qr');
gaussian_toc = toc(gaussian_tic);

%% Compute CountSketched tensor ID

CS_tic = tic;
Xk_CS = CS_tensor_ID(X, k, l, 'qr');
CS_toc = toc(CS_tic);

%% Print comparison of results

if compute_s_norm
    %gram_s_error = s_norm(X-Xk_gram, s_norm_tol, 'verbosity', s_norm_verbosity, 'init', s_norm_init);
    gram_s_error = 0;
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
fprintf('Speed gain: %.3f\n', gaussian_toc/CS_toc);
