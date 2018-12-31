%TENSOR_ID_COMPARISON A comparison between different tensor ID algorithms.
%
%   The purpose of this script is to compare different versions of tensor
%   ID algorithms.

%% Settings

N = 10;
I = 1000*ones(1, N);
R = 500;
k = R/2;
l = k + 10;
small_val = 1e-15;
lambda = [ones(1, R/2) small_val*ones(1, R/2)];

%% Generate test tensor

X = generate_dense_tensor(N, I, R, 'lambda_type', 'custom', 'lambda', lambda);

%% Compute gram tensor ID

gram_tic = tic;
Xk_gram = gram_tensor_ID(X, k);
gram_toc = toc(gram_tic);

%% Compute Gaussian tensor ID

gaussian_tic = tic;
Xk_gaussian = gaussian_tensor_ID(X, k, l, 'qr');
gaussian_toc = toc(gaussian_tic);

%% Compute CountSketched tensor ID

CS_tic = tic;
Xk_CS = CS_tensor_ID(X, k, l, 'qr');
CS_toc = toc(CS_tic);

%% Print comparison of results

fprintf('Time for gram tensor ID: %.2f s.\n', gram_toc);
fprintf('Time for Gaussian tensor ID: %.2f s.\n', gaussian_toc);
fprintf('Time for CountSketched tensor ID: %.2f s.\n', CS_toc);
