N = 5;
%Is = [10*1e+3 25*1e+3 50*1e+3 100*1e+3 250*1e+3 500*1e+3 1e+6];
%Is = [1e+3 2e+3];
Is = [1e+3 2.5e+3 5e+3 1e+4 2.5e+4 5e+4 1e+5 2.5e+5];
K = 1e+3;
L = K + 10;
mn = 8;
no_trials = 10;
K_mult = 10;
lambda = 10.^(-((0:K_mult*K-1)/(K_mult*K))*mn);
fac_mat_dens = .01;
s_norm_tol = 1e-12;
maxit = 1000;
results_matlab_file = 'matlab_output_exp_4';
verbosity = 1;
cnt = 1;

I = Is(6)*ones(1, N);

% Generate sparse tensor X
if verbosity >= 1
    fprintf('Generating sparse tensor... ');
end
X = generate_tensor(N, I, K_mult*K, 'lambda_type', 'custom', 'lambda', lambda, 'density', fac_mat_dens);
if verbosity >= 1
    fprintf('Done!\n')
end

% Compute Gaussian tensor ID
if verbosity >= 1
    fprintf('Running Gaussian tensor ID... ');
end
gaussian_tic = tic;
Xk_gaussian = gaussian_tensor_ID(X, K, L, 'qr');
gaussian_toc = toc(gaussian_tic);
if verbosity >= 1
    fprintf('Done!\n');
end

% Compute CountSketch tensor ID
if verbosity >= 1
    fprintf('Running CountSketch tensor ID... ');
end
CS_tic = tic;
Xk_CS = CS_tensor_ID(X, K, L, 'qr');
CS_toc = toc(CS_tic);
if verbosity >= 1
    fprintf('Done!\n');
end 