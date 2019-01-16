% RUN_EXPERIMENT_3 Run matrix ID experiment for real sparse matrix

%% Settings
% K: The target rank in the decomposition.
% L: This is K + oversampling parameter.
% SRFT_no_splits: Controls how the input matrix is split into pieces in the
%   execution of SRFT_matrix_ID. This helps avoid exceeding memory usage in
%   that algorithm when dealing with sparse matrices.
% no_rand_norm_vec: The number of random vectors used in the computation of
%   the randomized 2-norm.
% no_trials: The number of times each experiment is repeated.
% filename: This is the location/name of the file where the input data is
%   stored.
% results_matlab_file: This is the name of the mat file where the results
%   are stored 
% verbosity: Controls the verbosity (0 = least verbose, 1 = intermediate 
%   verbosity, 2 = max verbosity)

K = 1442;
L = K + 10;
SRFT_no_splits = 16;
no_rand_norm_vec = 10;
no_trials = 10;
filename = '../data/specular.mat';
results_matlab_file = 'matlab_output_exp_3';
verbosity = 1;

%% Main loop

load(filename);
A = Problem.A;
[I, R] = size(A);

% Create mat and set up matfile for saving results computed in Matlab
save_mat = matfile(results_matlab_file, 'Writable', true);
save_mat.trial = zeros(1, no_trials);
save_mat.time = zeros(3, no_trials);
save_mat.error = zeros(3, no_trials);

fprintf('Starting Experiment 3...\n');

for tr = 1:no_trials

    if verbosity >= 1
        fprintf('\nStarting trial %d\n', tr);
    end
    
    % Compute Gaussian matrix ID
    if verbosity >= 1
        fprintf('Running Gaussian matrix ID... ');
    end
    tic_gaussian = tic;
    [P_GA, J_GA] = Gaussian_matrix_ID(A, K, L, 'qr');
    toc_gaussian = toc(tic_gaussian);
    if verbosity >= 1
        fprintf('Done!\n');
    end
    
    % Compute SRFT matrix ID
    if verbosity >= 1
        fprintf('Running SRFT matrix ID... ');
    end
    tic_SRFT = tic;
    [P_SRFT, J_SRFT] = SRFT_matrix_ID(A, K, L, 'splits', SRFT_no_splits);
    toc_SRFT = toc(tic_SRFT);
    if verbosity >= 1
        fprintf('Done!\n');
    end
    
    % Compute CountSketch matrix ID
    if verbosity >= 1
        fprintf('Running CountSketched matrix ID... ');
    end
    tic_CS = tic;
    [P_CS, J_CS] = CS_matrix_ID(A, K, L, 'qr');
    toc_CS = toc(tic_CS);
    if verbosity >= 1
        fprintf('Done!\n');
    end
    
    % Compute errors
    if verbosity >= 1
        fprintf('Computing errors...\n');
    end
    X = randn(R, no_rand_norm_vec);
    X = X./sqrt(sum(X.^2,1));
    AX = A*X;
    
    error_GA = max(sqrt(sum(abs(AX - A(:,J_GA)*(P_GA*X)).^2, 1)));
    if verbosity >= 1
        fprintf('Gaussian matrix ID error: %.10e. Time: %.2f s.\n', error_GA, toc_gaussian);    
    end
    
    error_SRFT = max(sqrt(sum(abs(AX - A(:,J_SRFT)*(P_SRFT*X)).^2, 1)));
    if verbosity >= 1
        fprintf('SRFT matrix ID error: %.10e. Time: %.2f s.\n', error_SRFT, toc_SRFT);
    end

    error_CS = max(sqrt(sum(abs(AX - A(:,J_CS)*(P_CS*X)).^2, 1)));
    if verbosity >= 1
        fprintf('CountSketch matrix ID error: %.10e. Time: %.2f s.\n', error_CS, toc_CS);
    end
    
    % Saving errors and times
    if verbosity >= 1
        fprintf('Done computing errors. \nSaving to file... ')
    end
    save_mat.trial(1, tr)  = tr;   
    save_mat.error(1, tr)  = error_GA; 
    save_mat.error(2, tr)  = error_SRFT;
    save_mat.error(3, tr)  = error_CS;
    save_mat.time(1, tr)   = toc_gaussian;
    save_mat.time(2, tr)   = toc_SRFT;
    save_mat.time(3, tr)   = toc_CS;
    if verbosity >= 1
        fprintf('Done!\n')
    end
end