% RUN_EXPERIMENT_3 Run matrix ID experiment for real sparse matrix
%
%   RUN_EXPERIMENT_3 is a script that runs an experiment where a real-world
%   sparse matrix is decomposed using different versions of matrix
%   interpolative decomposition (ID). Depending on where and how the matrix
%   is stored, some of the code below may need to be adjusted.
%   
%   The following versions of matrix ID are used in the comparison:
%       1. Gaussian matrix ID [Ma11]
%       2. SRFT matrix ID [Wo08]
%       3. CountSketch matrix ID (proposal)
%
%   For further details on these methods, see the comment in
%   run_experiment_2.m.
%
%   All of the methods utilize column pivoted QR instead of the strongly
%   rank-revealing QR factorization of [Gu96].
%
% REFERENCES:
%
%   [Gu96]  M. Gu, and S. C. Eisenstat. Efficient algorithms for computing
%           a strong rank-revealing QR factorization. SIAM J. Sci. Comput.
%           17(1), pp. 848-869, 1996.
%
%   [Ma11]  P. G. Martinsson, V. Rokhlin, M. Tygert. A randomized algorithm
%           for the decomposition of matrices. Appl. Comput. Harmon. Anal.
%           30, pp. 47-68, 2011.
%
%   [Wo08]  F. Woolfe, E. Liberty, V. Rokhlin, M. Tygert. A fast randomized
%           algorithm for the approximation of matrices. Appl. Comput.
%           Harmon. Anal. 25, pp. 335-366, 2008.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

%% Settings
% K: The target rank in the decomposition.
% L: This is K + oversampling parameter.
% SRFT_no_splits: Controls how the input matrix is split into pieces in the
%   execution of SRFT_matrix_ID. This helps avoid exceeding memory usage in
%   that algorithm when dealing with sparse matrices.
% norm_type: Set to 'normest' to use Matlab's normest function to compute
%   error, and to 'normrand' to use randomized spectral norm from [Wo08].
% no_rand_norm_vec: The number of random vectors used in the computation of
%   the randomized 2-norm.
% normest_tol: Used to control tolerance of normest, if normest is used.
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
norm_type = 'normest';
no_rand_norm_vec = 10;
normest_tol = 1e-2;
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
    
    % If verbose, print which norm is used
    if verbosity >= 1
        if strcmp(norm_type, 'normrand')
            fprintf('Computing errors using randomized spectral norm...\n');
        elseif strcmp(norm_type, 'normest')
            fprintf('Computing errors using normest...\n');
        else
            error('Invalid norm type. Exiting.');
        end
    end
    
    % If randomized norm is used, compute A*X for random matrix X
    if strcmp(norm_type, 'normrand')    
        X = randn(R, no_rand_norm_vec);
        X = X./sqrt(sum(X.^2,1));
        AX = A*X;
    end
    
    % Compute error for Gaussian ID
    if strcmp(norm_type, 'normrand')
        error_GA = max(sqrt(sum(abs(AX - A(:,J_GA)*(P_GA*X)).^2, 1)));
    elseif strcmp(norm_type, 'normest')
        error_GA = normest(A - A(:,J_GA)*P_GA, normest_tol);
    end
    if verbosity >= 1
        fprintf('Gaussian matrix ID error: %.10e. Time: %.2f s.\n', error_GA, toc_gaussian);    
    end
    
    % Compute error for SRFT ID
    if strcmp(norm_type, 'normrand')
        error_SRFT = max(sqrt(sum(abs(AX - A(:,J_SRFT)*(P_SRFT*X)).^2, 1)));
    elseif strcmp(norm_type, 'normest')
        error_SRFT = normest(A - A(:,J_SRFT)*P_SRFT, normest_tol);
    end
    if verbosity >= 1
        fprintf('SRFT matrix ID error: %.10e. Time: %.2f s.\n', error_SRFT, toc_SRFT);
    end

    % Compute error for CS ID
    if strcmp(norm_type, 'normrand')
        error_CS = max(sqrt(sum(abs(AX - A(:,J_CS)*(P_CS*X)).^2, 1)));
    elseif strcmp(norm_type, 'normest')
        error_CS = normest(A - A(:,J_CS)*P_CS, normest_tol);
    end
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