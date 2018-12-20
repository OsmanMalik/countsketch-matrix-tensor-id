% RUN_EXPERIMENT_1 Run matrix ID experiment for dense or sparse matrices
%
%   RUN_EXPERIMENT_1 is a script that runs an experiment where dense 
%   matrices of known rank and with known singular values, OR a random
%   sparse matrices, are decomposed using different versions of matrix 
%   interpolative decomposition (ID). The following versions of matrix ID 
%   are used:
%       1. Matrix ID [1]
%       2. Gaussian matrix ID [3]
%       3. SRFT matrix ID [5]
%       4. CountSketch matrix ID utilizing SRRQR (proposal)
%       5. CountSketch matrix ID utilizing Matlab QR (proposal)
%   For 1 and 2, we utilize an existing implementation in the software
%   package RSVDPACK [4]. We run do this by saving the matrices to a text
%   file, and then reading them into the functions provided in RSVDPACK for 
%   matrix ID. 4 and 5 are two different versions of our proposal. The two
%   methods are identical, with the only difference that 4 utilizes the
%   strong rank-revealing QR factorization of [2], more precisely the
%   implementation in [6], whereas 5 utilizes Matlab's QR function. We also
%   compute the SVD for each given rank and use it as a comparison.
%
% REFERENCES:
%   [1] H. Cheng, Z. Gimbutas, P. G. Martinsson, and V. Rokhlin. On the
%       compression of low rank matrices. SIAM J. Sci. Comput. 26(4), pp.
%       1389-1404, 2005.
%
%   [2] M. Gu, and S. C. Eisenstat. Efficient algorithms for computing a
%       strong rank-revealing QR factorization. SIAM J. Sci. Comput. 17(1),
%       pp. 848-869, 1996.
%
%   [3] P. G. Martinsson, V. Rokhlin, M. Tygert. A randomized algorithm for
%       the decomposition of matrices. Appl. Comput. Harmon. Anal. 30, pp.
%       47-68, 2011.
%
%   [4] S. Voronin, and P. G. Martinsson. RSVDPACK: An implementation of 
%       randomized algorithms for computing the singular value, 
%       interpolative, and CUR decompositions of matrices on multi-core and
%       GPU architectures. arXiv:1502.05366v3 [math.NA], 2016.
%
%   [5] F. Woolfe, E. Liberty, V. Rokhlin, M. Tygert. A fast randomized
%       algorithm for the approximation of matrices. Appl. Comput. Harmon.
%       Anal. 25, pp. 335-366, 2008.
%   
%   [6] X. Xing. Interpolative Decomposition based on Strong RRQR. MATLAB
%       Central File Exchange. Retrieved November 23, 2018.

%% Settings
% n: The matrix will be of size n by n
% oversampling: The oversampling parameter to be used
% no_trials: The number of times each experiment is repeated
% bin_file: This is the file where the matrix A is stored for input into C
%   functions
% results_matlab_file: This is the name of the mat file where the results
%   are stored
% verbosity: Controls the verbosity (0 = least verbose, 1 = intermediate 
%   verbosity, 2 = max verbosity)
% matrix_type: Determines which kind of random matrix A is. Can be either
%   'dense' or 'sparse'
% ranks: Vector containing all ranks to use in experiment
% min_sv: The smallest nonzero singular vector will be of size 10^(-min_sv)
% sparse_matrix_density: Used to control the density of A when it is a
%   random sparse matrix

n = 10000;
oversampling = 10;
no_trials = 10;
bin_file = 'data/A_mat.bin';
results_matlab_file = 'matlab_output';
verbosity = 1;
matrix_type = 'dense';
ranks = [50 100 200 500 1000 2000 5000]; % Only used if matrix_type = 'dense'
min_sv = 8; % Only used if matrix_type = 'dense'
sparse_matrix_density = 0.05; % Only used if matrix_type = 'sparse'

%% Main loop

% Create mat and set up matfile for saving results computed in Matlab
save_mat = matfile(results_matlab_file, 'Writable', true);
save_mat.rank = zeros(1, length(ranks)*no_trials);
save_mat.trial = zeros(1, length(ranks)*no_trials);
save_mat.time = zeros(6, length(ranks)*no_trials);
save_mat.error = zeros(6, length(ranks)*no_trials);
cnt = 1;

fprintf('Starting Experiment 1...\n\n');

for tr = 1:no_trials
    for j = 1:length(ranks)
        k = ranks(j);
        l = k + oversampling;

        if verbosity >= 1
            fprintf('Starting trial %d for rank %d\n', tr, k);
        end

        % Generate matrix A
        if verbosity >= 1
            fprintf('Generating %s %d by %d matrix... ', matrix_type, n, n);
        end
        if strcmp(matrix_type, 'dense')
            A = generate_dense_matrix(n, k, min_sv);
        elseif strcmp(matrix_type, 'sparse')
            A = sprand(n, n, sparse_matrix_density);
        else
            error('Invalid matrix_type');
        end
        if verbosity >= 1
            fprintf('Done!\n');
        end

        % Save matrix A to file
        if verbosity >= 1
            fprintf('Saving matrix to file...\n');
        end
        save_matrix_to_file(A, bin_file, verbosity);
        if verbosity >= 1
            fprintf('Finished writing to file!\n');
        end

        % Compute matrix ID (RSVDPACK)
        if verbosity >= 1
            fprintf('Running RSVDPACK matrix ID... ');
        end
        [P_STD, J_STD, time_STD] = run_matrix_id_externally(k, n);
        if verbosity >= 1
            fprintf('Done!\n');
        end

        % Compute Gaussian matrix ID (RSVDPACK)
        if verbosity >= 1
            fprintf('Running RSVDPACK Gaussian matrix ID... ');
        end
        [P_GA, J_GA, time_GA] = run_gaussian_matrix_id_externally(k, oversampling, n);
        if verbosity >= 1
            fprintf('Done!\n');
        end

        % Compute SRFT matrix ID (Matlab qr)
        if verbosity >= 1
            fprintf('Running SRFT matrix ID... ');
        end
        tic_SRFT = tic;
        [P_SRFT, J_SRFT] = SRFT_matrix_ID(A, k, l);
        toc_SRFT = toc(tic_SRFT);
        if verbosity >= 1
            fprintf('Done!\n');
        end

        % Compute CountSketch matrix ID (utilizing SRRQR of [2])
        if verbosity >= 1
            fprintf('Running CountSketched matrix ID utilizing SRRQR... ');
        end
        tic_CS_SRRQR = tic;
        [P_CS_SRRQR, J_CS_SRRQR] = CS_matrix_ID(A, k, l, 'srrqr');
        toc_CS_SRRQR = toc(tic_CS_SRRQR);
        if verbosity >= 1
            fprintf('Done!\n');
        end

        % Compute CountSketch matrix ID (utilizing Matlab qr)
        if verbosity >= 1
            fprintf('Running CountSketched matrix ID utilizing Matlab QR... ');
        end
        tic_SC_QR = tic;
        [P_CS_QR, J_CS_QR] = CS_matrix_ID(A, k, l, 'qr');
        toc_SC_QR = toc(tic_SC_QR);
        if verbosity >= 1
            fprintf('Done!\n');
        end

        % Compute SVD for comparison
        if verbosity >= 1
            fprintf('Computing SVD... ');
        end
        tic_SVD = tic;
        [U, S, V] = svd(A);
        toc_SVD = toc(tic_SVD);
        if verbosity >= 1
            fprintf('Done!\n');
        end

        % Computing errors
        if verbosity >= 1
            fprintf('Computing errors...\n');
        end
        error_STD = norm(A - A(:, J_STD)*P_STD);
        if verbosity >= 1
            fprintf('RSVDPACK matrix ID error: %.10e. Time: %.2f s.\n', error_STD, time_STD);
        end
        error_GA = norm(A - A(:, J_GA)*P_GA);
        if verbosity >= 1
            fprintf('RSVDPACK Gaussian matrix ID error: %.10e. Time: %.2f s.\n', error_GA, time_GA);    
        end
        error_SRFT = norm(A - A(:, J_SRFT)*P_SRFT);
        if verbosity >= 1
            fprintf('SRFT matrix ID error: %.10e. Time: %.2f s.\n', error_SRFT, toc_SRFT);
        end
        error_CS_SRRQR = norm(A - A(:, J_CS_SRRQR)*P_CS_SRRQR);
        if verbosity >= 1
            fprintf('CountSketch [SRRQR] matrix ID error: %.10e. Time: %.2f s.\n', error_CS_SRRQR, toc_CS_SRRQR);
        end
        error_CS_QR = norm(A - A(:, J_CS_QR)*P_CS_QR);
        if verbosity >= 1
            fprintf('CountSketch [QR] matrix ID error: %.10e. Time: %.2f s.\n', error_CS_QR, toc_SC_QR);
        end
        error_SVD = norm(A - U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)');
        if verbosity >= 1
            fprintf('SVD error: %.10e. Time: %.2f s.\n', error_SVD, toc_SVD);
        end

        % Saving errors and times
        if verbosity >= 1
            fprintf('Done computing errors. \nSaving to file... ')
        end
        save_mat.rank(1, cnt) = k;
        save_mat.trial(1, cnt) = tr;
        save_mat.error(1, cnt) = error_SVD; 
        save_mat.error(2, cnt) = error_STD;
        save_mat.error(3, cnt) = error_GA;
        save_mat.error(4, cnt) = error_SRFT; 
        save_mat.error(5, cnt) = error_CS_SRRQR; 
        save_mat.error(6, cnt) = error_CS_QR; 
        save_mat.time(1, cnt) = toc_SVD;
        save_mat.time(2, cnt) = time_STD;
        save_mat.time(3, cnt) = time_GA;
        save_mat.time(4, cnt) = toc_SRFT;
        save_mat.time(5, cnt) = toc_CS_SRRQR;
        save_mat.time(6, cnt) = toc_SC_QR;
        cnt = cnt + 1;
        if verbosity >= 1
            fprintf('Done!\n')
        end
    end
end