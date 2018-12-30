%MATRIX_ID_COMPARISON   A comparison between three different matrix ID
%                       algorithms.
%
%   The purpose of this script is to compare three versions of matrix ID:
%       1. Deterministic, according to [1];
%       2. Randomized using a Gaussian random matrix, according to [3];
%       3. Randomized using subsampled randomized Fourier transform (SRFT),
%           according to [4];
%       4. Randomized using CountSketch (our proposal).
%
%   The deterministic matrix ID, which is used for 1, and as a component in
%   2, 3 and 4, uses the method described in [1], which incorporates the
%   strong rank-revealing QR (SRRQR) factorization of [2]. We use the 
%   implementation [5] of SRRQR based matrix ID.
%
% NOTES:
%   1.  The practical algorithm suggested in Theorem 2 of [1] corresponds to
%       choosing f = sqrt(n) in [2], where n is the number of columns of A,
%       the matrix to be decomposed.
%   2.  I think the practical algorithm for computing matrix ID of the
%       product G*A, where G is Gaussian, in [3] uses f = 2 in [2].
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
%   [4] F. Woolfe, E. Liberty, V. Rokhlin, M. Tygert. A fast randomized
%       algorithm for the approximation of matrices. Appl. Comput. Harmon.
%       Anal. 25, pp. 335-366, 2008.
%   
%   [5] X. Xing. Interpolative Decomposition based on Strong RRQR. MATLAB
%       Central File Exchange. Retrieved November 23, 2018.

%% Settings 

m = 2000;
n = 2000;
k = 200; % Numerical rank of matrix
l = k+8; % Oversampled l > k
noise_level = 1e-5;
f = 2;
matrix_type = 3;
mn = 8;

%% Generate test matrix

if matrix_type == 1
    [Q1, ~] = qr(randn(m));
    [Q2, ~] = qr(randn(n));
    mn = min(m, n);
    S = diag([mn:-1:mn-k+1 noise_level*randn(1,mn-k)]);
    if m > n 
        S = [S; zeros(m-n, n)];
    elseif m < n
        S = [S zeros(m, n - m)];
    end
    A = Q1*S*Q2.';
elseif matrix_type == 2
    A = full(sprand(m, n, .01));
elseif matrix_type == 3
    A = generate_dense_matrix(n, k, mn);
elseif matrix_type == 4
    A = sprand(m, n, .01);
end

%% Compute deterministic matrix ID
% This is matrix ID according to [1], implemented in [5]

time_srrqr_id_tic = tic;
    [P1, J1] = ID(A', 'rank', k, f);
    P1 = P1';
time_srrqr_id = toc(time_srrqr_id_tic);

%% Compute deterministic matrix ID using Matlab QR
% This is matrix ID according to [1], but using Matlab QR

time_qr_id_tic = tic;
    [Q1b, R1b, e1b] = qr(A, 0);
    T1b = R1b(1:k, 1:k) \ R1b(1:k, k+1:end);
    P1b = [eye(k) T1b];
    clear pvec1b;
    pvec1b(e1b) = 1:length(e1b);
    P1b = P1b(:, pvec1b);
    J1b = e1b(1:k)';
time_qr_id = toc(time_qr_id_tic);

%% Compute randomized Gaussian matrix ID
% This is matrix ID according to [3]

time_srrqr_gaussian_id_tic = tic;
    G = randn(l, m);
    R = G*A;
    [P2, J2] = ID(R', 'rank', k, f);
    P2 = P2';
time_srrqr_gaussian_id = toc(time_srrqr_gaussian_id_tic);

%% Compute randomized SRFT matrix ID
% This is matrix ID according to [4]. Note that [5] cannot handle complex
% numbers (indeed, [2] doesn't even mention complex numbers), so I just do
% a matrix ID using QR with column pivoting here.

Atemp = A;
time_qr_srft_id_tic = tic;
    [P3, J3] = SRFT_matrix_ID(A, k, l);
time_qr_srft_id = toc(time_qr_srft_id_tic);

%% Compute randomized CountSketch matrix ID (our proposal)

time_srrqr_cs_id_tic = tic;
    [P4, J4] = CS_matrix_ID(A, k, l, 'srrqr');
time_srrqr_cs_id = toc(time_srrqr_cs_id_tic);
    
%% Compute randomized CountSketch matrix ID using Matlab QR (our proposal)

time_qr_cs_id_tic = tic;
    [P5, J5] = CS_matrix_ID(A, k, l, 'qr');
time_qr_cs_id = toc(time_qr_cs_id_tic);

%% Do SVD for comparison

time_svd_tic = tic;
    [U5, S5, V5] = svd(A);
time_svd = toc(time_svd_tic);

%% Print comparison of results

fprintf('Error for rank %d SVD: %.10e. Time %.2f s.\n', k, S5(k+1, k+1), time_svd)
fprintf('Error for rank %d matrix ID [SRRQR]: %.10e. Time %.2f s.\n', k, norm(A - A(:,J1)*P1), time_srrqr_id)
fprintf('Error for rank %d matrix ID [Matlab QR]: %.10e. Time %.2f s.\n', k, norm(A - A(:,J1b)*P1b), time_qr_id)
fprintf('Error for rank %d matrix random Gaussian ID [SRRQR]: %.10e. Time %.2f s.\n', k, norm(A - A(:,J2)*P2), time_srrqr_gaussian_id)
fprintf('Error for rank %d matrix random SRFT ID [Matlab QR]: %.10e. Time %.2f s.\n', k, norm(A - A(:,J3)*P3), time_qr_srft_id)
fprintf('Error for rank %d matrix random CountSketch ID [SRRQR]: %.10e. Time %.2f s.\n', k, norm(A - A(:,J4)*P4), time_srrqr_cs_id)
fprintf('Error for rank %d matrix random CountSketch ID [Matlab QR]: %.10e. Time %.2f s.\n\n', k, norm(A - A(:,J5)*P5), time_qr_cs_id)