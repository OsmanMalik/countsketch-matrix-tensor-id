function [P, J] = Gaussian_matrix_ID(A, k, l, QR_type)
% GAUSSIAN_MATRIX_ID Computes Gaussian matrix ID
%
%   [P, J] = Gaussian_matrix_ID(A, k, l, QR_type) computes the rank-k
%   matrix ID
%       A \approx A(:, J)*P,
%   where A(:, J) contains k columns of A. The argument l > k is k plus the
%   oversampling parameter. QR_type controls the type of QR factorization
%   used: Set it to 'srrqr' to use the strong rank-revealing QR
%   factorization of [Gu96] (uses the implementation [Xi18]), or set it to
%   'qr' to use Matlab's built-in QR function. Before applying the matrix
%   ID, we sketch the matrix using a Gaussian matrix. This code is based on
%   Algorithm 6 in the appendix of [Vo15].
%
%   REFERENCES:
%
%   [Gu96]  M. Gu, and S. C. Eisenstat. Efficient algorithms for computing 
%           a strong rank-revealing QR factorization. SIAM J. Sci. Comput. 
%           17(1), pp. 848-869, 1996.
%
%   [Vo15]  S. Voronin, and P.-G. Martinsson. RSVDPACK: An implementation
%           of randomized algorithms for computing the singular value,
%           interpolative, and CUR decompositions of matrices on multi-core
%           and GPU architectures. arXiv:1502.05366 [cs, math].
%
%   [Xi18]  X. Xing. Interpolative Decomposition based on Strong RRQR. 
%           MATLAB Central File Exchange. Retrieved November 23, 2018.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

[m, ~] = size(A);

% Compute Gaussian sketch
G = randn(l, m);
Y = full(G*A);

% Proceed by using appropriate QR factorization
if strcmp(QR_type, 'srrqr')
    f = 2;
    [P, J] = ID(Y', 'rank', k, f);
    P = P';
elseif strcmp(QR_type, 'qr')
    [~, R, e] = qr(Y, 0);
    T = R(1:k, 1:k) \ R(1:k, k+1:end);
    P = [eye(k) T];
    pvec(e) = 1:length(e);
    P = P(:, pvec);
    J = e(1:k)';
end

end