function [P, J] = SRFT_matrix_ID(A, k, l)
%SRFT_MATRIX_ID Computes subsampled randomized FFT matrix ID
%
%   [P, J] = SRFT_MATRIX_ID(A, k, l) computes the rank-k matrix ID
%       A \approx A(:, J)*P,
%   where A(:, J) contains k columns of A. The argument l > k is k plus the
%   oversampling parameter. Before applying the matrix ID, we sketch the 
%   matrix using a subsampled randomized FFT as proposed in [1].
%
% REFERENCES:
%   [1] F. Woolfe, E. Liberty, V. Rokhlin, M. Tygert. A fast randomized
%       algorithm for the approximation of matrices. Appl. Comput. Harmon.
%       Anal. 25, pp. 335-366, 2008.

[m, ~] = size(A);

% Compute subsampled randomized FFT of A
no_sign_switch = binornd(m, .5);
switch_id = randsample(m, no_sign_switch);
A(switch_id, :) = -A(switch_id, :);
FDA = fft(full(A));
Y = FDA(randsample(m, l), :);

% Construct matrix ID of Y
[~, R, e] = qr(Y, 0);
T = R(1:k, 1:k) \ R(1:k, k+1:end);
P = [eye(k) T];
pvec(e) = 1:length(e);
P = P(:, pvec);
J = e(1:k)';

end