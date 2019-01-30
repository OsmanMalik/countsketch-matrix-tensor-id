function [P, J] = SRFT_matrix_ID(A, k, l, varargin)
% SRFT_MATRIX_ID Computes subsampled randomized FFT matrix ID
%
%   [P, J] = SRFT_MATRIX_ID(A, k, l) computes the rank-k matrix ID
%       A \approx A(:, J)*P,
%   where A(:, J) contains k columns of A. The argument l > k is k plus the
%   oversampling parameter. Before applying the matrix ID, we sketch the 
%   matrix using a subsampled randomized FFT as proposed in [Wo08].
%
%   [P, J] = SRFT_MATRIX_ID(___, 'splits', no_splits) splits up the
%   computation by only computing the SRFT sketch for a subset of 
%   size(A, 2)/no_splits columns of A at a time, thus avoiding make the
%   entire matrix A into a dense matrix in one go, in the case that A is
%   sparse. Note that no_splits must be a positive integer which divides
%   the number of columns of A.
%
% REFERENCES:
%
%   [Wo08]  F. Woolfe, E. Liberty, V. Rokhlin, M. Tygert. A fast randomized
%           algorithm for the approximation of matrices. Appl. Comput.
%           Harmon. Anal. 25, pp. 335-366, 2008.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

[m, n] = size(A);

% Handle optional inputs
params = inputParser;
addParameter(params, 'splits', 1, @(x) isscalar(x) & x > 0 & mod(n, x) == 0);
parse(params, varargin{:});

no_splits = params.Results.splits;

% Compute subsampled randomized FFT of A
no_sign_switch = binornd(m, .5);
switch_id = randsample(m, no_sign_switch);
A(switch_id, :) = -A(switch_id, :);
Y = zeros(l, n);
split_sz = n/no_splits;
sub_sample = randsample(m, l);
for split = 1:no_splits
    FDA = fft(full(A(:, 1+(split-1)*split_sz:split*split_sz)));
    Y(:, 1+(split-1)*split_sz:split*split_sz) = FDA(sub_sample, :);
end

% Construct matrix ID of Y
[~, R, e] = qr(Y, 0);
T = R(1:k, 1:k) \ R(1:k, k+1:end);
P = [eye(k) T];
pvec(e) = 1:length(e);
P = P(:, pvec);
J = e(1:k)';

end