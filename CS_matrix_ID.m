function [P, J] = CS_matrix_ID(A, k, l, QR_type)
%CS_MATRIX_ID Computes CountSketched matrix ID
%
%   [P, J] = CS_MATRIX_ID(A, k, l, QR_type) computes the rank-k matrix ID
%       A \approx A(:, J)*P,
%   where A(:, J) contains k columns of A. The argument l > k is k plus the
%   oversampling parameter. QR_type controls the type of QR factorization
%   used: Set it to 'srrqr' to use the strong rank-revealing QR
%   factorization of [1] (uses the implementation [2]), or set it to 'qr'
%   to use Matlab's built-in QR function. Before applying the matrix ID, we
%   sketch the matrix using a CountSketch matrix.
%
%   REFERENCES:
%   [1] M. Gu, and S. C. Eisenstat. Efficient algorithms for computing a
%       strong rank-revealing QR factorization. SIAM J. Sci. Comput. 17(1),
%       pp. 848-869, 1996.
%
%   [2] X. Xing. Interpolative Decomposition based on Strong RRQR. MATLAB
%       Central File Exchange. Retrieved November 23, 2018.

[m, ~] = size(A);

% Define hash function
h = [1:l randi(l, 1, m-l)].';
h = h(randperm(m));
s = randi(2, m, 1)*2-3;

% Compute CountSketch of A
if issparse(A)
    Z = countSketch_sparse(A.', int64(h), l, s).';
else
    Z = countSketch(A.', int64(h), l, s, 1).';
end

% Old code
%{
% Switch signs of some rows of A
no_sign_switch = binornd(m, .5);
switch_id = randsample(m, no_sign_switch);
A(switch_id, :) = -A(switch_id, :);

% Make A smaller by adding rows
randvec = [randsample(l,l); randi(l, m-l, 1)];
randvec = randvec(randperm(m));
CS = sparse(randvec, 1:m, ones(m,1), l, m);
Z = CS*A;
%}

% Proceed by using appropriate QR factorization
if strcmp(QR_type, 'srrqr')
    f = 2;
    [P, J] = ID(Z', 'rank', k, f);
    P = P';
elseif strcmp(QR_type, 'qr')
    [~, R, e] = qr(Z, 0);
    T = R(1:k, 1:k) \ R(1:k, k+1:end);
    P = [eye(k) T];
    pvec(e) = 1:length(e);
    P = P(:, pvec);
    J = e(1:k)';
end

end