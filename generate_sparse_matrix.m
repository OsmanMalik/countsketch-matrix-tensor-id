function A = generate_sparse_matrix(n, k, mn, s1, s2)
% GENERATE_SPARSE_MATRIX Generate sparse matrix of certain numerical rank
%
%   A = GENERATE_SPARSE_MATRIX(n, k, mn, s1, s2) returns a randomly
%   generated sparse matrix A of size n by n and of rank 2k. The 2-norm of
%   A is 1, and the smallest positive singular value is 10^(-mn). The
%   variable s1 controls number of nonzero elements in each column of U,
%   and s2 does the same for V. The sparsity of A is nnz(A) = 2*k*s1*s2.
%   This code constructs the matrix A via A = U*S*V.' 

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

if n < 2*k*max(s1, s2)
    error('Must have n >= 2*k*max(s1, s2)')
end

% Create S
S = sparse(1:2*k, 1:2*k, 10.^[-((0:k-1)/k)*mn repmat(-mn,1,k)], 2*k, 2*k);

% Create U
rs = randsample(n, 2*s1*k);
if s1 == 1
    z = randi(2, 1, 2*k)*2-3;
else
    z = randn(s1, 2*k); 
    z = z./sqrt(sum(z.^2));
end
colId = repelem((1:2*k).', s1, 1);
U = sparse(rs, colId, z, n, 2*k);

% Create V
rs = randsample(n, 2*s2*k);
if s2 == 1
    z = randi(2, 1, 2*k)*2-3;
else
    z = randn(s2, 2*k); 
    z = z./sqrt(sum(z.^2));
end
colId = repelem((1:2*k).', s2, 1);
V = sparse(rs, colId, z, n, 2*k);

% Compute A
A = U*S*V.';

end