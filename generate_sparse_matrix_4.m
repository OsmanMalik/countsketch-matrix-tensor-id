function A = generate_sparse_matrix_4(m, n, k, mn, rho)
%GENERATE_SPARSE_MATRIX_4 Generates a sparse matrix
%
%   A = GENERATE_SPARSE_MATRIX_4(m, n, k, mn, rho) returns a randomly
%   generated sparse matrix A of size m by n. A is constructed as A =
%   U*S*V.', where U and V are random sparse matrices, with each column a
%   random sparse vector with density rho normalized to have unit 2-norm. S
%   is chose so that the first k elements along the diagonal of S decay
%   exponentially to size 1e-mn, and the following k elements constant of
%   size 1e-mn.

U = sparse(m, 2*k);
V = sparse(n, 2*k);
for idx = 1:2*k
    U(:,idx) = sprand(m, 1, rho);
    V(:,idx) = sprand(n, 1, rho);
end
U(U~=0) = U(U~=0)*2-1;
V(V~=0) = V(V~=0)*2-1;
U = U./sqrt(sum(U.^2,1));
V = V./sqrt(sum(V.^2,1));
S = sparse(1:2*k, 1:2*k, 10.^[-((0:k-1)/k)*mn repmat(-mn,1,k)], 2*k, 2*k);

A = U*S*V.';

end