function A = generate_dense_matrix(n, k, mn)
%GENERATE_DENSE_MATRIX Generate dense matrix of certain numerical rank
%  
%   A = GENERATE_DENSE_MATRIX(n, k, mn) returns a randomly generated dense
%   matrix A of size n by n and of rank 2k. The 2-norm of A is 1, and the
%   smallest positive singular value is 10^(-mn). The first k singular 
%   values decay log-linearly. The following k singular values are constant
%   and of magnitude 10^(-mn). This code constructs the output matrix A via
%   A = U*S*V.'

% Construct S, U and V
S = sparse(1:2*k, 1:2*k, 10.^[-((0:k-1)/k)*mn repmat(-mn,1,k)]);
[U, ~] = qr(randn(n, 2*k), 0);
[V, ~] = qr(randn(n, 2*k), 0);

% Compute A
A = U*S*V.';

end

