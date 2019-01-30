function nrm = rand_norm(A, n)
% RAND_NORM Computes the randomized 2-norm of a matrix A
%
%   nrm = RAND_NORM(A, n) computes the randomized 2-norm of the matrix A,
%   using n random vectors. For details, including theoretical guarantees,
%   please see Section 3.4 of [Wo08].
%
% REFERENCES:
%
%   [Wo08]  F. Woolfe, E. Liberty, V. Rokhlin, M. Tygert. A fast randomized
%           algorithm for the approximation of matrices. Appl. Comput.
%           Harmon. Anal. 25, pp. 335-366, 2008.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

X = randn(size(A, 2), n);
X = X./sqrt(sum(X.^2));
Y = A*X;
nrm = max(sqrt(sum(Y.^2)));

end