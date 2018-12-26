function Xk = gaussian_tensor_ID(X, k, l)
%GAUSSIAN_TENSOR_ID Computes Gaussian tensor ID
%
%   Xk = GAUSSIAN_TENSOR_ID(X, k, l)
%
% REFERENCES:
%   [1] B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%       Version 2.6, Available online, February 2015. 
%       URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.
%   
%   [2] D. J. Biagioni, D. Beylkin, G. Beylkin. Randomized interpolative
%       decomposition of separated representations. J. Comput. Phys. 281, 
%       pp. 116-134, 2015.

N = ndims(X);
I = size(X);
R = ncomponents(X);
f = 2;

% Compute projection Y (Steps 1 and 2 in Alg. 3 of [2])
Y = ones(l, R);
for n = 1:N
    G = randn(l, I(n));
    Y = Y .* (G*X.U{n});
end
Y = repmat(X.lambda.', l, 1) .* Y;

% Compute rank-k matrix ID of Y (Step 3 in Alg. 3 of [2])
[P, J] = ID(Y.', 'rank', k, f);
P = P.';

% Form the rank-k tensor ID Xk (Step 4 in Alg. 3 of [2])
A = cell(N, 1);
for n = 1:N
    A{n} = X.U{n}(:, J);
end
alpha = P*X.lambda;
Xk = ktensor(alpha, A);

end