function Xk = gram_tensor_ID(X, k)
%GRAM_TENSOR_ID Compute tensor ID via gram matrix
%
%   This function requires Tensor Toolbox [Ba15] version 2.6.
%
%   Xk = GRAM_TENSOR_ID(X, k) returns a rank-k tensor ID of the input
%   tensor X computed using the gram matrix described in equation (1.12) in
%   [Bi15]. Throughout the code, all QR factorizations are done using
%   Matlab's QR function.
%
%   Please note that the supplement of [Bi15] contains some minor typos, so
%   the code in this file does not correspond to that paper exactly.
%
% REFERENCES:
%   [Ba15]  B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%           Version 2.6, Available online, February 2015. 
%           URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.
%   
%   [Bi15]  D. J. Biagioni, D. Beylkin, G. Beylkin. Randomized 
%           interpolative decomposition of separated representations. J. 
%           Comput. Phys. 281, pp. 116-134, 2015.

N = ndims(X);
R = ncomponents(X);
Id = eye(R);

% Compute gram matrix G (Eq (1.12) in [Bi15])
G = ones(R, R);
for n = 1:N
    G = G .* (X.U{n}.'*X.U{n});
end
G = G .* (X.lambda*X.lambda.');

% Compute relevant portions of matrix ID of G
[~, ~, p] = qr(G, 0);
Pc = Id(:, p);
J = p(1:k).';

% Compute unpivoted QR of Gc'*Pc
[~, Rt] = qr(G(:,J)'*Pc, 0);

% Compute S
S = Rt(1:k, 1:k) \ Rt(1:k, k+1:end);

% Compute P
P = [eye(k) S] * Pc';

% Form the rank-k tensor ID Xk (see e.g. Step 4 in Alg. 3 of [Bi15])
A = cell(N, 1);
for n = 1:N
    A{n} = X.U{n}(:, J);
end
alpha = P*X.lambda;
Xk = ktensor(alpha, A);

end