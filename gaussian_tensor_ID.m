function Xk = gaussian_tensor_ID(X, k, l, QR_type)
% GAUSSIAN_TENSOR_ID Computes Gaussian tensor ID
%
%   This function requires Tensor Toolbox version 2.6 [Ba15].
%   
%   Xk = GAUSSIAN_TENSOR_ID(X, k, l, QR_type) returns a rank-k tensor ID of
%   the input tensor X computed using an oversampling parameter l-k. Note
%   that we therefore require l >= k. The computation is done according to
%   Algorithm 3 of [Bi15]. QR_type controls which type of QR factorization
%   is used: Set it to 'srrqr' to use the strong rank-revealing QR
%   factorization of [Gu96] (uses the implementation [Xi18]), or set it to
%   'qr' to use Matlab's built-in QR function.
%
% REFERENCES:
%
%   [Ba15]  B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%           Version 2.6, Available online, February 2015. 
%           URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.
%   
%   [Bi15]  D. J. Biagioni, D. Beylkin, G. Beylkin. Randomized
%           interpolative decomposition of separated representations. J.
%           Comput. Phys. 281, pp. 116-134, 2015.
%
%   [Gu96]  M. Gu, and S. C. Eisenstat. Efficient algorithms for computing
%           a strong rank-revealing QR factorization. SIAM J. Sci. Comput.
%           17(1), pp. 848-869, 1996.
%
%   [Xi18]  X. Xing. Interpolative Decomposition based on Strong RRQR.
%           MATLAB Central File Exchange. Retrieved November 23, 2018.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

N = ndims(X);
I = size(X);
R = ncomponents(X);

% Compute projection Y (Steps 1 and 2 in Alg. 3 of [Bi15])
Y = ones(l, R);
for n = 1:N
    
    %fprintf('MODIFIED GAUSSIAN\n')
    %G = zeros(l, I(n));
    %nnzidx = sum(abs(X.U{n}), 2) ~= 0;
    %G(:, nnzidx) = randn(l, sum(nnzidx));
    
    %fprintf('OLD GAUSSIAN\n')
    G = randn(l, I(n));
    
    Y = Y .* (G*X.U{n});
end
Y = repmat(X.lambda.', l, 1) .* Y;

% Compute rank-k matrix ID of Y (Step 3 in Alg. 3 of [Bi15])
if strcmp(QR_type, 'srrqr')
    f = 2;
    [P, J] = ID(Y.', 'rank', k, f);
    P = P.';
elseif strcmp(QR_type, 'qr')
    [~, R, e] = qr(Y, 0);
    k = min(k, rank(R));
    T = R(1:k, 1:k) \ R(1:k, k+1:end);
    P = [eye(k) T];
    pvec(e) = 1:length(e);
    P = P(:, pvec);
    J = e(1:k)';
end

% Form the rank-k tensor ID Xk (Step 4 in Alg. 3 of [Bi15])
A = cell(N, 1);
for n = 1:N
    A{n} = X.U{n}(:, J);
end
%alpha = P*X.lambda;
alpha = X.lambda(J).*sum(P, 2);
Xk = ktensor(alpha, A);

end