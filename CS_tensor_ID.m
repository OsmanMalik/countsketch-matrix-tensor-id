function Xk = CS_tensor_ID(X, k, l, QR_type)
% CS_TENSOR_ID Computes CountSketched tensor ID
%
%   This function requires Tensor Toolbox version 2.6 [Ba15].
%
%   Xk = CS_TENSOR_ID(X, k, l, QR_type) returns a rank-k tensor ID of the
%   input tensor X computed using an oversampling parameter l-k. Note that
%   we therefore require l >= k. The computation is done by applying
%   TensorSketch. Tensor ID was first proposed in [Bi15].
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

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 29, 2019

% Get dimensions of X
N = ndims(X);
R = ncomponents(X);

% Compute TensorSketch of matrix corresponding to X
A = X.U;
A{1} = A{1}*sparse(1:R, 1:R, X.lambda, R, R);
Z = TensorSketch(A, l);

% Proceed by using appropriate QR factorization
if strcmp(QR_type, 'srrqr')
    f = 2;
    [P, J] = ID(Z.', 'rank', k, f);
    P = P';
elseif strcmp(QR_type, 'qr')
    [~, R, e] = qr(Z, 0);
    k = min(k, rank(R)); % Added this line to avoid issues when computing T when rank(Z) = rank(R) < k 
    T = R(1:k, 1:k) \ R(1:k, k+1:end);
    P = [eye(k) T];
    pvec(e) = 1:length(e);
    P = P(:, pvec);
    J = e(1:k)';
end

% Form the rank-k tensor ID Xk
B = cell(N, 1);
for n = 1:N
    B{n} = X.U{n}(:, J);
end
%alpha = P*X.lambda;
alpha = X.lambda(J).*sum(P, 2);
Xk = ktensor(alpha, B);

end