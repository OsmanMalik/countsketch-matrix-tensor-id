function Xk = CS_tensor_ID(X, k, l, QR_type)
%CS_TENSOR_ID Computes CountSketched tensor ID
%
%   This function requires Tensor Toolbox [1] version 2.6.
%
%   Xk = CS_TENSOR_ID(X, k, l, QR_type) returns a rank-k tensor ID of the
%   input tensor X computed using an oversampling parameter l-k. Note that
%   we therefore require l >= k. The computation is done by applying
%   TensorSketch. 
%
% REFERENCES:
%   [1] B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%       Version 2.6, Available online, February 2015. 
%       URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.

% Get dimensions of X
N = ndims(X);
R = ncomponents(X);

% Define hash functions
h = cell(N, 1);
s = cell(N, 1);
for n = 1:N
    h{n} = randi(l, size(X, n), 1);
    s{n} = randi(2, size(X, n), 1)*2-3;
end

% Compute TensorSketch of matrix corresponding to X
A = X.U;

A{1} = A{1}*sparse(1:R, 1:R, X.lambda, R, R);
%A{1} = repmat(X.lambda.', size(A{1}, 1), 1) .* A{1};
Z = TensorSketch(A, h, s, l);

% Proceed by using appropriate QR factorization
if strcmp(QR_type, 'srrqr')
    f = 2;
    [P, J] = ID(Z.', 'rank', k, f);
    P = P';
elseif strcmp(QR_type, 'qr')
    [~, R, e] = qr(Z, 0);
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
alpha = P*X.lambda;
Xk = ktensor(alpha, B);

end