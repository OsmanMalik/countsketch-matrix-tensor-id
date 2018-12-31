function [P, J] = gram_matrix_ID(A, k)
%GRAM_MATRIX_ID Compute matrix ID via gram matrix
%
%   [P, J] = GRAM_MATRIX_ID(A, k) computes the rank-k matrix ID
%       A \approx A(:, J)*P.
%   This is done by computing a symmetric ID of B = A'*A as explained in
%   Theorem 3.1 of [Bi15] and the related supplementary material. 
%   Throughout the code, all QR factorizations are done using Matlab's QR
%   function.
%
%   Please note that the supplement of [Bi15] contains some minor typos, so
%   the code in this file does not correspond to that paper exactly.
%
% REFERENCES:
%   [Bi15]  D. J. Biagioni, D. Beylkin, G. Beylkin. Randomized 
%           interpolative decomposition of separated representations. J. 
%           Comput. Phys. 281, pp. 116-134, 2015.

% Form gram matrix B and identity I
B = A'*A;
I = eye(size(B));

% Compute relevant portions of matrix ID of B
[~, ~, p] = qr(B, 0);
Pc = I(:, p);
J = p(1:k).';

% Compute unpivoted QR of Bc'*Pc
[~, Rt] = qr(B(:,J)'*Pc, 0);

% Compute S
S = Rt(1:k, 1:k) \ Rt(1:k, k+1:end);

% Compute P
P = [eye(k) S] * Pc';

end