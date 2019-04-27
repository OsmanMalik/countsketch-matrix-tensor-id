function Y = ktensor_square_iteration(X, max_iter)
% KTENSOR_SQUARE_ITERATION Repeatedly square and normalize the elements of
% the ktensor X.
%
%   This function requires Tensor Toolbox version 2.6 [Ba15].
%
%   Y = KTENSOR_SQUARE_ITERATION(X, max_iter) repeatedly squares the
%   elements of X and normalizes them, and then returns the result in Y
%   when the squared tensor has become rank 1. The value max_iter sets the
%   maximum number of iterations. This is Algorithm 2 in [Re17], but with a
%   slightly different normalization.
%
% REFERENCES:
%
%   [Ba15]  B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%           Version 2.6, Available online, February 2015. 
%           URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.
%
%   [Re17]  M. J. Reynolds, G. Beylkin, A. Doostan. Optimization via
%           separated representations and the canonical tensor
%           decomposition. J. Comput. Phys. 348, pp. 220-230, 2017.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     April 27, 2019

Y = X;
Y = normalize(Y);
Y = (1/max(Y.lambda))*Y;
K = ncomponents(Y);
L = K;

for it = 1:max_iter
    Q = ktensor_hprod(Y, Y);
    Y = CS_tensor_ID(Q, K, L, 'qr');
    Y = normalize(Y);
    Y = (1/max(Y.lambda))*Y;
    fprintf('Finished iteration %d. Maximum lambda is %.3e.\n', it, max(Y.lambda));
    if ncomponents(Y) == 1
        fprintf('Y is rank 1. Breaking...\n\n');
        break
    end
end

fprintf('Max iterations reached. Breaking...\n\n');

end