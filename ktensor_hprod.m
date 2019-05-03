function C = ktensor_hprod(A, B)
% KTENSOR_HPROD Compute Hadamard product of two ktensors with sparse factor
% matrices.
%
%   This function requires Tensor Toolbox version 2.6 [Ba15].
%
%   C = KTENSOR_HPROD(A, B) returns a ktensor C such that C is the Hadamard
%   product of A and B, with 
%   ncomponents(C) = ncomponents(A)*ncomponents(B).
%
% REFERENCES:
%
%   [Ba15]  B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%           Version 2.6, Available online, February 2015. 
%           URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     April 27, 2019

N1 = ndims(A);
N2 = ndims(B);
if N1 ~= N2
    error('Input tensors must have same number of dimensions');
end

lambda = kron(A.lambda, B.lambda);
U = cell(N1, 1);
for n = 1:N1
    if issparse(A.U{n}) && issparse(B.U{n})
        U{n} = sparse_khatrirao(A.U{n}.', B.U{n}.').';
    else
        U{n} = khatrirao(A.U{n}.', B.U{n}.').';
    end
end

C = ktensor(lambda, U);

end