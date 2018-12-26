function X = generate_dense_tensor(N, I, R, varargin)
%GENERATE_DENSE_TENSOR Generate dense tensor with known structure
%
%   This function requires Tensor Toolbox [1] version 2.6. 
%
%   X = GENERATE_DENSE_TENSOR(N, I, R) returns an N-way Tensor Toolbox
%   ktensor, where the n-th factor matrix is created by first creating a
%   matrix of size I{n} by R with iid standard normal entries, and then 
%   normalizing the each column so that each has a 2-norm of 1. The lambda
%   values of X is a vector of length R with r-th entry equal to exp(-r/2).
%   Tensors with these properties are e.g. used in the example of Section
%   5.1 of [2].
%
%   X = GENERATE_DENSE_TENSOR(___, 'k', k) returns the same tensor as
%   above, but with the r-th lambda entry now equal to exp(-r/k).
%
%   X = GENERATE_DENSE_TENSOR(___, 'repeat', repeat) returns the same
%   tensor as above, but with the last repeat columns of each factor matrix
%   drawn uniformly at random from the first R-repeat columns. Note that we
%   therefore must have repeat<R. This can be used to create the tensors
%   used in Section 5.1.1 of [2].
%   
% REFERENCES:
%   [1] B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%       Version 2.6, Available online, February 2015. 
%       URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.
%   
%   [2] D. J. Biagioni, D. Beylkin, G. Beylkin. Randomized interpolative
%       decomposition of separated representations. J. Comput. Phys. 281, 
%       pp. 116-134, 2015.

%% Handle optional inputs

params = inputParser;
addParameter(params, 'k', 2);
addParameter(params, 'repeat', 0, @(x) isscalar(x) & x < R);
parse(params, varargin{:});

k = params.Results.k;
repeat = params.Results.repeat;

%% Create tensor

A = cell(N,1);
idx = randi(R-repeat, repeat, 1);
for n = 1:N
    A{n} = randn(I{n}, R-repeat);
    A{n} = A{n}./sqrt(sum(A{n}.^2, 1));
    A{n} = [A{n} A{n}(:, idx)];
end

lambda = exp(-(1:R)/k).';
X = ktensor(lambda, A);

end