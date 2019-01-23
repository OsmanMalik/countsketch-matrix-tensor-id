function X = generate_tensor(N, I, R, varargin)
%GENERATE_TENSOR Generate dense tensor with known structure
%
%   This function requires Tensor Toolbox [1] version 2.6. 
%
%   X = GENERATE_TENSOR(N, I, R) returns an N-way Tensor Toolbox
%   ktensor, where the n-th factor matrix is created by first creating a
%   matrix of size I(n) by R with iid standard normal entries, and then 
%   normalizing the each column so that each has a 2-norm of 1. The lambda
%   values of X is a vector of length R with r-th entry equal to exp(-r/2).
%   Tensors with these properties are e.g. used in the example of Section
%   5.1 of [2].
%
%   X = GENERATE_TENSOR(___, 'k', k) returns the same tensor as
%   above, but with the r-th lambda entry now equal to exp(-r/k).
%
%   X = GENERATE_TENSOR(___, 'repeat', repeat) returns the same
%   tensor as above, but with the last repeat columns of each factor matrix
%   drawn uniformly at random from the first R-repeat columns. Note that we
%   therefore must have repeat<R. This can be used to create the tensors
%   used in Section 5.1.1 of [2].
%
%   X = GENERATE_TENSOR(___, 'lambda_type', str, 'lambda', lambda) 
%   returns the same tensor as above, but allows us to provide our own
%   vector containing the lambda values. If str is 'exp', then the standard
%   exponentially decaying weights from above are used; if str is 'custom',
%   then the vector lambda is used instead.
%
%   X = GENERATE_TENSOR(___, 'density', d) 
%   returns the same tensor as above, but each random column of each factor
%   matrix is a sparse vector with density d.
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
addParameter(params, 'lambda_type', 'exp')
addParameter(params, 'lambda', ones(R,1));
addParameter(params, 'density', 1, @(x) x > 0 & x <= 1);
parse(params, varargin{:});

k = params.Results.k;
repeat = params.Results.repeat;
lambda_type = params.Results.lambda_type;
lambda_custom = params.Results.lambda;
density = params.Results.density;

%% Create tensor

A = cell(N,1);
idx = randi(R-repeat, repeat, 1);
for n = 1:N
    if density == 1
        A{n} = randn(I(n), R-repeat);
    else
        t = ceil(I(n)*density);
        rid = randi(I(n), t*(R-repeat), 1);
        cid = repelem((1:R-repeat).', t, 1);
        val = randn(t*(R-repeat), 1); 
        A{n} = sparse(rid, cid, val, I(n), R-repeat);
        %for r = 1:R-repeat
        %    A{n}(:, r) = sprand(I(n), 1, density);
        %end
    end
    A{n} = A{n}./sqrt(sum(A{n}.^2, 1));
    A{n} = [A{n} A{n}(:, idx)];
end

if strcmp(lambda_type, 'exp')
    lambda = exp(-(1:R)/k).';
elseif strcmp(lambda_type, 'custom')
    if iscolumn(lambda_custom)
        lambda = lambda_custom;
    else
        lambda = lambda_custom.';
    end
end
X = ktensor(lambda, A);

end