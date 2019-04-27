function X = ktensor_to_double(X)
% KTENSOR_TO_DOUBLE Transform a ktensor to a double array.
%
%   This function requires Tensor Toolbox version 2.6 [Ba15].
%   
%   X = KTENSOR_TO_DOUBLE(X) transforms the input ktensor X into a double
%   array.
%
% REFERENCES:
%
%   [Ba15]  B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox 
%           Version 2.6, Available online, February 2015. 
%           URL: http://www.sandia.gov/~tgkolda/TensorToolbox/.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     April 27, 2019

for n = 1:ndims(X)
    X.U{n} = full(X.U{n});
end

X = double(X);

end