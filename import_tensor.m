% IMPORT_TENSOR Import tns data and save a smaller version.
%
%   Before running this script, download one of the tns files from
%   http://frostt.io/, e.g. enron.tns or delicious4d.tns.
%
%   IMPORT_TENSOR is a script that loads the data in a tns file and then
%   creates and saves a smaller version of this tensor. The reduced size
%   tensor is saved in the format compatible with SPLATT [Sh19].
%
%   Note that the writematrix function was introduced in Matlab R2019a.
%
% REFERENCES:
%
%   [Sh19]  Shaden Smith. https://github.com/ShadenSmith/splatt, accessed
%           April 28, 2019.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     April 28, 2019

%tensor_size = [6066 5699 244268 1176];
tensor_size = [532924 17262471 2480308 1443];

fprintf('Loading data... ');

%A = importdata('enron.tns', ' ');
A = importdata('delicious4d.tns', ' ');

fprintf('Done!\n');
fprintf('Creating sparse tensor... ');

X = sptensor(A(:, 1:4), A(:, 5), tensor_size);

fprintf('Done!\n')
fprintf('Sparsity of full tensor is %.4e\n', nnz(X)/prod(tensor_size));
fprintf('Reducing size of tensor... ');

X = X(:, 1:1e+6, 1:1e+6, :);

fprintf('Done!\n')
fprintf('Sparsity of reduced tensor is %.4e\n', nnz(X)/prod(tensor_size));
fprintf('Writing reduced tensor to file... ');

writematrix([X.subs X.vals], 'delicious4d_small.tns', 'delimiter', ' ', 'filetype', 'text');

fprintf('Done!\n');