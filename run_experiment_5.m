% RUN_EXPERIMENT_5 Run tensor ID experiment for real sparse tensor
%
%   This function requires Tensor Toolbox version 2.6 [Ba15].
%
%   RUN_EXPERIMENT_5 is a script that runs an experiment involving a CP
%   tensor that comes from decomposing a real-world data tensor. The goal
%   of the experiment is to find the maximum magnitude element in the CP
%   tensor using Algorithm 2 in [Re17], with the rank reduction step being
%   done using the following three versions of tensor ID:
%       1.  Tensor ID using Gram matrix [Bi15].
%       2.  Gaussian tensor ID [Bi15].
%       3.  CountSketch tensor ID (proposal).
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
%   [Re17]  M. J. Reynolds, G. Beylkin, A. Doostan. Optimization via
%           separated representations and the canonical tensor
%           decomposition. J. Comput. Phys. 348, pp. 220-230, 2017.

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     May 3, 2019

%% Settings
% data_path: Location of the files containing the decomposition
%   information.
% threshold: Any factor matrix entries of magnitude less than threshold are
%   set to zero
% max_iter: The maximum number of iterations used in the max magnitude
%   finding algorithm.
% no_trial: Number of times the experiment is repeated.
% gaussian_fullrandom: Set to false to avoid generating unnecessary columns
%   of the sketch matrices in the Gaussian tensor ID experiment.

data_path = 'C:\Users\Osman\Desktop\enron_results_10\';
threshold = 1e-6; % We use 1e-6 in our experiment
max_iter = 20;
no_trial = 1;
gaussian_fullrandom = false; % false in our experiment

%% Load the decomposed data
% Note that the "format long e" command is required for Matlab to properly
% read the modeN.mat files, which stores number in this format.

A = cell(4, 1);
format long e
A{1} = readmatrix([data_path 'mode1.mat'], 'filetype', 'text');
A{2} = readmatrix([data_path 'mode2.mat'], 'filetype', 'text');
A{3} = readmatrix([data_path 'mode3.mat'], 'filetype', 'text');
A{4} = readmatrix([data_path 'mode4.mat'], 'filetype', 'text');
for n = 1:4
    A{n} = A{n}(:, 1:end-1);
end
format short
X = ktensor(A);
X = normalize(X);

%% Threshold the data

Y = X;
for k = 1:4
    Y.U{k}(abs(Y.U{k}) < threshold) = 0;
    Y.U{k} = sparse(Y.U{k});
end

rel_error = norm(Y-X)/norm(X);
fprintf('Relative error by thresholding: %.6e\n', full(rel_error));

%% Main loop

addpath(genpath('help_functions'));

STATS = nan(3, 3, no_trial); 
% STATS stores experiment results as follows:
%   - Mode 1 contains 3 statistics: total run time, total sketch time, and
%     total number of iterations.
%   - Mode 2 contains the three methods: Gram, Gaussian and CountSketch.
%   - Mode 3 contains results for each of the no_trials trials.

MAX_ELEMENT = nan(4, 3, no_trial);
% MAX_ELEMENT stores experiment results as follows:
%   - Mode 1 contains each of the four indices of the maximum magnitude
%     element.
%   - Mode 2 contains the three methods: Gram, Gaussian and CountSketch.
%   - Mode 3 contains results for each of the no_trials trials.

for tr = 1:no_trial
    fprintf('Running trial %d of %d... \n', tr, no_trial);
    
    %gram_tic = tic;
    %[Y_gram, sketch_tocs_gram] = ktensor_square_iteration(Y, max_iter, 'gram');
    %for id = 1:4
    %    [~, MAX_ELEMENT(id, 1, tr)] = max(Y_gram.U{id});
    %end
    %gram_toc = toc(gram_tic);
    %STATS(1, 1, tr) = gram_toc;
    %STATS(2, 1, tr) = sum(sketch_tocs_gram(~isnan(sketch_tocs_gram)));
    %STATS(3, 1, tr) = sum(~isnan(sketch_tocs_gram));
    
    gaussian_tic = tic;
    [Y_gaussian, sketch_tocs_gaussian] = ktensor_square_iteration(Y, max_iter, 'gaussian', 'fullrandom', gaussian_fullrandom);
    for id = 1:4
        [~, MAX_ELEMENT(id, 2, tr)] = max(Y_gaussian.U{id});
    end
    gaussian_toc = toc(gaussian_tic);
    STATS(1, 2, tr) = gaussian_toc;
    STATS(2, 2, tr) = sum(sketch_tocs_gaussian(~isnan(sketch_tocs_gaussian)));
    STATS(3, 2, tr) = sum(~isnan(sketch_tocs_gaussian));

    CS_tic = tic;
    [Y_CS, sketch_tocs_CS] = ktensor_square_iteration(Y, max_iter, 'CS');
    for id = 1:4
        [~, MAX_ELEMENT(id, 3, tr)] = max(Y_CS.U{id});
    end
    CS_toc = toc(CS_tic);
    STATS(1, 3, tr) = CS_toc;
    STATS(2, 3, tr) = sum(sketch_tocs_CS(~isnan(sketch_tocs_CS)));
    STATS(3, 3, tr) = sum(~isnan(sketch_tocs_CS)); 

end