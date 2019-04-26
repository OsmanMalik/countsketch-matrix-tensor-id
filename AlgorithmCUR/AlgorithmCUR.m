function [C, U, R] = AlgorithmCUR(A, k, c, r)
%
% [C, U, R] = AlgorithmCUR(A, k, c, r)
%
% This function implements a practical version of the CUR algorithm descri-
% bed in Ref. [1]. The main difference between the description of the algo-
% rithm of [1] and the implementation described here is the fact that c and 
% r are given as inputs from the user while e is considered as a constant, 
% thus ignored from the input of the algorithm. The rest of the implementa-
% tion follows the description given in [1].
% -------------------------------------------------------------------------
%
% Input: 
%    - A: m x n matrix.
%    - k: rank parameter, with k << min(m,n). 
%    - c: number of columns that we want to select from A.  
%    - r: number of rows that we want to select from A. 
%      
% Ouput:
%    - C: m x c' matrix with c' columns from A. E(c') <= c. 
%    - R: r' x n matrix with r' rows form A.    E(r') <= r.
%    - U: c' x r' matrix.
%
% Outcome: 
%    The rank-t matrix T = CUR, with t = min(c',r'), is a ''good''approxima-
%    tion of the input matrix A, i.e. ||A - T||_F is ''close'' to ||A - Ak||_F,
%    where Ak is the best rank-k approximation of A computed with the Singular 
%    Value Decomposition. Since e is not part of the input of the algorithm 
%    and c and r are arbitrary chosen numbers, we do not expect that eqn. 5 
%    of [1] holds for any matrix A. Instead, in the implementation described 
%    here the user specifies the number of columns and rows to be selected 
%    and the algorithm selects the ''best'' such columns and rows in order 
%    to make the error ||A - T||_F  as small as possible.
%
% References:
% [1] Michael W. Mahoney and Petros Drineas
%     CUR matrix decompositions for improved data analysis
%     PNAS 2009 106: 697-702.
%
% -------------------------------------------------------------------------
% Author: Christos Boutsidis. March 2009. 
% E-mail: christos.boutsidis@gmail.com
% -------------------------------------------------------------------------

% Compute the right singular vectors of A and A'
% - v: n x k matrix with the top-k right singular vectors of A.
% - u: m x k matrix with the top-k right singular vectors of A'.
[u s v]         = svds(A, k) ;  

% The CUR Algorithm
C = ColumnSelect(A, k, c, v)        ;    % Choose c' columns from A. 
R = ColumnSelect(A', k, r, u)'      ;    % Choose r' rows from A.
U = pinv(C, .05) * A * pinv(R, .05) ;    % Compute U.         
return
%--------------------------------------------------------------------------

% Here is an implementation of the Algorithm ColumnSelect described in [1].
%--------------------------------------------------------------------------
function  C = ColumnSelect(A, k, c, v) 
%
% Input
%   - A: m x n matrix.
%   - k: rank parameter k.
%   - c: number of columns that we want to select from A.
%   - v: n x k matrix of the top-k right singular vectors of A.

% Output
%   - C: m x c' matrix with c' columns from A, E(c') <= c.

[m n] = size(A)  ; % the size of the input matrix A.

%------- Compute the normalized leverage scores of eqn. 3 of [1]. ---------
pi = zeros(1, n) ; 
for j=1:n
    pi(j) =  (norm(v(j,:))^2) / k  ;
end
%--------------------------------------------------------------------------

%---------------- randomized column selection -----------------------------

indexA = []; % indexA is initially empty. 

for j=1:n    % for every column of A
    
    % the j-th column of A is selected with probability prob_j.
    prob_j = min([1 c*pi(j)]);  % find the minimum of 1 and  c*pi(j)
    prob_j = prob_j(1);         % resolve the case where 1 = c*pi(j)
        
    if prob_j==1             % if prob_j=1 select the j-th column of A
        indexA = [indexA j];
    elseif  prob_j > rand    % if prob_j<1, generate a random number rand in [0,1] and then 
        indexA = [indexA j]; % if prob_j > rand, select the j-th column of A 
    end
    
end

% At the end of this process indexA contains the indices of the selected 
% columns of A, i.e. C = A(:, indexA);

C = A(:, indexA);
return
%--------------------------------------------------------------------------