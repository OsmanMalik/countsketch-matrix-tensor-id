% COMPARISON_ID_VS_CUR Compare conditioning of our proposed ID and CUR. 
%
%   COMPARISON_ID_VS_CUR is a script used to compare the conditioning
%   number of the coefficient matrix produced by our proposed CountSketched
%   matrix ID by the corresponding coefficient matrix produced by the CUR
%   algorithm from [Ma09], specifically the Matlab implementation by
%   [Ch09].
%
% REFERENCES:
%
%   [Ch09]  C. Boutsidis. http://www.boutsidis.org/software.html, accessed
%   25 April 2019.
%
%   [Ma09]  M. W. Mahoney, and P. Drineas. CUR matrix decompositions for
%   improved data analysis. PNAS 160(3), pp. 697–702, 2009.

% Add path to CUR algorithm
addpath('AlgorithmCUR');

% Generate a test matrix
mn = 8;
n = 1000;
k = 100;
A = generate_dense_matrix(n, n, k, mn);

% Compute CS matrix ID
[P, J] = CS_matrix_ID(A, k, k+10, 'qr');

% Compute CUR
[C, U, R] = AlgorithmCUR(A, n, k, n);

% Compute error
err_ID = norm(A - A(:, J)*P);
err_CUR = norm(A - C*U*R);

% Compute condition number of "coefficient matrix"
cond_ID = cond(P);
cond_CUR = cond(U*R);

% Print results
fprintf('Error for ID: %.4e\n', err_ID);
fprintf('Error for CUR: %.4e\n', err_CUR);
fprintf('Condition number for ID coefficient matrix: %.4e\n', cond_ID);
fprintf('Condition number for CUR coefficient matrix: %.4e\n', cond_CUR);