function SFA = ssfft(A, l, S)
%SSFFT Computes subsampled FFT.

[no_row, no_col] = size(A);
n = no_col;

% Round l to nearest power of 2
lp = 2^nextpow2(l);

% Add rows if necessary to make number of rows a power of 2
no_row_p = 2^nextpow2(no_row);
if no_row_p ~= no_row
    A = [A; zeros(no_row_p-no_row, n)];
end

% Compute variable mp such that number of rows in padded matrix is mp*lp
mp = no_row_p/lp;

% Store each column of padded A so that A = [a1, a2, ..., an], where each
% ai is a matrix of size (lp x mp) storing column i of A in row major order
VT = reshape(A, mp, lp*n);
V = zeros(lp, mp*n);
for id = 1:n
    V(:, 1+mp*(id-1):mp*id) = VT(:, 1+lp*(id-1):lp*id).';
end

% Step 1 of algorithm
W = fft(V);

% Step 2 of algorithm
J = repmat((1:lp).', 1, mp);
K = repmat(1:mp, lp, 1);
X = repmat(exp(-2*pi*1i*(J-1).*(K-1)/no_row), 1, n) .* W;

% Step 3 of algorithm
Y = zeros(mp, lp*n);
for id = 1:n
    Y(:, 1+lp*(id-1):lp*id) = X(:, 1+mp*(id-1):mp*id).';
end

if isempty(S)
    % Step 4 of algorithm
    Z = fft(Y);
    ZT = zeros(lp, mp*n);
    for id = 1:n
        ZT(:, 1+mp*(id-1):mp*id) = Z(:, 1+lp*(id-1):lp*id).';
    end

    % Reshape output and subsample if necessary
    SFA = reshape(ZT, no_row_p, no_col);
    if no_row_p ~= no_row
        SFA = SFA(1:no_row, :);
    end
else
    [Z_row, Z_col] = ind2sub([lp mp], S);
    SFA = zeros(length(S), no_col);
    dft_mat = dftmtx(mp);
    for z = 1:length(Z_row)
        SFA(z, :) = dft_mat(Z_col(z), :) * Y(:, (0:n-1)*lp + Z_row(z));
    end
end

end