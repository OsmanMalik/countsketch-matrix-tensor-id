function save_matrix_to_file(A, bin_file, verbosity)
%SAVE_MATRIX_TO_FILE Writes matrix to file
% 
%   SAVE_MATRIX_TO_FILE(A, bin_file) writes the matrix A to a file
%   specified in bin_file. It writes it in a format which is compatible
%   with RSVDPACK [1]. Note that the matrix is stored in ROW MAJOR order,
%   which is how the generated matrices in RSVDPACK are stored.
%
% REFERENCES:
%   [1] S. Voronin, and P. G. Martinsson. RSVDPACK: An implementation of 
%       randomized algorithms for computing the singular value, 
%       interpolative, and CUR decompositions of matrices on multi-core and
%       GPU architectures. arXiv:1502.05366v3 [math.NA], 2016.

[m, n] = size(A);

fp = fopen(bin_file, 'w');
fwrite(fp, m, 'int32'); % Write no. rows of A
fwrite(fp, n, 'int32'); % Write no. columns of A
for id1 = 1:m
    if verbosity >= 2
        fprintf('Writing row %d of %d... ', id1, n);
    end
    for id2 = 1:n
        fwrite(fp, A(id1, id2), 'double');
    end
    if verbosity >= 2
        fprintf('Done!\n');
    end
end
fclose(fp);

end

