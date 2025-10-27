% check_spd_strict.m
% Read "linsys_a.dat" (i j value triples), check strict symmetry (within tol),
% and only if symmetric run positive-definite test (Cholesky).
%
% This script WILL NOT symmetrize the matrix. It only checks and reports.

filename = 'linsys_a.dat';
tol = 1e-12;   % tolerance for symmetry check (set slightly above machine eps if desired)

% ---------- read file ----------
fid = fopen(filename, 'r');
if fid < 0
    error('Cannot open file %s', filename);
end
C = textscan(fid, '%d %d %f', 'CommentStyle', '#');
fclose(fid);

if isempty(C) || isempty(C{1})
    error('No data read from %s', filename);
end

i = C{1}(:);
j = C{2}(:);
v = C{3}(:);

if ~(numel(i) == numel(j) && numel(j) == numel(v))
    error('Input file parsing failed: column lengths differ (%d, %d, %d).', numel(i), numel(j), numel(v));
end

n = max(max(i), max(j));
A = sparse(i, j, v, n, n);

% ---------- quick basic checks ----------
fprintf('Matrix size: %d x %d, nnz = %d\n', n, n, nnz(A));

% If the matrix has tiny imaginary parts (shouldn't for data file), warn.
if ~isreal(A)
    max_im = max(abs(imag(A(:))));
    fprintf('Warning: matrix has imaginary parts (max imag = %.3e). This check requires a real matrix.\n', max_im);
    if max_im > (10 * tol)
        error('Imaginary parts too large; aborting.');
    else
        fprintf('Imag parts are tiny; continuing using real(A).\n');
        A = real(A);
    end
end

% ---------- symmetry check (strict, no symmetrization) ----------
% compute elementwise difference A - A'
D = A - A';
max_abs_diff = full(max(abs(D(:))));
norm_inf = norm(D, inf);

fprintf('Symmetry diagnostics:\n');
fprintf('  max |A - A''| = %.5e\n', max_abs_diff);
fprintf('  ||A - A''||_inf = %.5e\n', norm_inf);

if max_abs_diff > tol
    % Find location(s) of largest mismatch to help debugging
    [ri, rj, rv] = find(abs(D) == max_abs_diff);
    % If multiple, show up to first 10
    m = min(length(ri), 10);
    fprintf('Matrix is NOT symmetric within tol = %.3e. Showing up to %d largest mismatch locations (row col diff):\n', tol, m);
    for k = 1:m
        r = ri(k); c = rj(k);
        fprintf('  (%d, %d): A(%d,%d)=%.12g  A(%d,%d)=%.12g  diff=%.5e\n', ...
                r, c, r, c, full(A(r,c)), c, r, full(A(c,r)), full(A(r,c) - A(c,r)));
    end
    fprintf('\nRESULT: Original matrix is NOT symmetric -> cannot be positive definite in the strict sense.\n');
    return;   % stop here (no symmetrization, no PD test)
else
    fprintf('Matrix is symmetric within tol = %.3e. Proceeding to PD test.\n', tol);
end

% ---------- positive-definite test ----------
% Use Cholesky. In Octave: [R,p] = chol(A) returns p==0 on success.
try
    [R, p] = chol(A);
catch
    % chol might error for non-symmetric matrices; treat as failure
    p = 1;
end

if exist('p','var') && p == 0
    fprintf('Cholesky succeeded -> Matrix IS symmetric positive definite.\n');
else
    fprintf('Cholesky failed -> Matrix is NOT positive definite (p = %d).\n', p);
    % For diagnostics: try to return a few smallest eigenvalues if feasible
    try
        k = min(6, n);
        vals = eigs(A, k, 'SM');   % may fail for ill-conditioned or non-symmetric but we are symmetric here
        vals = sort(real(vals));
        fprintf('Approx smallest eigenvalues (via eigs):\n');
        disp(vals);
    catch
        if n <= 2000
            ev = eig(full(A));
            ev = sort(real(ev));
            m = min(6, numel(ev));
            fprintf('Exact smallest eigenvalues:\n');
            disp(ev(1:m));
        else
            fprintf('Could not compute eigenvalue diagnostics (eigs failed and matrix too large for full eig).\n');
        end
    end
end

% export nothing else (no symmetrized matrix). If you want A in workspace:
assignin('base', 'A_strict', A);
