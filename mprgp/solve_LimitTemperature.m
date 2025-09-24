[data_i] = load('../linsys_a.dat');
% data_i is [i j value] triplet per line
if size(data_i,2) < 3
	error('Expected 3 columns (i j value) in linsys_a.dat');
end
i = data_i(:,1);
j = data_i(:,2);
v = data_i(:,3);
n = max(max(i), max(j));
A = sparse(i, j, v, n, n);

tb = load('../linsys_b.dat');
if size(tb,2) == 1
	b = tb(:);
else
	b = tb(:,2);
end
if length(b) < n
	b(n,1) = 0;
end

% % Prepare Qpmprgp inputs





% c = -0.5*ones(n,1);
c = load('../lim.dat');


% niq = 1:n;          % mark all indices as unbounded (Qpmprgp sets c(niq) = -Inf)

% % Qpmprgp expects 'b' in workspace for its initial epsr calculation
% b = b0;
% Qpmprgp; x = u;


opts.precond = 'none';    % 'none', 'jacobi', or 'ichol'
opts.epsr = 1e-8;
opts.maxit = 500;
opts.Gamma = 1.0;
[u, info] = mprgp_solver(A, b, c, opts);
x = u;



data = csvread('output.csv', 1, 0);  % skip 1 header row
x_elmer = data(:,1);   % pick the value column
% x_direct = A \ b;

abs_err=norm(x - x_elmer);
rel_err=abs_err / norm(x_elmer);
disp(sprintf('Abs err %.3e, rel err %.3e\\n', abs_err, rel_err));

% abs_err_direct=norm(x - x_direct);
% rel_err_direct=abs_err_direct / norm(x_direct);
% disp(sprintf('Abs err %.3e, rel err %.3e\\n', abs_err_direct, rel_err_direct));


% r = A*x - b; res_rel = norm(r)/norm(b); fprintf('res_rel = %.3e\n', res_rel);
% k = condest(A); fprintf('condest(A) = %.3e\n', k); fprintf('bound for rel error ≲ condest * res_rel = %.3e\n', k * res_rel);


% 1) compare with direct solve
% x_direct = A \ b;
% fprintf('direct residual rel = %.3e, direct-rel-diff = %.3e\n', norm(A*x_direct-b)/norm(b), norm(x_direct-x)/norm(x_direct));
