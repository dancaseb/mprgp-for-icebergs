output_my_solver = csvread('output_my_solver.csv', 1, 0);  % skip 1 header row
output_default_solver = csvread('output_default_solver.csv', 1, 0);  % skip 1 header row

x_elmer = output_default_solver(:,1); 
x_my_solver = output_my_solver(:,1);

abs_err=norm(x_my_solver - x_elmer);
rel_err=abs_err / norm(x_elmer);
disp(sprintf('Abs err %.3e, rel err %.3e\\n', abs_err, rel_err));