clear all;
dim = 9;
space = 1/dim;
u = zeros(dim, dim);
f = zeros(dim,dim);
%f(2:end-1, 2:end-1) = 1;
f = smoothRHS(dim);
f_loc = [0; 0; space; space];
u_loc = [0; 0; space; space];

save 'MyInput.mat' '-v4' ;