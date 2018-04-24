clear all;

dim = 20;
u = zeros(dim, dim);
f = zeros(dim, dim);
f(2:end-1, 2:end-1) = 1;
f_loc = [0; 0; 0.025; 0.025];
u_loc = [0; 0; 0.025; 0.025];

save 'MyInput.mat' '-v4' ;