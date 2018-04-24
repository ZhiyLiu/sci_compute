clear all;

dimList = [129];
NList = length(dimList);
na = 3; nb = 3;
omega = 0.8;


for i = 1:NList
    rand('seed', 23423);
    
    dim = dimList(i);
    u = zeros(dim, dim);
    f = zeros(dim, dim);
    f(2:end-1, 2:end-1) = 1;%rand(dim-2, dim-2);
    space = 1/dim;
    f_loc = [0; 0; space; space];
    u_loc = [0; 0; space; space];

    save 'MyInput.mat' '-v4' ;
    s1 = './mgv';
    command = char(strcat(s1, {' '}, num2str(NList), {' '}, num2str(na),{' '}, num2str(nb), {' '}, num2str(omega), {' 1'}));
    system(command);
%    !./mgv $NList $na $nb $omega $i
    load('OutputAll.mat');
    plot(r);
end