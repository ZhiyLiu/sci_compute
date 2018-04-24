clear all;

% generate grid data
dim = 1025;
space = 1/dim;
u = zeros(dim, dim);
f = randomRHS(dim);
%f = smoothRHS(dim);
f_loc = [0; 0; space; space];
u_loc = [0; 0; space; space];

save 'MyInput.mat' '-v4' ;

% vary omega
omegaList = [0.2,0.66,0.8, 1, 0];
NList = length(omegaList);
na = 3; nb = 3;
for i = 1:NList
    omega = omegaList(i);

    s1 = './mgv';
    command = char(strcat(s1, {' '}, num2str(NList), {' '}, num2str(na),{' '}, num2str(nb), {' '}, num2str(omega), {' '}, num2str(i)));
    system(command);
%    !./mgv $NList $na $nb $omega $i
    load('OutputAll.mat');
end

plot(log10(r1));hold on;
plot(log10(r2));hold on;
plot(log10(r3));hold on;
plot(log10(r4));hold on;
plot(log10(r5));hold on;
maxLength = max([length(r1), length(r2), length(r3), length(r4), length(r5)]);
set(gca, 'xtick', 0:maxLength+10);
%ylim([-1,-0.2]);
legend( ...
num2str(omegaList(1)),...
num2str(omegaList(2)), ...
num2str(omegaList(3)), ...
num2str(omegaList(4)), ...
num2str(omegaList(5))...
);