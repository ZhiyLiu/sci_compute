clear all;
dim = 20;
space = 1/dim;
u = zeros(dim, dim);
%f = randomRHS(dim);
f = smoothRHS(dim);
f_loc = [0; 0; space; space];
u_loc = [0; 0; space; space];

save 'MyInput.mat' '-v4' ;

%vary sweep number
sweepList = [1,2;2,1;3,3;5,5;7,7];
NList = size(sweepList, 1);
omega = 0.8;
for i = 1:NList

    nb = sweepList(i,1);
    na = sweepList(i,2);
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
char(strcat(num2str(sweepList(1,1)),{','}, num2str(sweepList(1,2)))),...
char(strcat(num2str(sweepList(2,1)),{','}, num2str(sweepList(2,2)))),...
char(strcat(num2str(sweepList(3,1)),{','}, num2str(sweepList(3,2)))),...
char(strcat(num2str(sweepList(4,1)),{','}, num2str(sweepList(4,2)))),...
char(strcat(num2str(sweepList(5,1)),{','}, num2str(sweepList(5,2))))...
);