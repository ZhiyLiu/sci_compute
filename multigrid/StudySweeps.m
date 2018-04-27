clear all;
close all;
dim = 513;
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
    command = char(strcat(s1, {' '}, num2str(NList), {' '}, num2str(nb),{' '}, num2str(na), {' '}, num2str(omega), {' '}, num2str(i), {' 25'}));
    system(command);
%    !./mgv $NList $na $nb $omega $i
    load('OutputAll.mat');
end
time1_list = zeros(length(r1),1);
for i = 1: length(r1)
    time1_list(i) = i * time1;
end
time2_list = zeros(length(r2),1);
for i = 1: length(r2)
    time2_list(i) = i * time2;
end
time3_list = zeros(length(r3),1);
for i = 1: length(r3)
    time3_list(i) = i * time3;
end
time4_list = zeros(length(r4),1);
for i = 1: length(r4)
    time4_list(i) = i * time4;
end
time5_list = zeros(length(r5),1);
for i = 1: length(r5)
    time5_list(i) = i * time5;
end
% figure scaled by time
subplot(3,1,1);
plot(time1_list, log10(r1));hold on;
plot(time2_list, log10(r2));hold on;
plot(time3_list, log10(r3));hold on;
plot(time4_list, log10(r4));hold on;
plot(time5_list, log10(r5));hold on;

title('residual with time');
ylabel('log10(residual_norm)', 'Interpreter','none');
xlabel('time(s)');
maxLength = max([length(r1), length(r2), length(r3), length(r4), length(r5)]);
%set(gca, 'xtick', 0:maxLength+10);
%ylim([-1,-0.2]);
legend( ...
char(strcat(num2str(sweepList(1,1)),{','}, num2str(sweepList(1,2)))),...
char(strcat(num2str(sweepList(2,1)),{','}, num2str(sweepList(2,2)))),...
char(strcat(num2str(sweepList(3,1)),{','}, num2str(sweepList(3,2)))),...
char(strcat(num2str(sweepList(4,1)),{','}, num2str(sweepList(4,2)))),...
char(strcat(num2str(sweepList(5,1)),{','}, num2str(sweepList(5,2))))...
);
%figure scaled by sweep number
subplot(3,1,2);
x_base = 1:length(r1);
plot((sweepList(1,1) + sweepList(1,2) )* x_base, log10(r1));hold on;
plot((sweepList(2,1) + sweepList(2,2) )* x_base, log10(r2));hold on;
plot((sweepList(3,1) + sweepList(3,2) )* x_base, log10(r3));hold on;
plot((sweepList(4,1) + sweepList(4,2) )* x_base, log10(r4));hold on;
plot((sweepList(5,1) + sweepList(5,2) )* x_base, log10(r5));hold on;
maxLength = max([length(r1), length(r2), length(r3), length(r4), length(r5)]);
%set(gca, 'xtick', 0:maxLength+10);
%ylim([-1,-0.2]);
legend( ...
char(strcat(num2str(sweepList(1,1)),{','}, num2str(sweepList(1,2)))),...
char(strcat(num2str(sweepList(2,1)),{','}, num2str(sweepList(2,2)))),...
char(strcat(num2str(sweepList(3,1)),{','}, num2str(sweepList(3,2)))),...
char(strcat(num2str(sweepList(4,1)),{','}, num2str(sweepList(4,2)))),...
char(strcat(num2str(sweepList(5,1)),{','}, num2str(sweepList(5,2))))...
);

title('scaled by sweep number');
ylabel('log10(residual_norm)', 'Interpreter','none');
xlabel('sweep counts');
subplot(3,1,3);
plot(log10(r1));hold on;
plot(log10(r2));hold on;
plot(log10(r3));hold on;
plot(log10(r4));hold on;
plot(log10(r5));hold on;
%set(gca, 'xtick', 0:maxLength+1);
legend( ...
char(strcat(num2str(sweepList(1,1)),{','}, num2str(sweepList(1,2)))),...
char(strcat(num2str(sweepList(2,1)),{','}, num2str(sweepList(2,2)))),...
char(strcat(num2str(sweepList(3,1)),{','}, num2str(sweepList(3,2)))),...
char(strcat(num2str(sweepList(4,1)),{','}, num2str(sweepList(4,2)))),...
char(strcat(num2str(sweepList(5,1)),{','}, num2str(sweepList(5,2))))...
);

title('residual with v_cycles', 'Interpreter','none');
ylabel('log10(residual_norm)', 'Interpreter','none');
xlabel('# v_cycles', 'Interpreter','none');