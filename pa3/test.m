close all;
clear all;
load Input.mat
load Output.mat
% Draw with something like
figure; plot(x,y); hold on;plot(xy(1,:),xy(2,:),'o'); 
 h =0.01;
xEval = x(2:end-1); % strip out the first and last
yD = (y(3:end)-y(1:end-2))/(2*h);
yDD = (y(3:end)-2*y(2:end-1)+y(1:end-2))/(h*h);
figure; plot(xEval, yD); hold on; plot(xEval, yDD);

% load Input.mat
% pp = csape(xy(1, :), xy(2,:), 'variational');