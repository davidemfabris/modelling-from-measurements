% clearvars -except net;
clear; close all; clc

%% 1.Train NN
load('data.mat');
instants = floor(size(data.input, 1)/1);
locations = 1;
input = data.input(1:instants, 1:locations^(-1):end);
output = data.output(1:instants,  1:locations^(-1):end);
time = data.time(1:instants);
space = data.space(1:locations^(-1):end);

%% Neural Network Definition
layers = [32 32 32 32];
net = feedforwardnet(layers);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'logsig';
net.layers{3}.transferFcn = 'logsig';
net.trainFcn = 'trainscg';
net.trainParam.epochs = 1000;
net.trainParam.max_fail = 1000;

%% Training
size(input.')
net = train(net,input.',output.','useGPU','yes','showResources','yes');

%%
x0 = input(1, :);
ynn(1,:)= x0;
for jj=2:size(input,1)
    y0=net(x0.');
    ynn(jj,:)=y0.'; x0=y0.';
end

%%
figure
surf(space,time,ynn),shading interp, colormap(hot)
title('Flame Front Evolution'), xlabel('width [m]'), ylabel('time [s]'), zlabel('flame front')
figure
surf(space,time,output),shading interp, colormap(hot)
title('Flame Front Evolution'), xlabel('width [m]'), ylabel('time [s]'), zlabel('flame front')