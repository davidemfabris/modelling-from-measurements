% clearvars -except net;
clear; close all; clc

%% 1.Train NN
load('data.mat');
instants = 10;%floor(size(data.input, 1)/1);
locations = 1;
input = data.input(1:instants, 1:locations^(-1):end);
output = data.output(1:instants,  1:locations^(-1):end);
time = data.time(1:instants);
space = data.space(1:locations^(-1):end);

%% Neural Network Definition
layers = [5 5 5];
net = feedforwardnet(layers);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas';
net.layers{3}.transferFcn = 'purelin';
net.trainParam.epochs = 1000;

%% Training
size(input.')
net.trainFcn = 'trainscg';
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
pcolor(space,time,ynn),shading interp, colormap(hot)