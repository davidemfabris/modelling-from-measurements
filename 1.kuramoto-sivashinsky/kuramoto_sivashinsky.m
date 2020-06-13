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
layers = [32 32 32];
net = feedforwardnet(layers);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas';
net.layers{3}.transferFcn = 'purelin';
net.trainParam.epochs = 1000;

%% Training
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
pcolor(space,time,ynn),shading interp, colormap(hot)