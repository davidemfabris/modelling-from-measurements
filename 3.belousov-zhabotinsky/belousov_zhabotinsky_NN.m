clearvars -except BZ_tensor; close all; clc;

% load ('C:\Users\davidemaria.fabris\OneDrive - Politecnico di Milano\PhD\PhD Courses\Modelling from Measurements\BZ.mat');

[m,n,k] = size(BZ_tensor);

% for j=1:1
%     A=BZ_tensor(:,:,j);
%     pcolor(A), shading interp, pause(0.05)
% end

%% Data
input_tensor = permute(BZ_tensor(:,:,1:end-1),[3, 1, 2]);
output_tensor = permute(BZ_tensor(:,:,2:end),[3, 1, 2]);

II =  1:30:size(BZ_tensor, 3)-1;
jj = 1;
for ii = II
    input_mat = input_tensor(ii, :, :);
    input(:,jj) = input_mat(:);
    output_mat = output_tensor(ii, :, :);
    output(:,jj) = output_mat(:);
    jj = jj+1;
end

%% Network
layers = [16 16 16];
net = feedforwardnet(layers);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas';
net.layers{3}.transferFcn = 'purelin';
net.performParam.normalization = 'standard';
net.performParam.regularization = 0.5;
net.trainParam.max_fail = 50;

%% Training
net.trainFcn = 'traingdm';
net = train(net,input,output,'useGPU','yes','showResources','yes');

%% Testing
in = input(:,1);
toshow = reshape(in, 351, 451);
pcolor(toshow), shading interp, pause(0.5)
for dt = 1:10
    out = net(in);
    toshow = reshape(out, 351, 451);
    pcolor(toshow), shading interp, pause(0.5)
    in = out;
end