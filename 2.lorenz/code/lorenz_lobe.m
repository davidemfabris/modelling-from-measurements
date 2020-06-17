clearvars -except input output; close all; clc

%% 1.Train NN
% Parameters Initialization
sigma = 10; rho=28; beta=8/3;

% Domain Definition
dt = 0.01;
tspan=dt:dt:10;

% Dataset Population
options = odeset('RelTol',1e-10,'AbsTol',1e-11);
input=[]; output=[];
data = 300;
alpha = 1;
tau = 1;

if size(input,1)==0 || size(output,1)==0
    for jj=1:data
        x0_nn=30*(rand(3,1)-0.5);
        [t,y] = ode45(@(t,x) lorenz_rhs(t,x, [], sigma,beta,rho),tspan,x0_nn,options);
        lobe=double(y(:,2)>=-y(:,1));
        for win=1:length(y)-tau-alpha
            input=[input; [y(win:win+alpha-1,:)]];
            output=[output; lobe(win+alpha-1+tau)];
        end
    end
end

%% Neural Network Definition
layers = [16 16 16];
net = feedforwardnet(layers);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas';
net.layers{3}.transferFcn = 'purelin';

%% Training
size(input.');
net.trainFcn = 'trainscg';
net = train(net,input.',output.','useGPU','yes','showResources','yes');

%% Validation
x0_ode=30*(rand(3,1)-0.5);
options = odeset('RelTol',1e-10,'AbsTol',1e-11);
[t,yode]=ode45(@(t,x) lorenz_rhs(t,x, [], sigma,beta,rho),tspan,x0_ode,options);

for win=1:length(t)-alpha-tau
    x0_nn=[yode(win:win+alpha-1,:)];
    yyy(win)=net(x0_nn');
    ynn(win)=round(yyy(win))';
end

%% 
figure
subplot 211
plot(t, lobe, 'LineWidth', 2)
xlabel('time'), ylabel('lobe')
axis([0 max(t) -0.5 1.5])
legend('real')
subplot 212
plot(t(1+alpha+tau:end), ynn, 'LineWidth', 2)
xlabel('time'), ylabel('lobe')
axis([0 max(t) -0.5 1.5])
legend('forecast')