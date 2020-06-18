clear; close all; clc

%% 1.Train NN
% Parameters Initialization
sigma = 10; rho=28; beta=8/3;

% Domain Definition
dt = 0.01;
tspan=dt:dt:10;

% Dataset Population
options = odeset('RelTol',1e-10,'AbsTol',1e-11);
input=[]; output=[];
data = 20;
alpha = 5;
tau = 10;

if size(input,1)==0 || size(output,1)==0
    for jj=1:data
        x0_nn=30*(rand(3,1)-0.5);
        [t,y] = ode45(@(t,x) lorenz_rhs(t,x, [], sigma,beta,rho),tspan,x0_nn,options);
        lobe=double(y(:,2)>=-y(:,1));
        for win=1:length(y)-tau-alpha+1
            this = [y(win:win+alpha-1,:), lobe(win:win+alpha-1)];
            input=[input; reshape(this', [1,numel(this)])];
            output=[output; lobe(win+alpha+tau-1)];
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
net.trainParam.max_fail=20;
net = train(net,input.',output.','useGPU','yes','showResources','yes');

%% Validation
x0_ode=30*(rand(3,1)-0.5);
options = odeset('RelTol',1e-10,'AbsTol',1e-11);
[t,yode]=ode45(@(t,x) lorenz_rhs(t,x, [], sigma,beta,rho),tspan,x0_ode,options);
lobe=double(yode(:,2)>=-yode(:,1));

for win=1:length(t)-alpha-tau
    this = [yode(win:win+alpha-1,:), lobe(win:win+alpha-1)];
    x0_nn=[reshape(this', [1,numel(this)])];
    yyy(win)=net(x0_nn');
    ynn(win)=double(yyy(win)>=0.5)';
end

l = length(lobe(1+alpha+tau:end));

figure
plot3(yode(:,1), yode(:,2), yode(:,3))

figure
subplot 311
plot(dt:dt:l*dt, yode(1+alpha+tau:end, 1), 'LineWidth', 2)
xlabel('time'), ylabel('x')
legend('x')
subplot 312
plot(dt:dt:l*dt, yode(1+alpha+tau:end, 2), 'LineWidth', 2)
xlabel('time'), ylabel('y')
legend('y')
subplot 313
plot(dt:dt:l*dt, yode(1+alpha+tau:end, 3), 'LineWidth', 2)
xlabel('time'), ylabel('z')
legend('z')

figure
subplot 211
plot(dt:dt:l*dt, lobe(1+alpha+tau:end), 'LineWidth', 2)
xlabel('time'), ylabel('lobe')
legend('real')
subplot 212
plot(dt:dt:l*dt, ynn, 'LineWidth', 2)
xlabel('time'), ylabel('lobe')
legend('forecast')