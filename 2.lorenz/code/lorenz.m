clearvars -except net; close all; clc

%% 1.Train NN
% Parameters Initialization
sigma = 10; rho=28; beta=8/3;

% Initial Condition
x0=[-10; 10; 0];

% Domain Definition
dt = 0.01;
tspan=dt:dt:50;

% Solution
options = odeset('RelTol',1e-10,'AbsTol',1e-11);
[t,x]=ode45(@(t,x) lorenz_rhs(t,x, [], sigma,beta,rho),tspan,x0,options);

% Plotting

% figure
plot3(x(:,1),x(:,2),x(:,3));
hold on
% plot3(x(1,1),x(1,2),x(1,3),'o');
% plot3(x(end,1),x(end,2),x(end,3),'ro');
hold off
title('Lorenz Attractor')
xlabel('x'), ylabel('y'), zlabel('z')
grid on

% Data Population
r_train = [10, 28, 40];
input=[]; output=[];
data = 100;
for jj=1:data
    x0=30*(rand(3,1)-0.5);
    ii = ceil(length(r_train)*rand(1));
    rho_train = r_train(ii);
    [t,y] = ode45(@(t,x) lorenz_rhs(t,x, [], sigma,beta,rho_train),tspan,x0,options);
    input=[input; y(1:end-1,:)];
    output=[output; y(2:end,:)];
end

%% Neural Network Definition
layers = [16 16 16];
net = feedforwardnet(layers);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas';
net.layers{3}.transferFcn = 'purelin';
net.trainParam.epochs = 1000;

%% Training
size(input.')
net.trainFcn = 'trainscg';
net = train(net,input.',output.','useGPU','yes','showResources','yes');

% Naming
net_filename = strcat('data_',num2str(data),'_rho_',num2str(r_train),'_layers_',num2str(layers),'_net');
input_filename = strcat('data_',num2str(data),'_rho_',num2str(r_train),'_layers_',num2str(layers),'_input');
output_filename = strcat('data_',num2str(data),'_rho_',num2str(r_train),'_layers_',num2str(layers),'_output');

% Saving
% save(input_filename, 'input')
% save(output_filename, 'output')
% save(net_filename, 'net')

%% Validation
x0=30*(rand(3,1)-0.5);
r_val = r_train;
for rho_val =r_val
    options = odeset('RelTol',1e-10,'AbsTol',1e-11);
    [t,yode]=ode45(@(t,x) lorenz_rhs(t,x, [], sigma,beta,rho_val),tspan,x0,options);
    ynn(1,:)=x0;
    for jj=2:length(t)
        y0=net(x0);
        ynn(jj,:)=y0.'; x0=y0;
    end
    
    figure
    plot3(yode(:,1),yode(:,2),yode(:,3),'Linewidth',2)
    title(strcat('Validation \rho=',num2str(rho_val)))
    hold on
    plot3(ynn(:,1),ynn(:,2),ynn(:,3),':','Linewidth',2)
    grid on
    xlabel('x'), ylabel('y'), zlabel('z')
    hold off
    
    figure
    subplot 311
    plot(t, yode(:,1), t, ynn(:,1))
    title(strcat('Validation \rho=',num2str(rho_val)))
    xlabel('time'), ylabel('x')
    legend('ode', 'nn')
    subplot 312
    plot(t, yode(:,2), t, ynn(:,2))
    xlabel('time'), ylabel('y')
    legend('ode', 'nn')
    subplot 313
    plot(t, yode(:,3), t, ynn(:,3))
    xlabel('time'), ylabel('z')
    legend('ode', 'nn')
    hold off
end


%% Test
x0=30*(rand(3,1)-0.5);
r_test = [17, 35];
for rho_test =r_test
    options = odeset('RelTol',1e-10,'AbsTol',1e-11);
    [t,yode]=ode45(@(t,x) lorenz_rhs(t,x, [], sigma,beta,rho_test),tspan,x0,options);
    ynn(1,:)=x0;
    for jj=2:length(t)
        y0=net(x0);
        ynn(jj,:)=y0.'; x0=y0;
    end
    
    figure
    plot3(yode(:,1),yode(:,2),yode(:,3),'Linewidth',2)
    title(strcat('Test \rho=',num2str(rho_test)))
    hold on
    plot3(ynn(:,1),ynn(:,2),ynn(:,3),':','Linewidth',2)
    grid on
    xlabel('x'), ylabel('y'), zlabel('z')
    hold off
    
    figure
    subplot 311
    plot(t, yode(:,1), t, ynn(:,1))
    title(strcat('Test \rho=',num2str(rho_test)))
    xlabel('time'), ylabel('x')
    legend('ode', 'nn')
    subplot 312
    plot(t, yode(:,2), t, ynn(:,2))
    xlabel('time'), ylabel('y')
    legend('ode', 'nn')
    subplot 313
    plot(t, yode(:,3), t, ynn(:,3))
    xlabel('time'), ylabel('z')
    legend('ode', 'nn')
    hold off
end