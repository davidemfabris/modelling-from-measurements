clear; close all; clc;
warning off;

%% Data
augment = 1;
Years = (1845:2:1903)';
Time = Years-Years(1);
dt = 2;
load('X', 'X')
if augment==1
    load X_augmented
    X = X_augmented;
    Years = (1845:1:1903)';
    Time = Years-Years(1);
    dt = 1;
end

% Plotting
figure
subplot 211
bar(Years, X(1,:), 'FaceColor', [0    0.4470    0.7410])
title('Time Series Data')
legend('Data Prey')
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
subplot 212
bar(Years, X(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
xlabel('Time [y]')
ylabel('Population [#]')
legend('Data Predator')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on

%% 1.DMD
clearvars -except augment X X1 X2 Years Time dt;
close all

% Data Arrangement
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);
x_0 = X1(:,1);

% DMD Algorithm
% SVD of X1
[U, Sigma, V]=svd(X1, 'econ');
r=2;
Ur=U(:,1:r);
Sigmar=Sigma(1:r,1:r);
Vr=V(:,1:r);

% Full order A
A = X2*Vr/Sigmar*Ur';
X2_pred = A*X1;

figure
subplot 211
plot(Years(2:end), X2(1,:), 'Color', [0    0.4470    0.7410], 'LineWidth', 2)
hold on
plot(Years(2:end), X2_pred(1,:), '--', 'Color', [0    0.4470    0.7410], 'LineWidth', 2)
title('Full Order Matrix Prediction')
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
hold off
legend('Data Prey', 'Predicted Prey')
subplot 212
plot(Years(2:end), X2(2,:), 'Color', [ 0.8500    0.3250    0.098], 'LineWidth', 2)
hold on
plot(Years(2:end), X2_pred(2,:), '--', 'Color', [ 0.8500    0.3250    0.098], 'LineWidth', 2)
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
legend('Data Predator', 'Predicted Predator')
hold off

% figure
% subplot 211, plot(Ur,'o','Linewidth',2), title('POD Modes - Elements of Ur')
% subplot 212, plot(Years(1:end-1),Vr,'Linewidth',2), title('Time Modes - Elements of Vr')
% xlabel('time [y]'), ylabel('population [#]')

% Reduced order A
Atilde = Ur'*X2*Vr/Sigmar;
X2_pred_red = Atilde*X1;

figure
subplot 211
plot(Years(2:end), X2(1,:), 'Color', [0    0.4470    0.7410], 'LineWidth', 2)
hold on
plot(Years(2:end), X2_pred_red(1,:), '--', 'Color', [0    0.4470    0.7410], 'LineWidth', 2)
title('Reduced Order Matrix Predition')
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
hold off
legend('Data Prey', 'Predicted Prey')
subplot 212
plot(Years(2:end), X2(2,:), 'Color', [ 0.8500    0.3250    0.098], 'LineWidth', 2)
hold on
plot(Years(2:end), X2_pred_red(2,:), '--', 'Color', [ 0.8500    0.3250    0.098], 'LineWidth', 2)
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
legend('Data Predator', 'Predicted Predator')
hold off

% Spectral decomposition of Atilde
[W,Lambda] = eig(Atilde);

% Spatial modes of A
Phi=X2*(Vr/Sigmar)*W;

% Frequencies
mu=diag(Lambda);
omega=log(mu)/dt;

% Initial Conditions
b = Phi\x_0;

% Reconstruction
x_dmd = zeros(r,length(Time));
for instant = 1:length(Time)
    x_dmd(:,instant) =b.*exp(omega*(Time(instant)));
end
x_rec = Phi*x_dmd;

% figure
% subplot(2,2,1)
% waterfall(1:2,Years,X.')
% colormap([0 0 0])
% subplot(2,2,2)
% waterfall(1:2,1:r,abs(Phi).')
% colormap([0 0 0])
% xlabel('x','Fontsize',14)
% ylabel('modes','Fontsize',14)
% subplot(2,2,3)
% plot(diag(Sigma)/sum(diag(Sigma)),'ko')
% subplot(2,2,4)
% waterfall(1:2,Years,abs(x_rec).')
% colormap([0 0 0])
% xlabel('x','Fontsize',14)
% ylabel('time','Fontsize',14)
% text(-40,8,3,'|u|','Fontsize',14)
% text(-20,3,7,'DMD','Fontsize',14)

% Plotting
figure
plot(Time, x_rec(1, :).', 'Color', [0    0.4470    0.7410], 'Linewidth', 2)
grid on
hold on
plot(Time, x_rec(2, :).', 'Color', [0.8500    0.3250    0.098], 'Linewidth', 2)
grid on
title('Dynamic Mode Decomposition')
legend('Prey - Snowshoe Hare', 'Predator - Canadian Lynx')
xlabel('Time [y]')
ylabel('Population [#]')

% 1.Koopman
% clearvars -except X X1 X2 Years Time dt;
% close all
%
% XK1 = [X1; X1(1,:).*X1(2,:)];
% XK2 = [X2; X2(1,:).*X2(2,:)];
% x_0K = XK1(:,1);
%
% % Koopman Algorithm
% % SVD of X1
% [UK, SigmaK, VK]=svd(XK1, 'econ');
% rK=1;
% UKr=UK(:,1:rK);
% SigmaKr=SigmaK(1:rK,1:rK);
% VKr=VK(:,1:rK);
%
% % figure
% % subplot 211, plot(UKr,'o','Linewidth',2), title('POD Modes - Elements of UKr')
% % subplot 212, plot(Years(1:end-1),VKr,'Linewidth',2), title('Time Modes - Elements of VKr')
% % xlabel('time [y]'), ylabel('population [#]')
%
% % Reduced order A
% AKtilde = UKr'*XK2*VKr/SigmaKr;
% X2_pred_K = AKtilde*X1;
%
% figure
% subplot 211
% plot(Years(2:end), X2(1,:), 'Color', [0    0.4470    0.7410], 'LineWidth', 2)
% hold on
% plot(Years(2:end), X2_pred_K(1,:), '--', 'Color', [0    0.4470    0.7410], 'LineWidth', 2)
% title('Koopman')
% xlabel('Time [y]')
% ylabel('Population [#]')
% axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
% grid on
% hold off
% legend('Real Prey', 'Predicted Prey')
% subplot 212
% plot(Years(2:end), X2(2,:), 'Color', [ 0.8500    0.3250    0.098], 'LineWidth', 2)
% hold on
% plot(Years(2:end), X2_pred_K(2,:), '--', 'Color', [ 0.8500    0.3250    0.098], 'LineWidth', 2)
% xlabel('Time [y]')
% ylabel('Population [#]')
% axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
% grid on
% legend('Real Predator', 'Predicted Predator')
% hold off
%
% % Spectral decomposition of Atilde
% [WK,LambdaK] = eig(AKtilde);
%
% % Spatial modes of A
% PhiK=XK2*(VKr/SigmaKr)*WK;
%
% % Frequencies
% muK=diag(LambdaK);
% omegaK=log(muK)/dt;
%
% % Initial Conditions
% bK = PhiK\x_0K;
%
% % Reconstruction
% x_dmdK = zeros(rK,length(Time));
% for instant = 1:length(Time)
%     x_dmdK(:,instant) =bK.*exp(omegaK*(Time(instant)));
% end
% x_recK = PhiK*x_dmdK;

% Plotting
% figure
% plot(Time, x_recK(1, :).', 'Color', [0    0.4470    0.7410], 'Linewidth', 2)
% grid on
% hold on
% plot(Time, x_recK(2, :).', 'Color', [0.8500    0.3250    0.098], 'Linewidth', 2)
% grid on
% title('DMD - Koopman')
% legend('hare', 'lynx')
% xlabel('Time [years]')
% ylabel('Mode [#]')

% 1.Least Squares Fit
% X1 = X(:, 1:end-1);
% X2 = X(:, 2:end);
%
% A_pinv = X2*pinv(X1);
% X2_pinv = A_pinv*X1;
% X3_pinv = A_pinv*X2;
%
% figure
% subplot 211
% plot(Years(2:end), X2_pinv(1,:), 'Color', [0    0.4470    0.7410], 'LineWidth', 2)
% hold on
% plot(Years(2:end), X2(1,:), '--', 'Color', [0    0.4470    0.7410], 'LineWidth', 2)
% title('Least Squares')
% xlabel('Time [y]')
% ylabel('Population [#]')
% axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
% grid on
% hold off
% legend('Real Prey', 'Predicted Prey')
% subplot 212
% plot(Years(2:end), X2_pinv(2,:), 'Color', [ 0.8500    0.3250    0.098], 'LineWidth', 2)
% hold on
% plot(Years(2:end), X2(2,:), '--', 'Color', [ 0.8500    0.3250    0.098], 'LineWidth', 2)
% xlabel('Time [y]')
% ylabel('Population [#]')
% axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
% grid on
% legend('Real Predator', 'Predicted Predator')
% hold off

%% 2.Time-Delay Embedding
clearvars -except augment X X1 X2 Years Time dt;
close all

% Time Delay
l=length(Time);
if augment == 1
    mc=40;
else
    mc=20;
end

mo=l-mc;
H=[];
for shift=1:mo
    H = [H;X(:,shift:shift+mc-1)];
end
XH1 = H(:, 1:end-1);
XH2 = H(:, 2:end);

% Singular Value Decomposition
[UH, SigmaH, VH] = svd(XH1, 'econ');

% Plotting
figure
semilogy(diag(SigmaH)/sum(diag(SigmaH)), 'o', 'Linewidth', 2)
grid on
title('Normalized Singular Values')
xlabel('Mode Number [#]')
ylabel('Normalized Value [-]')

if augment == 1
    rH = 23;%size(SigmaH, 1);
else
    rH = 11;%size(SigmaH, 1);
end

figure
subplot 211, plot(UH(:,1:rH),'Linewidth',2), title('POD Modes - Elements of U'), grid on
xlabel('time [y]'), ylabel('population [#]'), legend('1', '2', '3')
subplot 212, plot(Time(1:mc-1),VH(:,1:rH),'Linewidth',2), title('POD Modes - Elements of V'), grid on
xlabel('time [y]'), ylabel('population [#]'), legend('1', '2', '3')

UHr = UH(:,1:rH);
VHr = VH(:,1:rH);
SigmaHr = SigmaH(1:rH, 1:rH);

% Reduced order A
AHtilde = UHr'*XH2*VHr/SigmaHr;

% Spectral decomposition of Atilde
[WH,LambdaH] = eig(AHtilde);

% Spatial modes of A
PhiH=XH2*(VHr/SigmaHr)*WH;

% Frequencies
muH=diag(LambdaH);
omegaH=log(muH)/dt;

figure
plot(real(omegaH), imag(omegaH), 'o')
title('Omega')
xlabel('Real [-]'), ylabel('Imaginary [-]'), grid on

% Initial Conditions
x_0 = H(:,1);
bH = PhiH\x_0;

% Reconstruction
xH_dmd = zeros(rH,length(Time));
for instant = 1:length(Time)
    xH_dmd(:,instant) =bH.*exp(omegaH*(Time(instant)));
end
xH_rec = PhiH*xH_dmd;

% Plotting
figure
subplot 211
bar(Years, X(1,:), 'FaceColor', [0    0.4470    0.7410])
title('Time Series Data')
legend('Data Prey')
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
subplot 212
bar(Years, X(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
xlabel('Time [y]')
ylabel('Population [#]')
legend('Data Predator')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on

% Plotting
figure
subplot 211
bar(Years, xH_rec(1, :).', 'FaceColor', [0    0.4470    0.7410])
grid on
title('Time-Delay Embedding Reconstruction')
legend('Reconstructed Prey')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
xlabel('Time [y]')
ylabel('Population [#]')
subplot 212
bar(Years, xH_rec(2, :).', 'FaceColor', [0.8500    0.3250    0.098])
grid on
legend('Reconstructed Predator')
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])

%% 3.Lotka-Volterra
clearvars -except X X1 X2 Years Time dt;
close all

% Differenciating
Xdot = X*0;
for ii=3:size(X,2)-3
    Xdot(:, ii) = 1/(12*dt)*(X(:, ii-2)-8*X(:, ii-1)+8*X(:, ii+1)-X(:, ii+2));
end
Xdot=Xdot(:,3:end-3);

% Fitting with pinv
Coeff=[];
Coeff1 = Xdot(1,:)*pinv([X(1,3:end-3);X(1,3:end-3).*X(2,3:end-3)]);
Coeff2 = Xdot(2,:)*pinv([X(2,3:end-3);X(1,3:end-3).*X(2,3:end-3)]);
Coeff(1, 1:2:3) = Coeff1;
Coeff(2, 2:3) = Coeff2;

% Fitting with lasso
B=[];
XX1 = [X(1,3:end-3); X(1,3:end-3).*X(2,3:end-3)]';
YY1 = Xdot(1,:);
XX2 = [X(2,3:end-3); X(1,3:end-3).*X(2,3:end-3)]';
YY2 = Xdot(2,:);
lambda=0:0.01:0.5;
[B1,STATS1] = lasso(XX1, YY1, 'Lambda', lambda);
[B2,STATS2] = lasso(XX2, YY2, 'Lambda', lambda);
B(1, 1:2:3, :) = B1(:,:);
B(2, 2:3, :) = B2(:, :);

% Plotting
% figure
% subplot 211
% bar(Xdot(1,:), 'FaceColor', [0    0.4470    0.7410])
% title('hare')
% xlabel('time [years]')
% ylabel('growth [#/y]')
% subplot 212
% bar(Xdot(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
% title('lynx')
% xlabel('time [years]')
% ylabel('growth [#/y]')

% ode45 solution using pinv coefficients
x_0 = X1(:,1);
t=0:dt:Time(end);
[time,sol_pinv]=ode45('lotka_volterra_rhs',t,x_0,[],Coeff);
sol_pinv=sol_pinv';

% ode45 solution using lasso coefficients
for ii=1:length(lambda)
    [time,sol_lasso(:,:,ii)]=ode45('lotka_volterra_rhs',t,x_0,[],B(:, :, ii));
end
sol_lasso = permute(sol_lasso, [2 1 3]);

% minimum discrepancies lambdas for lasso
rows = 1;
columns = 2;
error = sqrt(sum((sol_lasso-X).^2, columns));
[val1, idx1] = min(error(1, :));
[val2, idx2] = min(error(2, :));
B_min = [B(1, :, idx1);B(2, :, idx2)];

string = {'b', 'p', 'r', 'd'};
figure
plot(Coeff(Coeff~=0), 'b*')
hold on
plot(B_min(B_min~=0), 'ro')
legend('Pseudo-Inverse', 'Lasso')
text([1:length(Coeff(Coeff~=0))]+0.05, Coeff(Coeff~=0), string)
text([1:length(B_min(B_min~=0))]+0.05, B_min(B_min~=0), string)
xlabel('Coefficient [#]'), ylabel('Value [-]')
grid on
hold off
title('Lotka-Volterra Coefficients')

if size(error, 3)>1
    figure
    for ii = 1:size(error,3)
        plot(lambda(ii), error(1, :, ii), 'ro')
        hold on
        plot(lambda(ii), error(2, :, ii), 'bo')
    end
    hold off
else
    plot(lambda, error, 'o')
end

% pinv coefficients ode45 solutions plotting
figure
subplot 211
bar(Years, sol_pinv(1,:), 'FaceColor', [0    0.4470    0.7410]), grid on
title('Snowshoe Hare - Pseudo Inverse')
xlabel('Time [y]')
ylabel('Population [#]')
subplot 212
bar(Years, sol_pinv(2,:), 'FaceColor', [ 0.8500    0.3250    0.098]), grid on
title('Canadian Lynx - Pseudo Inverse')
xlabel('Time [y]')
ylabel('Population [#]')

% lasso coefficients ode45 solutions plotting
figure
subplot 211
bar(Years, sol_lasso(1,:,idx1), 'FaceColor', [0    0.4470    0.7410]), grid on
title('Snowshoe Hare - Lasso')
xlabel('Time [y]')
ylabel('Population [#]')
subplot 212
bar(Years, sol_lasso(2,:,idx2), 'FaceColor', [ 0.8500    0.3250    0.098]), grid on
title('Canadian Lynx - Lasso')
xlabel('Time [y]')
ylabel('Population [#]')

%% 4.SINDy
% clearvars -except X Years Time dt t;
% close all
% 
% % Data Organisation
% x1=X(1, :);
% x2=X(2, :);
% x_0 = X(:,1);
% 
% % Differenciating
% x1dot = x1*0;
% x2dot = x2*0;
% 
% for ii=3:size(X,2)-3
%     x1dot(:, ii) = 1/(12*dt)*(x1(ii-2)-8*x1(ii-1)+8*x1(ii+1)-x1(ii+2));
%     x2dot(:, ii) = 1/(12*dt)*(x2(ii-2)-8*x2(ii-1)+8*x2(ii+1)-x2(ii+2));
% end
% x1s = x1(3:end-3)';
% x1dot = x1dot(3:end-3)';
% x2s=x2(3:end-3)';
% x2dot = x2dot(3:end-3)';
% 
% A_handle=@(x1s, x2s) [x1s x2s x1s.^2 x1s.*x2s x2s.^2 x1s.^3 (x2s.^2).*x1s...
%                       x2s.^3 cos(x1s) cos(x2s) sin(x1s) sin(x2s) cos(x1s).*cos(x2s)...
%                       cos(x1s).*sin(x1s) cos(x2s).*sin(x2s) cos(x1s).*sin(x2s) cos(x2s).*sin(x1s)];
%                   
% A = A_handle(x1s, x2s);
% 
% % xi1=A\x1dot;
% % xi2=A\x2dot;
% % xi1=pinv(A)*x1dot;
% % xi2=pinv(A)*x2dot;
% xi1=lasso(A,x1dot.','Lambda',0.1);
% xi2=lasso(A,x2dot.','Lambda',0.1);
% 
% figure
% subplot(2,1,1), bar(xi1, 'FaceColor', [0    0.4470    0.7410]), grid on
% subplot(2,1,2), bar(xi2, 'FaceColor', [ 0.8500    0.3250    0.098]), grid on
% 
% indexes_1 = abs(xi1)<0;
% indexes_2 = abs(xi2)<0;
% 
% xi1(indexes_1) = 0;
% xi2(indexes_2) = 0;
% 
% param = [xi1'; xi2'];
% 
% [time,sol_sindy]=ode45('lotka_volterra_sindy_rhs',t,x_0,[], A_handle, param);
% sol_sindy = sol_sindy';
% 
% figure
% subplot 211
% bar(sol_sindy(1,:), 'FaceColor', [0    0.4470    0.7410])
% title('Snowshoe Hare - SINDy')
% xlabel('Time [y]')
% ylabel('Population [#]')
% subplot 212
% bar(sol_sindy(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
% title('Canadian Lynx - SINDy')
% xlabel('Time [y]')
% ylabel('Population [#]')

%% 4.SINDy - Discrete Time
clearvars -except A_handle X X1 X2 Years Time dt t;
close all

Xk = X(:,1:end-1)';
Xk1 = X(:,2:end)';

A_handle=@(x1s, x2s) [x1s x2s x1s.^2 x1s.*x2s x2s.^2 x1s.^3 (x2s.^2).*x1s...
                      x2s.^3 cos(x1s) cos(x2s) sin(x1s) sin(x2s) cos(x1s).*cos(x2s)...
                      cos(x1s).*sin(x1s) cos(x2s).*sin(x2s) cos(x1s).*sin(x2s) cos(x2s).*sin(x1s)];

A_dt=A_handle(Xk(:,1), Xk(:,2));

xi_dt = pinv(A_dt)*Xk1;
xi1_dt_pinv = xi_dt(:,1);
xi2_dt_pinv = xi_dt(:,2);
xi1_dt_lasso = lasso(A_dt, Xk1(:,1), 'Lambda', 2);
xi2_dt_lasso = lasso(A_dt, Xk1(:,2), 'Lambda', 2);

figure
subplot(2,1,1), bar(xi1_dt_pinv, 'FaceColor', [0    0.4470    0.7410]), grid on, title('Pseudo Inverse')
subplot(2,1,2), bar(xi2_dt_pinv, 'FaceColor', [0.8500    0.3250    0.098]), grid on

figure
subplot(2,1,1), bar(xi1_dt_lasso, 'FaceColor', [0    0.4470    0.7410]), grid on, title('Lasso')
subplot(2,1,2), bar(xi2_dt_lasso, 'FaceColor', [0.8500    0.3250    0.098]), grid on

threshold = 0;
indexes_1_pinv = abs(xi1_dt_pinv)<threshold;
indexes_2_pinv = abs(xi2_dt_pinv)<threshold;
indexes_1_lasso = abs(xi1_dt_lasso)<threshold;
indexes_2_lasso = abs(xi2_dt_lasso)<threshold;

xi1_dt_pinv(indexes_1_pinv)=0;
xi2_dt_pinv(indexes_2_pinv)=0;
xi1_dt_lasso(indexes_1_lasso)=0;
xi2_dt_lasso(indexes_2_lasso)=0;

Xpred_pinv = [A_dt*xi1_dt_pinv, A_dt*xi2_dt_pinv];
Xpred_lasso = [A_dt*xi1_dt_lasso, A_dt*xi2_dt_lasso];

figure
subplot 211
bar(Xk1(:,1), 'FaceColor', [0    0.4470    0.7410])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Prey')
title('Pseudo Inverse')
subplot 212
bar(Xpred_pinv(:,1), 'FaceColor', [0    0.4470    0.7410])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Predicted Prey')
figure
subplot 211
bar(Xk1(:,2), 'FaceColor', [0.8500    0.3250    0.098])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Predator')
title('Pseudo Inverse')
subplot 212
bar(Xpred_pinv(:,2),  'FaceColor', [0.8500    0.3250    0.098])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Predicted Predator')

figure
subplot 211
bar(Xk1(:,1), 'FaceColor', [0    0.4470    0.7410])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Prey')
title('Lasso')
subplot 212
bar(Xpred_lasso(:,1), 'FaceColor', [0    0.4470    0.7410])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Predicted Prey')
figure
subplot 211
bar(Xk1(:,2), 'FaceColor', [0.8500    0.3250    0.098])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Predator')
title('Lasso')
subplot 212
bar(Xpred_lasso(:,2),  'FaceColor', [0.8500    0.3250    0.098])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Predicted Predator')