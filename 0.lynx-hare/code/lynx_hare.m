clear; close all; clc;
warning off;

%% Data Plotting
augment = 0;
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

%% 1.Dynamic Mode Decomposition
clearvars -except augment X X1 X2 Years Time dt;
close all

% Data Arrangement
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);
x_0 = X1(:,1);

% DMD Algorithm
% SVD of X1
[U, Sigma, V]=svd(X1, 'econ');

% Plotting
figure
semilogy(diag(Sigma)/sum(diag(Sigma)), 'o', 'Linewidth', 2)
grid on
title('Normalized Singular Values')
xlabel('Mode Number [#]')
ylabel('Normalized Value [-]')


r=2;
Ur=U(:,1:r);
Sigmar=Sigma(1:r,1:r);
Vr=V(:,1:r);

% Full order A
A = X2*Vr/Sigmar*Ur';
X2_pred = A*X1;

figure(70)
subplot 411
bar(Years(2:end), X2(1,:), 'FaceColor', [0    0.4470    0.7410])
title('Full Order Matrix Fitting')
xlabel('Time [y]')
ylabel('Population [#]')
grid on
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
legend('Real Prey')
subplot 412
bar(Years(2:end), X2_pred(1,:), 'FaceColor', [0    0.4470    0.7410])
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
legend('Predicted Prey')
subplot 413
bar(Years(2:end), X2(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
legend('Real Predator')
subplot 414
bar(Years(2:end), X2_pred(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
legend('Predicted Predator')

% figure
% subplot 211, plot(Ur,'o','Linewidth',2), title('POD Modes - Elements of Ur')
% subplot 212, plot(Years(1:end-1),Vr,'Linewidth',2), title('Time Modes - Elements of Vr')
% xlabel('time [y]'), ylabel('population [#]')

% Reduced order A
Atilde = Ur'*X2*Vr/Sigmar;
X2_pred_red = Atilde*X1;

figure(80)
subplot 411
bar(Years(2:end), X2(1,:), 'FaceColor', [0    0.4470    0.7410])
title('Reduced Order Matrix Fitting')
xlabel('Time [y]')
ylabel('Population [#]')
grid on
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
legend('Real Prey')
subplot 412
bar(Years(2:end), X2_pred_red(1,:), 'FaceColor', [0    0.4470    0.7410])
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
legend('Predicted Prey')
subplot 413
bar(Years(2:end), X2(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
legend('Real Predator')
subplot 414
bar(Years(2:end), X2_pred_red(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
xlabel('Time [y]')
ylabel('Population [#]')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
legend('Predicted Predator')

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
    rH = 17;%size(SigmaH, 1);
end

% figure
% subplot 211, plot(UH(:,1:rH),'Linewidth',2), title('POD Modes - Elements of U'), grid on
% xlabel('time [y]'), ylabel('population [#]'), legend('1', '2', '3')
% subplot 212, plot(Time(1:mc-1),VH(:,1:rH),'Linewidth',2), title('POD Modes - Elements of V'), grid on
% xlabel('time [y]'), ylabel('population [#]'), legend('1', '2', '3')

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
subplot 411
bar(Years, X(1,:), 'FaceColor', [0    0.4470    0.7410])
title('Time Delay Embedding')
legend('Data Prey')
xlabel('Time [y]')
ylabel('Population [#]')
grid on
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
subplot 412
bar(Years, xH_rec(1, :).', 'FaceColor', [0    0.4470    0.7410])
legend('Reconstructed Prey')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
xlabel('Time [y]')
ylabel('Population [#]')
grid on

% Plotting
% figure
subplot 413
bar(Years, X(2,:), 'FaceColor', [ 0.8500    0.3250    0.098])
xlabel('Time [y]')
ylabel('Population [#]')
legend('Data Predator')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
grid on
subplot 414
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
lambda=0.25;
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
% text([1:length(Coeff(Coeff~=0))]+0.05, Coeff(Coeff~=0), string)
% text([1:length(B_min(B_min~=0))]+0.05, B_min(B_min~=0), string)
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

figure(33)
subplot 311
bar(Years, X(1,:), 'FaceColor', [0    0.4470    0.7410]), grid on
title('Lotka-Volterra')
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Prey')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(max(sol_lasso)))])
figure(44)
subplot 311
bar(Years, X(2,:), 'FaceColor', [ 0.8500    0.3250    0.098]), grid on
title('Lotka-Volterra')
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Predator')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(max(sol_lasso)))])

% pinv coefficients ode45 solutions plotting
figure(33)
subplot 312
bar(Years, sol_pinv(1,:), 'FaceColor', [0    0.4470    0.7410]), grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Pseudo Inverse')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(max(sol_lasso)))])
figure(44)
subplot 312
bar(Years, sol_pinv(2,:), 'FaceColor', [ 0.8500    0.3250    0.098]), grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Pseudo Inverse')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(max(sol_lasso)))])

% lasso coefficients ode45 solutions plotting
figure(33)
subplot 313
bar(Years, sol_lasso(1,:,idx1), 'FaceColor', [0    0.4470    0.7410]), grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Lasso')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(max(sol_lasso)))])
figure(44)
subplot 313
bar(Years, sol_lasso(2,:,idx2), 'FaceColor', [ 0.8500    0.3250    0.098]), grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Lasso')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(max(sol_lasso)))])

%% 4.SINDy
clearvars -except A_handle X X1 X2 Years Time dt t;
close all

Xk = X(:,1:end-1)';
Xk1 = X(:,2:end)';

A_handle=@(x1s, x2s) [x1s x2s x1s.^2 x1s.*x2s x2s.^2 x1s.^3 (x2s.^2).*x1s...
                      x2s.^3 cos(x1s) cos(x2s) sin(x1s) sin(x2s) cos(x1s).*cos(x2s)...
                      cos(x1s).*sin(x1s) cos(x2s).*sin(x2s) cos(x1s).*sin(x2s) cos(x2s).*sin(x1s)];

A=A_handle(Xk(:,1), Xk(:,2));

xi_dt = pinv(A)*Xk1;
xi1_dt_pinv = xi_dt(:,1);
xi2_dt_pinv = xi_dt(:,2);
xi1_dt_lasso = lasso(A, Xk1(:,1), 'Lambda', .2);
xi2_dt_lasso = lasso(A, Xk1(:,2), 'Lambda', .2);

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

Xpred_pinv = [A*xi1_dt_pinv, A*xi2_dt_pinv];
Xpred_lasso = [A*xi1_dt_lasso, A*xi2_dt_lasso];

figure
subplot 411
bar(Years(2:end),Xk1(:,1), 'FaceColor', [0    0.4470    0.7410])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Prey')
title('SINDy - Pseudo Inverse')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
subplot 412
bar(Years(2:end),abs(Xpred_pinv(:,1)), 'FaceColor', [0    0.4470    0.7410])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Predicted Prey')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])

% figure
subplot 413
bar(Years(2:end),Xk1(:,2), 'FaceColor', [0.8500    0.3250    0.098])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Predator')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
subplot 414
bar(Years(2:end),abs(Xpred_pinv(:,2)),  'FaceColor', [0.8500    0.3250    0.098])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Predicted Predator')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])

figure
subplot 411
bar(Years(2:end), Xk1(:,1), 'FaceColor', [0    0.4470    0.7410])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Prey')
title('SINDy - Lasso')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
subplot 412
bar(Years(2:end), abs(Xpred_lasso(:,1)), 'FaceColor', [0    0.4470    0.7410])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Predicted Prey')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])

% figure
subplot 413
bar(Years(2:end), Xk1(:,2), 'FaceColor', [0.8500    0.3250    0.098])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Real Predator')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
subplot 414
bar(Years(2:end), abs(Xpred_lasso(:,2)),  'FaceColor', [0.8500    0.3250    0.098])
grid on
xlabel('Time [y]')
ylabel('Population [#]')
legend('Predicted Predator')
axis([Years(1)-dt, max(Years)+dt, 0, 1.1*max(max(X))])
