clearvars -except BZ_tensor U V Sigma b; close all; clc;

if exist('BZ_tensor', 'var')==0
    try
        load ('C:\Users\davidemaria.fabris\OneDrive - Politecnico di Milano\PhD\PhD Courses\Modelling from Measurements\BZ.mat');
    catch
        load ('C:\Users\david\OneDrive - Politecnico di Milano\PhD\PhD Courses\Modelling from Measurements\BZ.mat');
    end
end

[m,n,o] = size(BZ_tensor);
kk = 1;

%% Data Plotting
for j=1:1:1
    A=BZ_tensor(:,:,j);
    pcolor(A), shading interp, pause(0.001)
    title('Chemical Oscillator Frame')
    xlabel('x [pixel]')
    ylabel('y [pixel]')
    figure
    [x, y] = meshgrid(1:n, 1:m);
    surf(x, y, A, 'FaceAlpha', 1.0, 'EdgeColor', [.25 .25 .25])
    title('Chemical Oscillator 3D Frame')
    xlabel('x [pixel]'), ylabel('y [pixel]'), zlabel('Intensity [-]')
end

%% Data
BZ = permute(BZ_tensor, [3,1,2]);
X = gpuArray(BZ(:,:));
X = X';
dt = 1;
x_0 = X(:,1);

%% Plotting
% figure
% plot(X(1:10000:m*n, :))
% grid on
% title('Pixel History')
% grid on, xlabel('Time'), ylabel('Intensity')

%% DMD
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);

% SVD of X1
if exist('U', 'var')==0 || exist('V', 'var')==0 || exist('Sigma', 'var')==0
    [U, Sigma, V] = svd(X1, 'econ');
end

% Plotting
figure
loglog(diag(Sigma)/sum(diag(Sigma)), 'o', 'Linewidth', 2)
grid on
title('Normalized Singular Values')
xlabel('Mode Number [#]')
ylabel('Normalized Value [-]')

% Reduced order
r=20;
Ur=U(:,1:r);
Sigmar=Sigma(1:r,1:r);
Vr=V(:,1:r);
Atilde = Ur'*X2*Vr/Sigmar;

% Spectral Decomposition
[W,Lambda] = eig(Atilde);

% Spatial Modes of A
Phi=X2*(Vr/Sigmar)*W;

% Frequencies
mu=diag(Lambda);
omega=log(mu)/dt;

% Initial Conditions
if exist('b', 'var')==0 || size(b,1)~=r
    b = Phi\x_0;
end

% for j=1:r
%     image = real(Phi(:,j));
%     A=reshape(image, m, n);
%     figure
%     [x, y] = meshgrid(1:n, 1:m);
%     surf(x, y, A, 'FaceAlpha', 1.0, 'EdgeColor', 'none')
%     title('Chemical Oscillator Modes')
%     xlabel('x [pixel]'), ylabel('y [pixel]'), zlabel('Intensity [-]')
% end

% Reconstruction
x_dmd = gpuArray(zeros(r,0));
for instant = 1:o
    x_dmd(:,instant) =b.*exp(omega*instant);
end
x_rec = real(Phi*x_dmd);

%%
figure
for j=1:5:o
    image = x_rec(:,j);
    A=reshape(image, m, n);
    pcolor(A), shading interp, pause(0.001)
    title('Chemical Oscillator Frame')
    xlabel('x [pixel]')
    ylabel('y [pixel]')
end

figure
for j=1:5:o
    image = x_rec(:,j);
    A=reshape(image, m, n);
    [x, y] = meshgrid(1:n, 1:m);
    surf(x, y, A, 'FaceAlpha', 1.0, 'EdgeColor', 'none'), pause(.001)
    title('Chemical Oscillator 3D')
    xlabel('x [pixel]'), ylabel('y [pixel]'), zlabel('Intensity [-]')
end