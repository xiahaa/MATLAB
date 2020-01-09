clc;close all;clear all;
%% homography SL3 filter

addpath libsl3
addpath filter

%% fist simulation
H0 = [1 0 0;0 1 0;0 0 1];
A = [0.1 0.1 0.1;0.05 -0.18 0.025;0.05 0 0.08];
dt = 0.01;
N = 1500;
Hreal = zeros(3,3,N);
Hmeas = zeros(3,3,N);
Areal = zeros(3,3,N);

Hreal(:,:,1) = H0 ./ (det(H0)^(1/3));
Hmeas(:,:,1) = H0 ./ (det(H0)^(1/3));
Areal(:,:,1) = A;
for i = 2:N
    Q = randn(3,3).*1;
    Q = tosl3(Q);
    H = Hreal(:,:,i-1)*expm((A).*dt);
    Hreal(:,:,i) = H ./ (det(H)^(1/3));
    H = Hreal(:,:,i-1)*expm((A+Q).*dt);
    Hmeas(:,:,i) = H ./ (det(H)^(1/3));
    Areal(:,:,i) = A;
end

Hest = zeros(3,3,N);
H1 = [5 1 2;1 0.4 1;1 0.4 2];
Hest(:,:,1) = H1./(det(H1)^(1/3));
Aest = zeros(3,3,N);
A0 = [0 0 0;0 0 0;0 0 0];
for i = 2:N
    [Hest(:,:,i), Aest(:,:,i)] = SL3filter1(Hest(:,:,i-1), Hmeas(:,:,i), Aest(:,:,i-1), dt);
end

figure
for i = 1:1:3
    for j = 1:1:3
        subplot(3,3,(i-1)*3+j);
        b11 = Areal(i,j,:);b11 = b11(:);
        a11 = Aest(i,j,:);a11 = a11(:);
        plot(b11,'r-','LineWidth',2);hold on;
        plot(a11,'b-.','LineWidth',2);hold on;
        legend({'real','est'});
    end
end
figure
for i = 1:1:3
    for j = 1:1:3
        subplot(3,3,(i-1)*3+j);
        b11 = Hreal(i,j,:);b11 = b11(:);
        a11 = Hest(i,j,:);a11 = a11(:);
        plot(b11,'r-','LineWidth',2);hold on;grid on;
        plot(a11,'b-.','LineWidth',2);hold on;
        legend({'real','est'});
    end
end

%% simulation-2
n=100;
K=eye(3,3);
p=randn(3,n);
p=p./vecnorm(p);

Hest = zeros(3,3,N);
H1 = [5 1 2;1 0.4 1;1 0.4 2];
Hest(:,:,1) = H1./(det(H1)^(1/3));
for i = 2:N
    q = inv(Hreal(:,:,i)) * p;
    q = q./vecnorm(q);
    [Hest(:,:,i)] = SL3filter2(Hest(:,:,i-1), q, p, Areal(:,:,i-1), dt);
end

figure
for i = 1:1:3
    for j = 1:1:3
        subplot(3,3,(i-1)*3+j);
        b11 = Hreal(i,j,:);b11 = b11(:);
        a11 = Hest(i,j,:);a11 = a11(:);
        plot(b11,'r-','LineWidth',2);hold on;grid on;
        plot(a11,'b-.','LineWidth',2);hold on;
        legend({'real','est'});
    end
end

function [Hest, tau] = SL3filter3(Hest, tau, dt, pref, pcur,omega)
% Hamel T, Mahony R, Trumpf J, et al.
% Homography estimation on the special linear group based on direct point correspondence[C]
% //2011 50th IEEE Conference on Decision and Control and European Control Conference. IEEE, 2011: 7902-7908.
% NOTE: do not know how to validate the observor.
    % compute point correction
    w = correctionByPoints(Hest, pref, pcur);
    % update tau
    KI = 1; KP = 3;

    Omega = toso3(omega);
    taudot = liebracket(tau,Omega)+KI*adjointSL3(Hest,w);

    tau = tau + taudot*dt;

    Hdot = Hest*(Omega+tau)+KP*w*Hest;
    Hest = Hest*(eye(3)+(Hdot)*dt);
end
