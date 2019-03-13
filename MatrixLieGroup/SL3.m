clc;close all;clear all;
%% homography SL3 filter

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
    Q = rand(3,3).*0.2;
    Q = tosl3(Q);
    H = Hreal(:,:,i-1)*expm((A).*dt);
    Hreal(:,:,i) = H ./ (det(H)^(1/3));
    H = Hreal(:,:,i-1)*expm((A+Q).*dt);
    Hmeas(:,:,i) = H ./ (det(H)^(1/3));
    Areal(:,:,i) = A;
end

Hest = zeros(3,3,N);
H1 = [3 1 2;1 0.4 1;1 0.4 2];
Hest(:,:,1) = H1./(det(H1)^(1/3));
Aest = zeros(3,3,N);
A0 = [0 0 0;0 0 0;0 0 0];
Hest(:,:,1) = H0;
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
    end
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SL3 libraray                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sl3 = tosl3(H)
% From SL3 to its corresponding sl3
    sl3 = (H - trace(H)/3*eye(3));
end

function n = SL3norm(H1,H2)
% Frobenius norm of SL3
    n = sqrt(trace(H1'*H2));
end

function AdH = adjointSL3(H,X)
    AdH = H*X*inv(H);
end

function res = liebracket(A,B)
    res = A*b-B*A;
end

function res = toso3(w)
    res = [0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
end

function [Hest, Aest] = SL3filter1(Hest, Hmeas, Aest, dt)
% Malis E, Hamel T, Mahony R, et al. Dynamic estimation of homography
% transformations on the special linear group for visual servo control[C]
% 2009 IEEE international conference on robotics and automation.
% IEEE, 2009: 1498-1503.
    Hpred = Hest;
    Htilde = inv(Hpred)*Hmeas;
    Hproj = tosl3(Htilde'*(eye(3) - Htilde));

    kh = 2; ka = 1;

    alpha = -kh*adjointSL3(Htilde,Hproj);
    beta  = -ka*Hproj;

    Adot = beta;

    Aest = Aest + Adot*dt;
    Hdot = (adjointSL3((Htilde),Aest) + alpha);
    Hest = Hest*(eye(3)+(Hdot)*dt);% alternatively, Hest = Hest*expm(Hdot*dt)
    Hest = Hest ./ (det(Hest)^(1/3));

end

function w = correctionByPoints(Hest, pref, pcur)
% be very careful here, Hest transforms current feature to refered
% feature. Another concern is the |e|.
%
    e = Hest*pcur;
    e = e ./ sqrt(e(1,:).^2+e(2,:).^2+e(3,:).^2);
    w = zeros(3,3);
    for i = 1:size(e,2)
        w = w + (eye(3) - e(:,i)*e(:,i)')*pref(:,i)*e(:,i)';
    end
end

% the benefit of this observor is that no direct measuremnt of H is
% necessary
function [Hest, tau] = SL3filter2(Hest, tau, dt, pref, pcur,omega)
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
