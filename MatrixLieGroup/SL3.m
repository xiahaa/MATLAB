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

Hreal(:,:,1) = H0;
Hmeas(:,:,1) = H0;
Areal(:,:,1) = A;
for i = 2:N
    Q = rand(3,3).*0.2;
    Q = tosl3(Q);
    Hreal(:,:,i) = Hreal(:,:,i-1)*expm((A).*dt);
    Hmeas(:,:,i) = Hreal(:,:,i-1)*expm((A+Q).*dt);
    Areal(:,:,i) = A;
end

Hest = zeros(3,3,N);
H1 = [3 1 2;1 0.4 1;1 0.4 2];
Hest(:,:,1) = H1./(det(H1)^(1/3));
Aest = zeros(3,3,N);
Aest(:,:,1) = [0 0 0;0 0 0;0 0 0];
for i = 2:N
    [Hest(:,:,i), Aest(:,:,i)] = SL3filter1(Hest(:,:,i-1), Hmeas(:,:,i), Aest(:,:,i-1), dt);
end

figure
for i = 1:1:3
    for j = 1:1:3
        subplot(3,3,(i-1)*3+j);
        b11 = Areal(i,j,:);b11 = b11(:);
        a11 = Aest(i,j,:);a11 = a11(:);
        plot(b11,'r-');hold on;
        plot(a11,'b--');hold on;
    end
end
figure
for i = 1:1:3
    for j = 1:1:3
        subplot(3,3,(i-1)*3+j);
        b11 = Hreal(i,j,:);b11 = b11(:);
        a11 = Hest(i,j,:);a11 = a11(:);
        plot(b11,'r-');hold on;
        plot(a11,'b--');hold on;
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

function [Hest, Aest] = SL3filter1(Hest, Hmeas, Aest, dt)
% Malis E, Hamel T, Mahony R, et al. Dynamic estimation of homography 
% transformations on the special linear group for visual servo control[C]
% 2009 IEEE international conference on robotics and automation. 
% IEEE, 2009: 1498-1503.
    
    Hpred = Hest;  
    Htilde = inv(Hpred)*Hmeas;
    Hproj = tosl3(Htilde'*(eye(3) - Htilde));
    
    kh = 2; ka = 1;
    
    alpha = -kh.*adjointSL3((Htilde),Hproj);
    beta  = -ka.*Hproj;
    
    Adot = beta;
    
    Aest = Aest + Adot*dt;
    Hdot = (adjointSL3((Htilde),Aest) + alpha);
    Hest = Hest*expm((Hdot)*dt);
    Hest = Hest ./ (det(Hest)^(1/3));
end

