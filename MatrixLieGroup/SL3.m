clc;close all;clear all;
%% homography SL3 filter

%% fist simulation
H0 = [1 0 0;0 1 0;0 0 1];
A = [0.1 0.1 0.1;0.05 -0.18 0.025;0.05 0 0.08];
dt = 0.1;
N = 200;
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
Hest(:,:,1) = [3 1 2;1 0 1;1 0.6 2];
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
a11 = Hreal(1,1,:);a11 = a11(:);
b11 = Hest(1,1,:);b11 = b11(:);
plot(a11,'r-');hold on;
plot(b11,'b--');hold on;













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SL3 libraray                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sl3 = tosl3(H)
% From SL3 to its corresponding sl3
    sl3 = (H - trace(H)/3.*eye(3));
end

function n = SL3norm(H1,H2)
% Frobenius norm of SL3
    n = sqrt(trace(H1'*H2));
end

function AdH = adjointSL3(H,X)
    AdH = H*X*inv(H);
end

function [Hest, Aest] = SL3filter1(Hest, Hmeas, Aest, dt)
    Hpred = Hest*expm(Aest.*dt);
    
    Htilde = inv(Hpred)*Hmeas;
    Hproj = tosl3(Htilde' - Htilde'*Htilde);
    
    kh = 2; ka = 1;
    
    alpha = -kh.*adjointSL3(Htilde,Hproj);
    beta  = -ka.*Hproj;
    
    Hdot = (adjointSL3(Htilde,Aest) + alpha);
    Adot = beta;
    
%     Hdot = tosl3(Hdot);
%     Adot = tosl3(Adot);
    
    Hest = Hpred*expm(Hdot.*dt);
    Aest = Aest+Adot.*dt;
end

