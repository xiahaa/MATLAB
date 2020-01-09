clc;close all;clear all;
%% homography SL3 filter

addpath libsl3
addpath filter
addpath data

sim3;

function sim1()
    N = 1500;
    dt = 0.01;

    [Hreal, Hmeas, Areal] = datagen_sim1(N,dt);

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
end

function sim2()
    N = 1500;
    dt = 0.01;

    [Hreal, Hmeas, Areal] = datagen_sim1(N,dt);

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
end

function sim3()
    N = 1500;
    dt = 0.01;
    [P,Hreal,Rs,ws,xi,Etareal] = datagen_sim3(N,dt);
    %% simulation-3
    K=eye(3,3);
    p = K*P;
    p=p./vecnorm(p);

    Hest = zeros(3,3,N);
    H1 = [5 1 2;1 0.4 1;1 0.4 2];
    Hest(:,:,1) = H1./(det(H1)^(1/3));
    Eta = zeros(3,3,N);
    for i = 2:N
        Q = Rs(:,:,i)'*(P - xi(:,i));
        q = K*Q + 0.01*randn(3,size(Q,2));
        q = q./vecnorm(q);
        [Hest(:,:,i),Eta(:,:,i)] = SL3filter3(Hest(:,:,i-1), Eta(:,:,i-1), dt, p, q, ws(:,i)+0.01*randn(3,1));
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
    
    figure
    for i = 1:1:3
        for j = 1:1:3
            subplot(3,3,(i-1)*3+j);
            b11 = Etareal(i,j,:);b11 = b11(:);
            a11 = Eta(i,j,:);a11 = a11(:);
            plot(b11,'r-','LineWidth',2);hold on;grid on;
            plot(a11,'b-.','LineWidth',2);hold on;
            legend({'real','est'});
        end
    end
end


