clc;close all;clear all;
%% homography SL3 filter

addpath libsl3
addpath filter
addpath data

% sim1;
% sim2;
% sim3;
% sim4;
% sim5;
sim6;

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

function sim4()
    N = 1500;
    dt = 0.01;
    [P,Hreal,Rs,ws,xi,~] = datagen_sim3(N,dt);
    
    %% simulation-3
    K=eye(3,3);
    p = K*P;
    Kinv = inv(K);
    pkinv = Kinv * p;
    % from point to line
    l0 = zeros(3,round(size(P,2)*0.5));
    k = 1;
    l0 = cross(pkinv(:,[1:2:end-1]), pkinv(:,[2:2:end]));
    l0 = l0 ./ vecnorm(l0);

    Hest = zeros(3,3,N);
    H1 = [5 1 2;1 0.4 1;1 0.4 2];
    Hest(:,:,1) = H1./(det(H1)^(1/3));
    Eta = zeros(3,3,N);
    for i = 2:N
        Q = Rs(:,:,i)'*(P - xi(:,i));
        q = K*Q + 0.01*randn(3,size(Q,2));
        qkinv = Kinv * q;
        l1 = cross(qkinv(:,[1:2:end-1]), qkinv(:,[2:2:end]));
        l1 = l1 ./ vecnorm(l1);
        [Hest(:,:,i),Eta(:,:,i)] = SL3filter4(Hest(:,:,i-1), Eta(:,:,i-1), dt, l0, l1, ws(:,i)+0.01*randn(3,1));
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

function sim5()
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
        q = inv(Hreal(:,:,i)) * p + randn(3,1)*0.02;
        q = q./vecnorm(q);
        [Hest(:,:,i)] = SL3filter5(Hest(:,:,i-1), q, p, Areal(:,:,i-1), dt);
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

function sim6()
    %% firstly, generate the ladybird fake image
    C1 = diag([1 1/4 -1]);
    C2 = diag([1 1 -1.5^2]);
    C3 = diag([1 1 -1]);
    C5 = [1 0 0.4;0 1 0;0.4 0 0.4^2-0.3^2];
    C8 = [1 0 -0.4;0 1 0;-0.4 0 0.4^2-0.3^2];
    
    C4 = [1 0 0.3;0 1 -0.6;0.3 -0.6 0.3^2+0.6^2-0.2^2];
    C7 = [1 0 -0.3;0 1 -0.6;-0.3 -0.6 0.3^2+0.6^2-0.2^2];
    
    C6 = [1 0 0.3;0 1 0.5;0.3 0.5 0.3^2+0.5^2-0.1^2];
    C9 = [1 0 -0.3;0 1 0.5;-0.3 0.5 0.3^2+0.5^2-0.1^2];
    
    %% just for visualization, check ok
%     syms x y real
%     f1=[x y 1]*C1*[x y 1]';
%     f2=[x y 1]*C2*[x y 1]';
%     f3=[x y 1]*C3*[x y 1]';
%     f5=[x y 1]*C5*[x y 1]';
%     f8=[x y 1]*C8*[x y 1]';
%     f4=[x y 1]*C4*[x y 1]';
%     f7=[x y 1]*C7*[x y 1]';
%     f6=[x y 1]*C6*[x y 1]';
%     f9=[x y 1]*C9*[x y 1]';
%     
%     figure
%     fimplicit(f1);hold on;xlim([-2,2]);
%     fimplicit(f2);
%     fimplicit(f3);
%     fimplicit(f5);
%     fimplicit(f8);
%     fimplicit(f4);
%     fimplicit(f7);
%     fimplicit(f6);
%     fimplicit(f9);
%     axis equal;
%     grid on;
    
    N = 1500;
    dt = 0.005;

    [Hreal, Hmeas, Areal] = datagen_sim1(N,dt);
    
    %% simulation-2
    n=100;
    K=eye(3,3);

    Hest = zeros(3,3,N);
    H1 = [5 1 2;1 0.4 1;1 0.4 2];
    Hest(:,:,1) = H1./(det(H1)^(1/3));
    
    CC = cat(3,C1,C2,C3,C4,C5,C6,C7,C8,C9);
    
    for i = 2:N
        CCmeas = CC;
        for j = 1:size(CC,3)
            CCtrue = Hreal(:,:,i)'*CC(:,:,j)*Hreal(:,:,i);
            CCmeas(:,:,j) = CCtrue;% + randn(3,3)*0.01;
        end
        [Hest(:,:,i)] = SL3filter6(Hmeas(:,:,i-1), CC, CCmeas, dt);
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

