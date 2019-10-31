clc;clear all;close all;
if ismac
    addpath './SOn_regression-master/'
    addpath './SOn_regression-master/STL'
    addpath './utils/'
    addpath './libso3'
    addpath './jiachao/'
else
    addpath './SOn_regression-master\SOn_regression-master/'
    addpath './jiachao/'
    addpath './utils/'
end
% Example 2: load from mat file
data = load('controlpoints.mat');

%% Define parameters of the discrete regression curve
% The curve has Nd points on SO(n)
Nd = 97;
% Nd = 50;

for i = 1:length(Nd)
    regression_comparison(data, Nd(i));
end

function [cost2] = regression_comparison(data, Nd)
    n = data.n;
    N = data.N;
    p = data.p;

    % For each control point, pick a weight (positive number). A larger value
    % means the regression curve will pass closer to that control point.
    w = ones(N, 1);

    % Each control point attracts one particular point of the regression curve.
    % Specifically, control point k (in 1:N) attracts curve point s(k).
    % The vector s of length N usually satsifies:
    % s(1) = 1, s(end) = Nd and s(k+1) > s(k).
    s = round(linspace(1, Nd, N));
    
    % Time interval between two discretization points of the regression curve.
    % This is only used to fix a scaling. It is useful in particular so that
    % other parameter values such as w, lambda and mu (see below) have the same
    % sense even when the discretization parameter Nd is changed.
    delta_tau = 1/(Nd-1);
    
    % Weight of the velocity regularization term (nonnegative). The larger it
    % is, the more velocity along the discrete curve is penalized. A large
    % value usually results in a shorter curve.
    lambda = 0;

    % Weight of the acceleration regularization term (nonnegative). The larger
    % it is, the more acceleration along the discrete curve is penalized. A
    % large value usually results is a 'straighter' curve (closer to a
    % geodesic.)
    mu = 1e-2;
    
    %% Pack all data defining the regression problem in a problem structure.
    problem.n = n;
    problem.N = N;
    problem.Nd = Nd;
    problem.p = p;
    problem.s = s;
    problem.w = w;
    problem.delta_tau = delta_tau;
    problem.lambda = lambda;
    problem.mu = mu;
    
    %% my part
    N1 = problem.N;
    N2 = problem.Nd;
    tau = problem.delta_tau;
    miu = problem.mu;
    indices = problem.s;
    
    %% Call the optimization procedure to compute the regression curve.
    % Compute an initial guess for the curve. If this step is omitted, digress
    % (below) will compute one itself. X0 is a 3D matrix of size n x n x Nd,
    % such that each slice X0(:, :, k) is a rotation matrix.
    X0 = initguess(problem);
    Rdata = X0(:,:,indices);
    Rdata = reshape(Rdata,3,[]);
    Rreg = reshape(X0,3,[]);

    tolthresh = [0.1,0.3,0.5,0.7,0.9];
    
    for i = 1:length(tolthresh)
        tic;
        [Rreg2,newcost] = non_optimization_on_so3(Rdata, Rreg, miu, lambda, indices, tau, 4, tolthresh(i));
        xi = data_term_error(Rdata,Rreg2,indices);
        v = numerical_diff_v(Rreg2);
        newcosts(i) = fcost(xi,v,tau,lambda,miu,0);
        times(i) = toc;
        Rregs{i} = Rreg2;
    end
    
    ts = ((1:N2)-1).*tau;
    figure(1);
    cmap = lines(length(tolthresh));
    for i = 1:length(Rregs)
        Rreg2 = Rregs{i};
        [speed2, acc2] = compute_profiles(problem, reshape(Rreg2,3,3,[]));
        subplot(2, 1, 1);
        h(i) = plot(ts,speed2,'LineWidth',2,'Color',cmap(i,:));hold on;grid on;
        subplot(2, 1, 2);
        plot(ts,acc2,'LineWidth',2,'Color',cmap(i,:));hold on;grid on;
    end
    
    subplot(2, 1, 1);
    title('Speed Profile');
    xlabel('Time');
    ylabel('Speed');
    legend(h,convertStringsToChars(string(tolthresh)));
    pbaspect([1.6, 1, 1]);
    grid on;
    subplot(2, 1, 2);
    title('Acceleration Profile');
    xlabel('Time');
    ylabel('Acceleration');
%     legend('Initial', 'Trust-Region', 'Proposed', 'Newton', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);
    grid on;
    ylim([0,50]);
    
    showSO3group(Rdata,Rregs,tolthresh,indices);
    
end

function showSO3group(Rdata,Rregs,tolthresh,indices)
    Rdata = reshape(Rdata,3,3,[]);
    for i = 1:size(Rdata,3)
        quatd(i,:) = rot2quat(Rdata(:,:,i))';
    end
    
    ts = linspace(0,1,size(Rregs{1},2)/3);
    td = ts(indices);
    t0 = ts;

    figure
    subplot(2,2,1);
    h1 = plot(td, quatd(:,1),'ks','MarkerSize',8, 'MarkerFaceColor','k');
    hold on;grid on;
    title('q0','FontName','Arial','FontSize',15);

    subplot(2,2,2);
    plot(td, quatd(:,2),'ks','MarkerSize',8, 'MarkerFaceColor','k');
    hold on;grid on;
    title('q1','FontName','Arial','FontSize',15);

    subplot(2,2,3);
    plot(td, quatd(:,3),'ks','MarkerSize',8, 'MarkerFaceColor','k');
    hold on;grid on;
    title('q2','FontName','Arial','FontSize',15);

    subplot(2,2,4);
    plot(td, quatd(:,4),'ks','MarkerSize',8, 'MarkerFaceColor','k');
    hold on;grid on;
    title('q3','FontName','Arial','FontSize',15);
    
    h2 = [];
    for i = 1:length(Rregs)
        Rreg = reshape(Rregs{i},3,3,[]);
        for j = 1:size(Rreg,3)
            quat(j,:) = rot2quat(Rreg(:,:,j))';
        end
    
        subplot(2,2,1);
        h3 = plot(ts, quat(:,1),'-','LineWidth',1.5);
        h2 = [h2 h3];

        subplot(2,2,2);
        plot(ts, quat(:,2),'-','LineWidth',1.5);
        

        subplot(2,2,3);
        plot(ts, quat(:,3),'-','LineWidth',1.5);
        

        subplot(2,2,4);
        plot(ts, quat(:,4),'-','LineWidth',1.5);
        
    end
    leg1 = convertStringsToChars(string(tolthresh));
    legend([h1,h2],{'Input',leg1{:}},'FontName','Arial','FontSize',15);
end