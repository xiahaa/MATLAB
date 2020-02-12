clc;clear all;close all;
addpath './jiachao/'
addpath './utils/'
addpath 'SOn_regression-master/'

t = linspace(0,6*pi,200);
roll = 10 * pi / 180 + sin(t) ;
pitch = 20 * pi / 180 + sin(t);
yaw = 45 * pi / 180 + sin(t);

for i = 1:length(t)
    Rreal(:,:,i) = angle2dcm(yaw(i), pitch(i), roll(i));
    Rnoisy(:,:,i) = angle2dcm(yaw(i)+ 0.1 * randn(1), pitch(i)+ 0.1 * randn(1), roll(i)+ 0.1 * randn(1));
end
tic
Rreg = regression_so3(Rnoisy,Rreal);
toc;

Rreg1 = boumel(Rnoisy);
tic
Rreg2 = traj_smoothing_via_jc(Rnoisy, 1:size(Rreg,3), 10, 10);
toc
Rreg2 = reshape(Rreg2,3,3,[]);

for i = 2:length(t)
    so3 = logSO3(Rreal(:,:,i-1)'*Rreal(:,:,i));
    anglereal(i) = norm(so3) * 180 / pi;
    so3 = logSO3(Rnoisy(:,:,i-1)'*Rnoisy(:,:,i));
    anglenoisy(i) = norm(so3) * 180 / pi;
    so3 = logSO3(Rreg(:,:,i-1)'*Rreg(:,:,i));
    anglereg(i) = norm(so3) * 180 / pi;
    so3 = logSO3(Rreg1(:,:,i-1)'*Rreg1(:,:,i));
    anglereg1(i) = norm(so3) * 180 / pi;
    so3 = logSO3(Rreg2(:,:,i-1)'*Rreg2(:,:,i));
    anglereg2(i) = norm(so3) * 180 / pi;
end
colormap=jet;
ts = linspace(0,1,length(anglereal));
figure
plot(anglereal,'LineWidth',1.5);hold on;grid on;
plot(anglenoisy,'LineWidth',1.5);
plot(anglereg,'LineWidth',1.5);
plot(anglereg1,'LineWidth',1.5);
plot(anglereg2,'LineWidth',1.5);
xlabel('Frame');
ylabel('Relative Rotation Angle: degree');
legend({'Ground Truth','Input','Our','Bou','Jia'},'Location','northwest');



function X1 = regression_so3(Rreg,Rreal)
    if ismac
        addpath './utils/'
        addpath './libso3'
    else
        addpath './utils/'
    end
    n = 3;
    N = size(Rreg,3);
    p = Rreg;
    
    % For each control point, pick a weight (positive number). A larger value
    % means the regression curve will pass closer to that control point.
    w = ones(N, 1);

    %% Define parameters of the discrete regression curve
    Nd = N;

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

    %% Call the optimization procedure to compute the regression curve.
    X0 = Rreg;
    
    %% my part
    N1 = problem.N;
    N2 = problem.Nd;
    indices =  problem.s;
    tau = problem.delta_tau;
    miu = problem.mu;
    
    Rdata = reshape(X0,3,[]);
    Rreg = reshape(X0,3,[]);
    
    % initialize with piecewise geodesic path using park's method
    % start optimization
    iter = 1;
    maxiter = 100;
    
    oldcost = -1e6;
    newcost = 1e6;
    
    tol1 = 1e-4;

    tr = 1;
    
    if 1
%     indices = [];
%     while iter < maxiter        
% %         sequential update
%         newcost = 0;
%         for j = 1:N2
%             id = j;%ids(j);
%             xi = data_term_error(Rdata,Rreg,indices,id);
%             v = numerical_diff_v(Rreg,id);
%             dxi = seq_sol(xi, v, indices, tau, lambda, miu, N2, id);
%             
%             if norm(dxi) > tr
%                 dxi = dxi ./ norm(dxi) .* tr;
%             end
%             dxis(:,id)=dxi;
%             Rreg(:,id*3-2:id*3) = Rreg(:,id*3-2:id*3) * expSO3(dxi);
%         end
%         
%         xi = data_term_error(Rdata,Rreg,indices);
%         v = numerical_diff_v(Rreg);
%         newcost = fcost(xi,v,tau,lambda,miu);
%         newcosts(iter) = newcost;
%         if abs(newcost - oldcost) < tol1
%             break;
%         end
%         oldcost = newcost;
%         
%         iter = iter + 1;
%         disp(iter);
%     end
    
    [Rreg,newcosts] = non_optimization_on_so3(Rdata, Rreg, miu, lambda, indices, tau, 4, 0.25);
    X1 = reshape(Rreg,3,3,[]);
    
    figure(7);
    plot(newcosts,'r-','LineWidth',1.5); grid on;
    xlabel('Iteration');
    ylabel('Objective Function Value');
    end
    
    showSO3(Rdata,Rreg);
    
    for i = 1:N2
        X1(:,:,i) = Rreg(:,i*3-2:i*3);
    end
    
    [speed0, acc0] = compute_profiles(problem, X0);
    [speed1, acc1] = compute_profiles(problem, X1);

    % Passage time of each point on the discrete curves.
    time = problem.delta_tau*( 0 : (problem.Nd-1) );

    figure(5);
    subplot(2, 1, 1);
    plot(1:N2,speed0,1:N2,speed1);
%     plot(time, speed0, time, speed1);
    title('Speed of initial curve and optimized curve');
    xlabel('Frame');
    ylabel('Speed');
    legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);

    subplot(2, 1, 2);
    plot(1:N2,acc0,1:N2,acc1);
%     plot(time, acc0, time, acc1);
    title('Acceleration of initial curve and optimized curve');
    xlabel('Frame');
    ylabel('Acceleration');
    legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);
    
    % Extend final sample to delay end of animation
    Rdata = reshape(Rdata,3,3,[]);
    Rreg = reshape(Rreg,3,3,[]);
    for i = 1:size(Rreg,3)
        quat(i,:) = rot2quat(Rreg(:,:,i))';
        pos(i,:) = (i-1).*[0.01,0.01,0.01];
    end
    
    for i = 1:size(Rdata,3)
        quatd(i,:) = rot2quat(Rdata(:,:,i))';
    end
    posd = pos(indices,:);

    for i = 1:size(Rreal,3)
        quatreal(i,:) = rot2quat(Rreal(:,:,i))';
    end
    posd = pos(indices,:);
    
    posPlot = pos;
    quatPlot = quat;

    
    extraTime = 1;
    samplePeriod = 1/100;
    onesVector = ones(extraTime*(1/samplePeriod), 1);
    posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
    quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];

    % Create 6 DOF animation
    SamplePlotFreq = 2;
                
    for i = 1:size(X0,3)
        quat0(i,:) = rot2quat(X0(:,:,i))';
    end
    ts = linspace(0,1,size(Rreg,3));
    t0 = ts;
    
    figure
    subplot(2,2,1);
    plot(1:N2, quat(:,1),'r-','LineWidth',2);
    hold on;grid on;
    line=plot(1:N2, quat0(:,1),'go','LineWidth',1);
    plot(1:N2, quatreal(:,1),'b-','LineWidth',1);
    alpha(line,0.1);
    legend({'Regression','Ground Truth','Input'},'FontName','Arial','FontSize',15);
    title('q0','FontName','Arial','FontSize',15);
    xlabel('Frame');
%     ylabel('q_0');
    
    subplot(2,2,2);
    plot(1:N2, quat(:,2),'r-','LineWidth',2);
    hold on;grid on;
    plot(1:N2, quatreal(:,2),'b-','LineWidth',1);
    plot(1:N2, quat0(:,2),'go','LineWidth',1);
    title('q1','FontName','Arial','FontSize',15);
    xlabel('Frame');
%     ylabel('q_1');
    
    subplot(2,2,3);
    plot(1:N2, quat(:,3),'r-','LineWidth',2);
    hold on;grid on;
    plot(1:N2, quatreal(:,3),'b-','LineWidth',1);
    plot(1:N2, quat0(:,3),'go','LineWidth',1);
    xlabel('Frame');
    title('q2','FontName','Arial','FontSize',15);
%     ylabel('q_2');
    
    subplot(2,2,4);
    plot(1:N2, quat(:,4),'r-','LineWidth',2);
    hold on;grid on;
    plot(1:N2, quatreal(:,4),'b-','LineWidth',1);
    plot(1:N2, quat0(:,4),'go','LineWidth',1);
    xlabel('Frame');
    title('q3','FontName','Arial','FontSize',15);
%     ylabel('q_3');
                
end




