clc;clear all;close all;
addpath './jiachao/'

t = linspace(0,6*pi,500);
roll = 10 * pi / 180 + sin(t) ;
pitch = 20 * pi / 180 + sin(t);
yaw = 45 * pi / 180 + sin(t);

for i = 1:length(t)
    Rreal(:,:,i) = angle2dcm(yaw(i), pitch(i), roll(i));
    Rnoisy(:,:,i) = angle2dcm(yaw(i)+ 0.1 * randn(1), pitch(i)+ 0.1 * randn(1), roll(i)+ 0.1 * randn(1));
end
tic
Rreg = regression_so3(Rnoisy);
toc;

% Rreg1 = boumel(Rnoisy);
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
%     so3 = logSO3(Rreg1(:,:,i-1)'*Rreg1(:,:,i));
%     anglereg1(i) = norm(so3) * 180 / pi;
    so3 = logSO3(Rreg2(:,:,i-1)'*Rreg2(:,:,i));
    anglereg2(i) = norm(so3) * 180 / pi;
end
colormap=jet;
figure
plot(anglereal,'LineWidth',1.5);hold on;grid on;
plot(anglenoisy,'LineWidth',1.5);
plot(anglereg,'LineWidth',1.5);
plot(anglereg2,'LineWidth',1.5);
xlabel('Frame');
ylabel('Relative Rotation Angle: degree');
legend({'Ground Truth','Input','Proposed','Jia'},'Location','northwest');
% plot(anglereg1);



function X1 = regression_so3(Rreg)
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
    lambda = 1000;

    % Weight of the acceleration regularization term (nonnegative). The larger
    % it is, the more acceleration along the discrete curve is penalized. A
    % large value usually results is a 'straighter' curve (closer to a
    % geodesic.)
    mu = 0;

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
    while iter < maxiter        
        % sequential update
        newcost = 0;
        for j = 1:N2
            id = j;%ids(j);
            xi = data_term_error(Rdata,Rreg,indices,id);
            v = numerical_diff_v(Rreg,id);
            dxi = seq_sol(xi, v, indices, tau, lambda, miu, N2, id);
            % 
            if norm(dxi) > tr
                dxi = dxi ./ norm(dxi) .* tr;
            end
            dxis(:,id)=dxi;
            Rreg(:,id*3-2:id*3) = Rreg(:,id*3-2:id*3) * expSO3(dxi);
        end
        
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        newcost = cost(xi,v,tau,lambda,miu);
        newcosts(iter) = newcost;
        if abs(newcost - oldcost) < tol1
            break;
        end
        oldcost = newcost;
        
        iter = iter + 1;
%         disp(iter);
    end
    X1 = reshape(Rreg,3,3,[]);
%     figure(7);
%     plot(newcosts,'r-','LineWidth',1.5); grid on;
%     xlabel('Iteration');
%     ylabel('Objective Function Value');
%     end
%     
%     showSO3(Rdata,Rreg);
%     
%     for i = 1:N2
%         X1(:,:,i) = Rreg(:,i*3-2:i*3);
%     end
%     
%     [speed0, acc0] = compute_profiles(problem, X0);
%     [speed1, acc1] = compute_profiles(problem, X1);
% 
%     % Passage time of each point on the discrete curves.
%     time = problem.delta_tau*( 0 : (problem.Nd-1) );
% 
%     figure(5);
%     subplot(2, 1, 1);
%     plot(1:N2,speed0,1:N2,speed1);
% %     plot(time, speed0, time, speed1);
%     title('Speed of initial curve and optimized curve');
%     xlabel('Frame');
%     ylabel('Speed');
%     legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
%     pbaspect([1.6, 1, 1]);
% 
%     subplot(2, 1, 2);
%     plot(1:N2,acc0,1:N2,acc1);
% %     plot(time, acc0, time, acc1);
%     title('Acceleration of initial curve and optimized curve');
%     xlabel('Frame');
%     ylabel('Acceleration');
%     legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
%     pbaspect([1.6, 1, 1]);
%     
%     % Extend final sample to delay end of animation
%     Rdata = reshape(Rdata,3,3,[]);
%     Rreg = reshape(Rreg,3,3,[]);
%     
%     for i = 1:size(Rreg,3)
%         quat(i,:) = rot2quat(Rreg(:,:,i))';
%         pos(i,:) = (i-1).*[0.01,0.01,0.01];
%     end
%     
%     for i = 1:size(Rdata,3)
%         quatd(i,:) = rot2quat(Rdata(:,:,i))';
%     end
%     posd = pos(indices,:);
% 
%     
%     posPlot = pos;
%     quatPlot = quat;
% 
%     
%     extraTime = 1;
%     samplePeriod = 1/100;
%     onesVector = ones(extraTime*(1/samplePeriod), 1);
%     posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
%     quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];
% 
%     % Create 6 DOF animation
%     SamplePlotFreq = 2;
%                 
%     for i = 1:size(X0,3)
%         quat0(i,:) = rot2quat(X0(:,:,i))';
%     end
%     ts = linspace(0,1,size(Rreg,3));
%     t0 = ts;
%     
%     figure
%     subplot(2,2,1);
%     plot(1:N2, quat(:,1),'r-','LineWidth',2);
%     hold on;grid on;
%     line=plot(1:N2, quat0(:,1),'go','LineWidth',1);
%     alpha(line,0.1);
%     legend({'Regression','Input'},'FontName','Arial','FontSize',15);
% %     title('q0','FontName','Arial','FontSize',15);
%     xlabel('Frame');
%     ylabel('q0');
%     
%     subplot(2,2,2);
%     plot(1:N2, quat(:,2),'r-','LineWidth',2);
%     hold on;grid on;
%     plot(1:N2, quat0(:,2),'go','LineWidth',1);
% %     title('q1','FontName','Arial','FontSize',15);
%     xlabel('Frame');
%     ylabel('q1');
%     
%     subplot(2,2,3);
%     plot(1:N2, quat(:,3),'r-','LineWidth',2);
%     hold on;grid on;
%     plot(1:N2, quat0(:,3),'go','LineWidth',1);
%     title('q2','FontName','Arial','FontSize',15);
%     xlabel('Frame');
%     ylabel('q2');
%     
%     subplot(2,2,4);
%     plot(1:N2, quat(:,4),'r-','LineWidth',2);
%     hold on;grid on;
%     plot(1:N2, quat0(:,4),'go','LineWidth',1);
%     title('q3','FontName','Arial','FontSize',15);
%     xlabel('Frame');
%     ylabel('q3');
%                 
    end
end


function dxi = seq_sol(xi, v, indices, tau, lambda, miu, N, id)
    lhs = zeros(3,3);
    rhs = zeros(3,1);
    
    if ~isempty(xi)
        Jr = rightJinv(xi);
        lhs = lhs + Jr'*Jr;
        rhs = rhs + Jr'*xi;
    end
        
    % second term
    % endpoints 
    c1 = lambda / tau;
    if lambda ~= 0
        if id == 1
            Jr = rightJinv(v(:,1));
            lhs = lhs + Jr'*Jr.*c1;
            rhs = rhs + Jr'*(v(:,1)).*c1;
        elseif id == N
            Jr = rightJinv(v(:,end));
            lhs = lhs + Jr'*Jr.*c1;
            rhs = rhs + Jr'*(v(:,end)).*c1;
        else
            if id == 2
                id1 = 1;
                id2 = 2;
            else
                id1 = 2;
                id2 = 3;
            end
                
            Jr1 = rightJinv(v(:,id1));
            Jr2 = rightJinv(v(:,id2));
            A1 = Jr1'*Jr1;
            b1 = Jr1'*v(:,id1);
            A2 = Jr2'*Jr2;
            b2 = Jr2'*(v(:,id2));
            lhs = lhs + (A1+A2).*c1;
            rhs = rhs + (b1+b2).*c1;
        end
    end
    
    % add angular velocity constraint
%     eps = 5;
%     if id == 1
%         Jr = rightJinv(v(:,1));
%         Aineq = [Jr./tau;-Jr./tau];
%         bineq = [eps-v(:,1)./tau;eps+v(:,1)./tau];
%     elseif id == N
%         Jr = rightJinv(v(:,end));
%         Aineq = [Jr./tau;-Jr./tau];
%         bineq = [eps-v(:,end)./tau;eps+v(:,end)./tau];
%     else
%         if id == 2
%             id1 = 1;id2 = 2;
%         else
%             id1 = 2;id2 = 3;
%         end
%         Jr1 = rightJinv(v(:,id1));
%         Jr2 = rightJinv(v(:,id2));
%         Aineq = [Jr1./tau;-Jr1./tau;Jr2./tau;-Jr2./tau];
%         bineq = [eps-v(:,id1)./tau;eps+v(:,id1)./tau;eps-v(:,id2)./tau;eps+v(:,id2)./tau];
%     end
    
    % third term
    c2 = miu / (tau^3);
    %% new, use parallel transport and unify all +/-
    ss = 1;
    
    if id == 1
        Jr = rightJinv(v(:,1));% * Rreg(:,1:3)';
        lhs = lhs + Jr'*Jr.*c2;
        rhs = rhs + Jr'*(v(:,1)+v(:,2).*ss).*c2;
    elseif id == N
        Jr = rightJinv(v(:,end));% * Rreg(:,end-2:end)';
        lhs = lhs + Jr'*Jr.*c2;
        rhs = rhs + Jr'*(v(:,end-1).*ss+v(:,end)).*c2;
    elseif id == 2
        % 2, two times
        Jr1 = rightJinv(v(:,1));% * Rreg(:,4:6)'; 
        Jr2 = rightJinv(v(:,2));% * Rreg(:,4:6)';
        A1 = Jr1+Jr2; 
        b1 = A1'*(v(:,2)+v(:,1));A1 = A1'*A1;
    
        A2 = Jr2'*Jr2;
        b2 = Jr2'*(v(:,3).*ss+v(:,2));

        lhs = lhs + (A1+A2).*c2;
        rhs = rhs + (b1+b2).*c2;
    elseif id == N-1
        % end - 1, two times
        Jr1 = rightJinv(v(:,end-1));% * Rreg(:,end-5:end-3)'; 
        Jr2 = rightJinv(v(:,end));% * Rreg(:,end-5:end-3)';
        A1 = Jr1+Jr2; 
        b1 = A1'*(v(:,end)+v(:,end-1));A1 = A1'*A1;

        A2 = Jr1'*Jr1;
        b2 = Jr1'*(v(:,end-2).*ss+v(:,end-1));

        lhs = lhs + (A1+A2).*c2;
        rhs = rhs + (b1+b2).*c2;
    else
        % 3 times
        Jr1 = rightJinv(v(:,2));% * Rreg(:,id*3-2:id*3)';
        Jr2 = rightJinv(v(:,3));% * Rreg(:,id*3-2:id*3)';
        A1 = Jr1+Jr2;
        b1 = A1'*(v(:,3) + v(:,2));A1 = A1'*A1;

        A2 = Jr1;
        b2 = A2'*(v(:,2)+v(:,1).*ss);A2 = A2'*A2;

        A3 = Jr2;
        b3 = A3'*(v(:,4).*ss+v(:,3));A3 = A3'*A3;

        lhs = lhs + (A1+A2+A3).*c2;
        rhs = rhs + (b1+b2+b3).*c2;
    end

    if c1 == 0 && c2 == 0
        index = find(indices == id,1);
        if isempty(index)
            lhs = eye(3);
        end
    end
    
    LHS = lhs;
    RHS = rhs;
    
    dxi = -LHS\RHS;% unconstrained optimization
    
%     %% what if I use constrained optimization.
%     Aineq(dummy1*3+1:end,:) = [];
%     bineq(dummy1*3+1:end,:) = [];
%     amax = 0.01;%sqrt(1000)/2*tau*tau;
%     Aineq2 = [Aineq;-Aineq];
%     bineq2 = [amax-bineq;amax+bineq];
% %     
%     dxi = quadprog(2.*LHS,2*RHS',Aineq,bineq,[],[],[],[],[],options);
end




function y = cost(xi,v,tau,lambda,miu)
    % cost term 1, data cost
    cost1 = sum(vecnorm(xi,2).^2.*2);
    
    % cost term 2, first order smooth cost, integrate with trapezoidal
    % rule, consistent with Boumal's paper. TODO change in paper.
    N = size(v,2)+1;
    wv = [1 ones(1,N-2)];
    cost2 = sum(vecnorm(v,2).^2.*(2/tau).*wv);
    
    % cost term 3, second order smooth cost, integrate with trapezoidal
    % rule
    a = zeros(3,N-2);
    for i = 2:N-1
        a(:,i-1)=v(:,i)-v(:,i-1);
    end
    cost3 = sum(vecnorm(a,2).^2.*(2/tau^3));
    
    y = cost1 * 0.5 + cost2 * 0.5 * lambda + cost3 * 0.5 * miu;
end

function X1 = boumel(Rreg)
    n = 3;
    N = size(Rreg,3);
    p = Rreg;

    % For each control point, pick a weight (positive number). A larger value
    % means the regression curve will pass closer to that control point.
    w = ones(N, 1);

    %% Define parameters of the discrete regression curve
    % The curve has Nd points on SO(n)
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
    lambda = 10;

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

    % Compute an initial guess for the curve. If this step is omitted, digress
    % (below) will compute one itself. X0 is a 3D matrix of size n x n x Nd,
    % such that each slice X0(:, :, k) is a rotation matrix.
    %
    X0 = Rreg;

    % Run the optimization procedure to compute X1, the discrete regression
    % curve. X1 is a 3D matrix of size n x n x Nd with each slice a rotation
    % matrix. The second output, info, is a struct-array containing information
    % returned by the optimization algorithm. The third output, optim_problem,
    % is the Manopt optimization problem structure used to produce X1. It can
    % be used to run another algorithm, e.g., for research purposes.
    %
    [X1, info, optim_problem] = digress(problem, X0);

    %% Plot optimization information
    figure(3);
    semilogy([info.time], [info.gradnorm], 'k.-');
    title('Gradient norm');
    xlabel('Computation time [s]');
    pbaspect([1.6, 1, 1]);


    %% Plot speed and acceleration of X0 and X1

    [speed0, acc0] = compute_profiles(problem, X0);
    [speed1, acc1] = compute_profiles(problem, X1);

    % Passage time of each point on the discrete curves.
    time = problem.delta_tau*( 0 : (problem.Nd-1) );

    figure(4);

    subplot(2, 1, 1);
    plot(time, speed0, time, speed1);
    title('Speed of initial curve and optimized curve');
    xlabel('Time');
    ylabel('Speed');
    legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);

    subplot(2, 1, 2);
    plot(time, acc0, time, acc1);
    title('Acceleration of initial curve and optimized curve');
    xlabel('Time');
    ylabel('Acceleration');
    legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);

%     ylim([0, 20]);
end
