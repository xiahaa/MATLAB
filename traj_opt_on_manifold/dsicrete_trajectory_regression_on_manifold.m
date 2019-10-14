function dsicrete_trajectory_regression_on_manifold
    if ismac
        addpath './SOn_regression-master/'
        addpath './SOn_regression-master/STL'
        addpath './utils/'
        addpath './libso3'
    else
        addpath './SOn_regression-master\SOn_regression-master/'
        addpath './jiachao/'
        addpath './utils/'
    end
    clc;close all;clear all;
    
    % Example 2: load from mat file
    data = load('controlpoints.mat');
    n = data.n;
    N = data.N;
    p = data.p;

    % For each control point, pick a weight (positive number). A larger value
    % means the regression curve will pass closer to that control point.
    w = ones(N, 1);

    %% Define parameters of the discrete regression curve

    % The curve has Nd points on SO(n)
    Nd = 97;

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
    lambda = 0.1;

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
    X0 = initguess(problem);
    

    %% my part
    N1 = problem.N;
    N2 = problem.Nd;
    indices =  problem.s;
    tau = problem.delta_tau;
    miu = problem.mu;
    
    Rdata = zeros(3,3*N1);
    Rreg = zeros(3,3*N2);
    
    % fill in data to Rdata
    for i = 1:N1
        Rdata(:,i*3-2:i*3) = X0(:,:,indices(i));
    end
    for i = 1:N2
        Rreg(:,i*3-2:i*3) = X0(:,:,i)*expSO3(0.1*rand(3,1));
    end
    
    % initialize with piecewise geodesic path using park's method
    
    % start optimization
    iter = 1;
    maxiter = 100;
    
    oldcost = -1e6;
    newcost = 1e6;
    
    tol1 = 1e-5;
    
    % seems without trust-region, parallel update will be oscillate.
    % try with sequential update
    % try with quasi-parallel update
    cheeseboard_id = ones(1,N2);
    cheeseboard_id(2:2:N2) = 0;
    cheeseboard_id = logical(cheeseboard_id);% todo, I think I need to use the parallel transport for covariant vector
        
    tr = 1;
    
%     Rreg = traj_opt_by_optimization(Rdata, Rreg, miu, indices, tau);
%     Rreg = seg2seg_seq_sol(Rdata, Rreg, indices, tau, lambda, miu, N2);
    if 1

    tr = 0.1;
    
%     Rreg = traj_opt_by_optimization(Rdata, Rreg, miu, indices, tau);
%     Rreg = traj_smoothing_via_jc(Rreg, indices, 100000, 100);
        
%     [speed0, acc0] = compute_profiles(problem, X0);
    options = optimoptions('quadprog','MaxIterations',100,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'Display','off');
    tic
    while iter < maxiter
%         xi = data_term_error(Rdata,Rreg,indices);
%         v = numerical_diff_v(Rreg);
%         newcost = cost(xi,v,tau,lambda,miu);
%         
%         if abs(newcost - oldcost) < tol1
%             break;
%         end
%         oldcost = newcost;
        
        % sequential update
        newcost = 0;
%         ids = randperm(N2,N2);
        for j = 1:N2
            id = j;%ids(j);
            xi = data_term_error(Rdata,Rreg,indices,id);
            v = numerical_diff_v(Rreg,id);
            dxi = seq_sol(xi, v, indices, tau, lambda, miu, N2, id,Rreg,options);
%             dxis = -LHS(id*3-2:id*3,id*3-2:id*3)\RHS(id*3-2:id*3);
            % 
            if norm(dxi) > tr
                dxi = dxi ./ norm(dxi) .* tr;
            end
            dxis(:,id)=dxi;
            Rreg(:,id*3-2:id*3) = Rreg(:,id*3-2:id*3) * expSO3(dxi);

%             if norm(dxis) > newcost
%                 newcost = norm(dxis);
%             end
        end
%         for j = 1:N2
%             id = j;
%             dxi = dxis(:,id).*cheeseboard_id(id);
%             Rreg(:,id*3-2:id*3) = Rreg(:,id*3-2:id*3) * expSO3(dxi);
%         end
%         cheeseboard_id = ~cheeseboard_id;
        
        % doesnot work
%         Rreg = opt_regression(Rdata, indices, tau, lambda, miu, N2);
        
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        newcost = cost(xi,v,tau,lambda,miu);
        newcosts(iter) = newcost;
        if abs(newcost - oldcost) < tol1
            break;
        end
        oldcost = newcost;
        
        % TODO, do we need to check the norm of the gradient and exit if
        % the norm of gradient is lower than a threshold.
        
        iter = iter + 1;
        disp(iter);
    end
    toc
    
    figure(7);
    plot(newcosts,'r-o','LineWidth',2);
    end
    
    showSO3(Rdata,Rreg);
    
    figure(1);
    plotrotations(X0(:, :, 1:8:Nd));
    view(0, 0);
    
    for i = 1:N2
        X1(:,:,i) = Rreg(:,i*3-2:i*3);
    end

    figure(2);
    plotrotations(X1(:, :, 1:4:Nd));
    view(0, 0);
    
    figure(3);
    plotrotations(X0(:, :, indices));
    view(0, 0);
    figure(4);
    plotrotations(X1(:, :, indices));
    view(0, 0);
    
    
    [speed0, acc0] = compute_profiles(problem, X0);
    [speed1, acc1] = compute_profiles(problem, X1);

    % Passage time of each point on the discrete curves.
    time = problem.delta_tau*( 0 : (problem.Nd-1) );

    figure(5);

    subplot(1, 2, 1);
    plot(1:N2,speed0,1:N2,speed1);
%     plot(time, speed0, time, speed1);
    title('Speed of initial curve and optimized curve');
    xlabel('Time');
    ylabel('Speed');
    legend('Initial curve', 'Optimized curve', 'Location', 'SouthEast');
    pbaspect([1.6, 1, 1]);

    subplot(1, 2, 2);
    plot(1:N2,acc0,1:N2,acc1);
%     plot(time, acc0, time, acc1);
    title('Acceleration of initial curve and optimized curve');
    xlabel('Time');
    ylabel('Acceleration');
    legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);

    ylim([0, 100]);
    
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

    
    posPlot = pos;
    quatPlot = quat;

    
    
    extraTime = 1;
    samplePeriod = 1/100;
    onesVector = ones(extraTime*(1/samplePeriod), 1);
    posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
    quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];

    % Create 6 DOF animation
    SamplePlotFreq = 2;
    Spin = 120;
    filename = 'gait.mp4';
    SixDOFanimation(posPlot, quat2rotm(quatPlot), posd, quat2rotm(quatd),...
                    'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All', ...
                    'Position', [9 39 1280 768], 'View', [(100:(Spin/(length(posPlot)-1)):(100+Spin))', 10*ones(length(posPlot), 1)], ...
                    'AxisLength', 0.1, 'ShowArrowHead', false, ...
                    'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, ...
                    'CreateAVI', false, 'AVIfileNameEnum', filename, 'AVIfps', ((1/samplePeriod) / SamplePlotFreq));

                
    for i = 1:size(X0,3)
        quat0(i,:) = rot2quat(X0(:,:,i))';
    end
    ts = linspace(0,1,size(Rreg,3));
    td = ts(indices);
    t0 = ts;
    
    figure
    subplot(2,2,1);
    plot(ts, quat(:,1),'r-','LineWidth',1.5);
    hold on;grid on;
    plot(td, quatd(:,1),'bs','MarkerSize',5);
    plot(t0, quat0(:,1),'g--','LineWidth',1.5);
    legend({'Regression','Anchor','Geodesics'},'FontName','Arial','FontSize',15);
    title('q0','FontName','Arial','FontSize',15);
    
    subplot(2,2,2);
    plot(ts, quat(:,2),'r-','LineWidth',1.5);
    hold on;grid on;
    plot(td, quatd(:,2),'bs','MarkerSize',5);
    plot(t0, quat0(:,2),'g--','LineWidth',1.5);
    title('q1','FontName','Arial','FontSize',15);

    
    subplot(2,2,3);
    plot(ts, quat(:,3),'r-','LineWidth',1.5);
    hold on;grid on;
    plot(td, quatd(:,3),'bs','MarkerSize',5);
    plot(t0, quat0(:,3),'g--','LineWidth',1.5);
    title('q2','FontName','Arial','FontSize',15);

    
    subplot(2,2,4);
    plot(ts, quat(:,4),'r-','LineWidth',1.5);
    hold on;grid on;
    plot(td, quatd(:,4),'bs','MarkerSize',5);
    plot(t0, quat0(:,4),'g--','LineWidth',1.5);
    title('q3','FontName','Arial','FontSize',15);

                
end

function Rreg = seg2seg_seq_sol(Rdata, Rreg, indices, tau, lambda, miu, N)
    Rreg = reshape(Rreg,3,3,[]);
    options = optimoptions('quadprog','MaxIterations',100,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'Display','off');
    cost2 = inf;
    for k = 1:10
        for i = 2:length(indices)-1
            Ncur = indices(i-1):indices(i+1);
            indicescur = [];
            for ii = 1:length(indices)
                indicescur = [indicescur find(Ncur==indices(ii),1)];
            end
            
            iter = 1;
            maxiter = 200;
            oldcost = inf;
            tol1 = 1e-6;
    %         tr = 1;
            N2 = length(Ncur);
            Rregcur = Rreg(:,:,Ncur);
            Rregcur = reshape(Rregcur,3,[]);
            while iter < maxiter
                newcost = -1e6;
                stidx = 1;
                if i > 2
                    stidx = N2 - (indices(i+1)-indices(i)) + 1;
                end
                for j = stidx:N2
                    id = j;
                    xi = data_term_error(Rdata,Rregcur,indicescur,id);
                    v = numerical_diff_v(Rregcur,id);
                    dxi = seq_sol(xi, v, indicescur, tau, lambda, miu, N2, id, Rregcur, options);
    %                 if norm(dxi) > tr
    %                     dxi = dxi ./ norm(dxi) .* tr;
    %                 end
                    Rregcur(:,id*3-2:id*3) = Rregcur(:,id*3-2:id*3) * expSO3(dxi);
                    if norm(dxi) > newcost
                        newcost = norm(dxi);
                    end
                end
                if abs(newcost - oldcost) < tol1
                    break;
                else
                    oldcost = newcost;
                end
                iter = iter + 1;
            end
            Rregcur = reshape(Rregcur,3,3,[]);
            [speed0, acc0] = compute_profiles_fast(tau,Rregcur);
            figure(1);
            plot(1:N2,speed0,1:N2,acc0);
            Rreg(:,:,Ncur) = Rregcur;
        end
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        cost1 = cost(xi,v,tau,lambda,miu);
        if abs(cost1 - cost2) < 1e-3
            break;
        else
            cost2 = cost1;
        end
    end
    Rreg = reshape(Rreg,3,[]);
end

function Rreg = coarse_to_fine_seq_sol(Rdata, Rreg, indices, tau, lambda, miu, N)
    % start from 30
    Ns = 20;
    N0 = Ns;
    Rreg = reshape(Rreg,3,3,[]);
    options = optimoptions('quadprog','MaxIterations',100,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'Display','off');
    while N0 <= N
        Ncur = round(linspace(1,N,N0));
        for i = 1:length(indices)
            if isempty(find(Ncur == indices(i),1))
                Ncur = [Ncur indices(i)];
            end
        end
        Ncur = sort(Ncur,'ascend');
        indicescur = indices;
        for i = 1:length(indices)
            indicescur(i) = find(Ncur==indices(i),1);
        end
        
        iter = 1;
        maxiter = 200;
        oldcost = inf;
        tol1 = 1e-6;
%         tr = 1;
        N2 = length(Ncur);
        Rregcur = Rreg(:,:,Ncur);
        Rregcur = reshape(Rregcur,3,[]);
        while iter < maxiter
            newcost = -1e6;
            for j = 1:N2
                id = j;
                xi = data_term_error(Rdata,Rregcur,indicescur,id);
                v = numerical_diff_v(Rregcur,id);
                dxi = seq_sol(xi, v, indicescur, tau, lambda, miu, N2, id, Rregcur, options);
%                 if norm(dxi) > tr
%                     dxi = dxi ./ norm(dxi) .* tr;
%                 end
                Rregcur(:,id*3-2:id*3) = Rregcur(:,id*3-2:id*3) * expSO3(dxi);
                if norm(dxi) > newcost
                    newcost = norm(dxi);
                end
            end
            if abs(newcost - oldcost) < tol1
                break;
            else
                oldcost = newcost;
            end
            iter = iter + 1;
        end
        Rregcur = reshape(Rregcur,3,3,[]);
%         [speed0, acc0] = compute_profiles_fast(tau,Rregcur);
%         figure(1);
%         plot(1:N2,speed0,1:N2,acc0);
        Rreg(:,:,Ncur) = Rregcur;
        if N0 >= N
            break;
        else
            N0 = min(N0+Ns,N);
        end
        
    end
    Rreg = reshape(Rreg,3,[]);
end

function [speed, acc] = compute_profiles_fast(dtau,X)
    Nd = size(X,3);
    speed = zeros(1, Nd);
    for k = 1 : Nd-1
        fw = logSO3(X(:, :, k)'*X(:, :, k+1));
        v = fw/dtau;
        speed(k) = sqrt(2)*norm(v);
    end
    % Backward difference for last point.
    speed(end) = speed(end-1);
    % Acceleration is NaN at first and last point. For all the others,
    % using a symmetric difference formula.
    acc = NaN(1, Nd);
    for k = 2 : Nd-1
        fw = logSO3(X(:, :, k)'*X(:, :, k+1));
        bw = logSO3(X(:, :, k)'*X(:, :, k-1));
        a = ( fw + bw ) / ( dtau^2 );
        acc(k) = sqrt(2)*norm(a);
    end

end



function dxi = seq_sol(xi, v, indices, tau, lambda, miu, N, id,Rreg,options)
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
