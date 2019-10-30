function varargout = regression(R0,R1)
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
    
    n = 3;
    N = 2;
    p(:,:,1) = R0;
    p(:,:,2) = R1;
    
    w = ones(N, 1);
    % The curve has Nd points on SO(n)
    Nd = 300;
    s = round(linspace(1, Nd, N));

    % Time interval between two discretization points of the regression curve.
    % This is only used to fix a scaling. It is useful in particular so that
    % other parameter values such as w, lambda and mu (see below) have the same
    % sense even when the discretization parameter Nd is changed.
    delta_tau = 1/(Nd-1);

    % Weight of the velocity regularization term (nonnegative). The larger it
    % is, the more velocity along the discrete curve is penalized. A large
    % value usually results in a shorter curve.
    lambda = 100;

    mu = 1e-0;
        
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
        Rreg(:,i*3-2:i*3) = X0(:,:,i);%*expSO3(0.1*rand(3,1));
    end
    
    % initialize with piecewise geodesic path using park's method
    
    % start optimization
    iter = 1;
    maxiter = 100;
    
    oldcost = -1e6;
    newcost = 1e6;
    
    tol1 = 1e-10;
    
    % seems without trust-region, parallel update will be oscillate.
    % try with sequential update
    % try with quasi-parallel update
    cheeseboard_id = ones(1,N2);
    cheeseboard_id(2:2:N2) = 0;
    cheeseboard_id = logical(cheeseboard_id);% todo, I think I need to use the parallel transport for covariant vector
        
    tr = 1;
    
%     Rreg = seg2seg_seq_sol(Rdata, Rreg, indices, tau, lambda, miu, N2);
    if 1

    tr = 0.1;
    
%     Rreg = traj_smoothing_via_jc(Rreg, indices, 100000, 100);
        
%     [speed0, acc0] = compute_profiles(problem, X0);
    options = optimoptions('quadprog','MaxIterations',100,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'Display','off');
    
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
%             Rreg(:,id*3-2:id*3) = Rreg(:,id*3-2:id*3) * expSO3(dxi);

%             if norm(dxis) > newcost
%                 newcost = norm(dxis);
%             end
        end
        for j = 1:N2
            id = j;
%             dxi = dxis(:,id).*cheeseboard_id(id);
            Rreg(:,id*3-2:id*3) = Rreg(:,id*3-2:id*3) * expSO3(dxi);
        end
%         cheeseboard_id = ~cheeseboard_id;
        
        % doesnot work
%         Rreg = opt_regression(Rdata, indices, tau, lambda, miu, N2);
        
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        newcost = fcost(xi,v,tau,lambda,miu);
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
    
    X1 = reshape(Rreg,3,3,[]);
    
    [speed0, acc0] = compute_profiles(problem, X0);
    [speed1, acc1] = compute_profiles(problem, X1);
    costret = sum(acc1(~isnan(acc1)).^2);
    
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
                
    end
    
    varargout{1} = costret;
    varargout{2} = X1;
    varargout{3} = v./tau;
    varargout{4} = acc1;
    varargout{5} = tau;
end






