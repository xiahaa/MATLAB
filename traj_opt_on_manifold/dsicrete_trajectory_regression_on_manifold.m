function dsicrete_trajectory_regression_on_manifold
    addpath './SOn_regression-master\SOn_regression-master/'
    
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
    lambda = 0;

    % Weight of the acceleration regularization term (nonnegative). The larger
    % it is, the more acceleration along the discrete curve is penalized. A
    % large value usually results is a 'straighter' curve (closer to a
    % geodesic.)
    mu = 1e-2;%1e-2;

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
        Rreg(:,i*3-2:i*3) = X0(:,:,i);
    end
    
    % initialize with piecewise geodesic path using park's method
    
    % start optimization
    iter = 1;
    maxiter = 50;
    
    oldcost = -1e6;
    newcost = 1e6;
    
    tol1 = 1e-6;
    
    % seems without trust-region, parallel update will be oscillate.
    % try with sequential update
    % try with quasi-parallel update
    cheeseboard_id = ones(1,N2);
%     cheeseboard_id(2:2:N2) = 0;
%     cheeseboard_id = logical(cheeseboard_id);% todo, I think I need to use the parallel transport for covariant vector
        
    tr = 0.1;
    
    tic
    while iter < maxiter
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        newcost = cost(xi,v,tau,lambda,miu);
%         
        if abs(newcost - oldcost) < tol1
            break;
        end
        oldcost = newcost;
%         
%         % compute gradient, batch optimization
%         % here, I will compute graduate for all so3s and then update as a
%         % batch. I think perhaps sequencial update could be another option
%         % since somehow this is a nonconvex optimization, no global minimum
%         % is guaranteed.
        [LHS, RHS] = batch_sol(xi, v, indices, tau, lambda, miu, N2, Rreg);
        dxis = -LHS\RHS;
%         % update
        for j = 1:N2
            dxi = dxis(j*3-2:j*3).*cheeseboard_id(j);
            if norm(dxi) > tr
                dxi = dxi ./ norm(dxi) .* tr;
            end
            Rreg(:,j*3-2:j*3) = Rreg(:,j*3-2:j*3) * expSO3(dxi);
        end
%         cheeseboard_id = ~cheeseboard_id;
        
        % sequential update
%         newcost = 0;
%         for id = 1:N2
%             xi = data_term_error(Rdata,Rreg,indices,id);
%             v = numerical_diff_v(Rreg,id);
%             [LHS, RHS] = seq_sol(xi, v, indices, tau, lambda, miu, N2, id);
% %             dxis = -LHS(id*3-2:id*3,id*3-2:id*3)\RHS(id*3-2:id*3);
%             dxis = -LHS\RHS;
%             Rreg(:,id*3-2:id*3) = Rreg(:,id*3-2:id*3) * expSO3(dxis);
%             if norm(dxis) > newcost
%                 newcost = norm(dxis);
%             end
%         end
%         xi = data_term_error(Rdata,Rreg,indices);
%         v = numerical_diff_v(Rreg);
%         newcost = cost(xi,v,tau,lambda,miu);
%         if abs(newcost - oldcost) < tol1
%             break;
%         end
%         oldcost = newcost;
        
        % TODO, do we need to check the norm of the gradient and exit if
        % the norm of gradient is lower than a threshold.
        
        iter = iter + 1;
        disp(iter);
    end
    toc
    
    figure(1);
    plotrotations(X0(:, :, 1:8:Nd));
    view(0, 0);
    
    for i = 1:N2
        X1(:,:,i) = Rreg(:,i*3-2:i*3);
    end

    figure(2);
    plotrotations(X1(:, :, 1:8:Nd));
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
    plot(time, speed0, time, speed1);
    title('Speed of initial curve and optimized curve');
    xlabel('Time');
    ylabel('Speed');
    legend('Initial curve', 'Optimized curve', 'Location', 'SouthEast');
    pbaspect([1.6, 1, 1]);

    subplot(1, 2, 2);
    plot(time, acc0, time, acc1);
    title('Acceleration of initial curve and optimized curve');
    xlabel('Time');
    ylabel('Acceleration');
    legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);

    ylim([0, 20]);
    
end

function [LHS, RHS] = seq_sol(xi, v, indices, tau, lambda, miu, N, id)
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
    
    if id == 1
        Jr = rightJinv(-v(:,1));
        lhs = lhs + Jr'*Jr.*c1;
        rhs = rhs + Jr'*(-v(:,1)).*c1;
    elseif id == N
        Jr = rightJinv(v(:,end));
        lhs = lhs + Jr'*Jr.*c1;
        rhs = rhs + Jr'*(v(:,end)).*c1;
    else
        Jr1 = rightJinv(v(:,1));
        Jr2 = rightJinv(-v(:,2));
        A1 = Jr1'*Jr1;
        b1 = Jr1'*v(:,1);
        A2 = Jr2'*Jr2;
        b2 = Jr2'*(-v(:,2));
        lhs = lhs + (A1+A2).*c1;
        rhs = rhs + (b1+b2).*c1;
    end
    
    % third term
    c2 = miu / (tau^3);
    % end points
    if id == 1
        Jr = rightJinv(-v(:,1));
        lhs = lhs + Jr'*Jr.*c2;
        rhs = rhs + Jr'*(-v(:,1)+v(:,2)).*c2;
    elseif id == N
        Jr = rightJinv(v(:,end));
        lhs = lhs + Jr'*Jr.*c2;
        rhs = rhs + Jr'*(-v(:,end-1)+v(:,end)).*c2;
    else
        Jr1 = rightJinv(v(:,1));
        Jr2 = rightJinv(-v(:,2));
        Jr = Jr1+Jr2;
        Jrtb = Jr'*(-v(:,2) + v(:,1));
        lhs = lhs + Jr'*Jr.*c2;
        rhs = rhs + Jrtb.*c2;
    end
    
    if c1 == 0 && c2 == 0
        index = find(indices == id,1);
        if isempty(index)
            lhs = eye(3);
        end
    end
    
    LHS = lhs;
    RHS = rhs;
end

function [LHS, RHS] = batch_sol(xi, v, indices, tau, lambda, miu, N, Rreg)
    lhs = zeros(3,3,N);
    rhs = zeros(3,N);
    
    % fist deal with term 1
    for i = 1:length(indices)
        Jr = rightJinv(xi(:,i));% * Rreg(:,indices(i)*3-2:indices(i)*3);
        lhs(:,:,indices(i)) = lhs(:,:,indices(i)) + Jr'*Jr;
        rhs(:,indices(i)) = rhs(:,indices(i)) + Jr'*xi(:,i);
    end
    
    % second term
    % endpoints 
    c1 = lambda / tau;
    if lambda ~= 0
        Jr = rightJinv(-v(:,1)) * Rreg(:,1:3);
        lhs(:,:,1) = lhs(:,:,1) + Jr'*Jr.*c1;
        rhs(:,1) = rhs(:,1) + Jr'*(-v(:,1)).*c1;

        Jr = rightJinv(v(:,end)) * Rreg(:,end-2:end);
        lhs(:,:,end) = lhs(:,:,end) + Jr'*Jr.*c1;
        rhs(:,end) = rhs(:,end) + Jr'*(v(:,end)).*c1;

        for i = 2:N-1
            Jr1 = rightJinv(v(:,i-1)) * Rreg(:,i*3-2:i*3);
            Jr2 = rightJinv(-v(:,i)) * Rreg(:,i*3-2:i*3);
            A1 = Jr1'*Jr1;
            b1 = Jr1'*v(:,i-1);
            A2 = Jr2'*Jr2;
            b2 = Jr2'*(-v(:,i));
            lhs(:,:,i) = lhs(:,:,i) + (A1+A2).*c1;
            rhs(:,i) = rhs(:,i) + (b1+b2).*c1;
        end
    end
    
    
    % third term
    c2 = miu / (tau^3);
    % end points
    Jr = rightJinv(-v(:,1));% * Rreg(:,1:3);
    lhs(:,:,1) = lhs(:,:,1) + Jr'*Jr.*c2;
    rhs(:,1) = rhs(:,1) + Jr'*(-v(:,1)+v(:,2)).*c2;
    
    Jr = rightJinv(v(:,end));% * Rreg(:,end-2:end);
    lhs(:,:,end) = lhs(:,:,end) + Jr'*Jr.*c2;
    rhs(:,end) = rhs(:,end) + Jr'*(-v(:,end-1)+v(:,end)).*c2;
    
    % 2, two times
    Jr1 = rightJinv(v(:,1));% * Rreg(:,1:3); 
    Jr2 = rightJinv(-v(:,2));% * Rreg(:,1:3);
    A1 = Jr1+Jr2; A1 = A1'*A1;
    b1 = A1'*(-v(:,2)+v(:,1));
    
    A2 = Jr2'*Jr2;
    b2 = Jr2'*(v(:,3)-v(:,2));
    
    lhs(:,:,end) = lhs(:,:,end) + (A1+A2).*c2;
    rhs(:,end) = rhs(:,end) + (b1+b2).*c2;
    % end - 1, two times
    Jr1 = rightJinv(v(:,end-1));% * Rreg(:,1:3); 
    Jr2 = rightJinv(-v(:,end));% * Rreg(:,1:3);
    A1 = Jr1+Jr2; A1 = A1'*A1;
    b1 = A1'*(-v(:,end)+v(:,end-1));
    
    A2 = Jr2'*Jr2;
    b2 = Jr2'*(v(:,end-2)+v(:,end-1));
    
    lhs(:,:,end) = lhs(:,:,end) + (A1+A2).*c2;
    rhs(:,end) = rhs(:,end) + (b1+b2).*c2;
    
    % 3 times
    for i = 3:N-2
        Jr1 = rightJinv(v(:,i-1));% * Rreg(:,i*3-2:i*3);
        Jr2 = rightJinv(-v(:,i));% * Rreg(:,i*3-2:i*3);
        
        Jr3 = rightJinv(-v(:,i));% * Rreg(:,i*3-2:i*3);
        Jr4 = rightJinv(-v(:,i));% * Rreg(:,i*3-2:i*3);
        
        Jr = Jr1+Jr2;
        Jrtb = Jr'*(-v(:,i) + v(:,i-1));
        lhs(:,:,i) = lhs(:,:,i) + Jr'*Jr.*c2;
        rhs(:,i) = rhs(:,i) + Jrtb.*c2;
    end
    
    if c1 == 0 && c2 == 0
        ii = 1:N;
        ii(indices) = [];
        for i = 1:length(ii)
            lhs(:,:,ii(i)) = eye(3);
        end
    end
    
    LHS = spblkdiag(lhs);
    RHS = rhs(:);
end

function xi = data_term_error(Rdata,Rreg,indices,varargin)
    if nargin == 3
        xi = zeros(3,length(indices));
        for i = 1:length(indices)
            ii = indices(i);
            xi(:,i) = logSO3(Rdata(:,(i*3-2):i*3)'*Rreg(:,(ii*3-2):ii*3));
        end
    else
        id = find(indices == varargin{1},1);
        if isempty(id) 
            xi = [];
        else
            ii = indices(id);
            xi = logSO3(Rdata(:,(id*3-2):id*3)'*Rreg(:,(ii*3-2):ii*3));
        end
    end
end

function v = numerical_diff_v(Rreg,varargin)
    N = round(size(Rreg,2)/3);
    if nargin == 1
        v = zeros(3,N-1);
        for i = 1:N-1
            v(:,i) = logSO3(Rreg(:,(i*3-2):i*3)'*Rreg(:,(i*3+1):(i*3+3)));
        end
    else
        i = varargin{1};
        if i == 1
            v(:,1) = logSO3(Rreg(:,(i*3-2):i*3)'*Rreg(:,(i*3+1):(i*3+3)));
            v(:,2) = logSO3(Rreg(:,(i*3+1):i*3+3)'*Rreg(:,(i*3+4):(i*3+6)));
        elseif i == N
            v(:,1) = logSO3(Rreg(:,(i*3-8):i*3-6)'*Rreg(:,(i*3-5):(i*3-3)));
            v(:,2) = logSO3(Rreg(:,(i*3-5):i*3-3)'*Rreg(:,(i*3-2):(i*3)));
        else
            v(:,1) = logSO3(Rreg(:,(i*3-5):i*3-3)'*Rreg(:,(i*3-2):(i*3)));
            v(:,2) = logSO3(Rreg(:,(i*3-2):i*3)'*Rreg(:,(i*3+1):(i*3+3)));
        end
    end
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

function vechat = hat(vec)
    if size(vec,1) == 3 
        vechat = [  0,     -vec(3),  vec(2);
                vec(3),   0    , -vec(1);
               -vec(2),  vec(1),   0    ];  
    elseif size(vec,1) == 6
        vechat = [ hat( vec(4:6,1) ) vec(1:3,1); zeros(1,4) ];    
    end   
end

function R = expSO3(r)
    angle = norm(r);
    tol = 1e-12;
    if angle < tol
        R = eye(3);
        xM = eye(3);
        cmPhi = hat(r);
        N = 10;% finite series
        for n = 1:N
            xM = xM * (cmPhi / n);
            R = R + xM;
        end
        [U,~,V] = svd(R);
        R = V*diag([1,1,det(V*U')])*U';% projection to SO3
    else
        so3 = hat(r);
        R = eye(3) + sin(angle)/angle*so3 + (1-cos(angle))/angle^2*so3^2;
    end
end

function r = logSO3(R)
    r = [0;0;0];
    angle = acos((trace(R)-1)*0.5);
    if angle > 1e-10
        so3 = angle / (2*sin(angle))*(R-R');
        r = [-so3(2,3);so3(1,3);-so3(1,2)];
    else
        so3 = zeros(3,3);
        for i = 1:2
            so3 = so3 + (-1)^(i-1)/i.*(R-eye(3))^(i);
        end
        so3 = 0.5.*(so3-so3');
        r = [-so3(2,3);so3(1,3);-so3(1,2)];
    end
end


function J = leftJ(r)
    tolerance = 1e-12;
    angle = norm(r);
    if angle < tolerance
        % If the angle is small, fall back on the series representation
        N = 10;
        J = eye(3);
        pxn = eye(3);
        px = hat(r);
        for n = 1:N
            pxn = pxn*px/(n + 1);    
            J = J + pxn;
        end
    else
        axis = r/angle;

        cph = (1 - cos(angle))/angle;
        sph = sin(angle)/angle;

        J = sph * eye(3) + (1 - sph) * axis * axis' + cph * hat(axis);
    end       
end

function J = rightJ(r)
    tolerance = 1e-12;
    angle = norm(r);
    if angle < tolerance
        % If the angle is small, fall back on the series representation
        N = 10;
        J = eye(3);
        pxn = eye(3);
        px = -hat(r);
        for n = 1:N
            pxn = pxn*px/(n + 1);    
            J = J + pxn;
        end
    else
        axis = r/angle;

        cph = (1 - cos(angle))/angle;
        sph = sin(angle)/angle;

        J = sph * eye(3) + (1 - sph) * axis * axis' - cph * hat(axis);
    end    
end

function Jinv = leftJinv(r)
    tolerance = 1e-12;
    phi = r;
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation        
        Jinv = eye(3);
        pxn = eye(3);
        px = hat(vec);
        for n = 1:N
            pxn = pxn * px/n;
            Jinv = Jinv + lut_bernoullinumber(n) * pxn;
        end
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        Jinv =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               - ph_2 * hat(axis);
    end   
end

function Jinv = rightJinv(r)
    tolerance = 1e-12;
    phi = r;
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation        
        Jinv = eye(3);
        pxn = eye(3);
        px = -hat(phi);
        for n = 1:10
            pxn = pxn * px/n;
            Jinv = Jinv + lut_bernoullinumber(n) * pxn;
        end
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        Jinv =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               + ph_2 * hat(axis);
    end   
end

function b = lut_bernoullinumber(n)
    lut = [1,-0.500000000000000,0.166666666666667,0,-0.0333333333333333,0, ...
           0.0238095238095238,0,-0.0333333333333334,0,0.0757575757575749,0, ...
           -0.253113553113554,0,1.16666666666667,0,-7.09215686274515,0,54.9711779448621,0,-529.124242424242];
    b = lut(n+1);
end