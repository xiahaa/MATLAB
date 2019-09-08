function varargout = traj_opt_by_optimization(Rdata, Rreg, miu, indices, tau)
    % N here is only used to resampling
    N = size(Rreg,2)/3;
    xreg = group_log(Rreg);
    xid = group_log(Rdata);
    % use SQP for trajectroy regression
    oldcost = Inf;
    oldnorm = Inf;
    maxiter = 50;
    costs = [];
    
    options = optimoptions('quadprog', ...
                'MaxIterations',100, ...
                'OptimalityTolerance',1e-3,'StepTolerance',1e-3,'Display','off');

    
    for i = 1:maxiter
        % qp 
        [xi,cost] = qp_sol(xreg, miu, N, xid, indices, tau, options);
%         cost = ncost(Rdata,Rreg,tau,0,miu,indices);
        costs = [costs cost];
        % update
        Rreg = group_update(xreg,xi);
        % converge
        if abs(cost - oldcost) < 1e-5 || abs(oldnorm - norm(xi)) < 1e-5
            break;
        end
        oldcost = cost;
        oldnorm = norm(xi);
        % expand
        xreg = group_log(Rreg);
    end
    figure(10)
    plot(costs,'r-o');
    varargout{1} = Rreg;
end

function xi = group_log(R)
    N = size(R,2)/3;
    xi = zeros(3,N);
    for i = 1:N
        xi(:,i) = logSO3(R(:,i*3-2:i*3));
    end
end

function R = group_update(xinit,xi)
    N = length(xi)/18;
    R = zeros(3,3,N+1);
    % 1
    x = xi(1:18,1);
    a = x(1:3);b = x(4:6);c=x(7:9);d=x(10:12);e=x(13:15);f=x(16:18);
    R(:,:,1) = expSO3(xinit(:,1)+f);
    for i = 1:N
        x = xi((i-1)*18+1:i*18,1);
        a = x(1:3);b = x(4:6);c=x(7:9);d=x(10:12);e=x(13:15);f=x(16:18);
        R(:,:,i+1) = expSO3(xinit(:,i)+a+b+c+d+e+f);
    end
    R = reshape(R,3,[]);
end

function y = ncost(Rdata,Rreg,tau,lambda,miu,indices)
    xi = data_term_error(Rdata,Rreg,indices);
    v = numerical_diff_v(Rreg);
    N = size(v,2)+1;
    % cost term 1, data cost
    cost1 = sum(vecnorm(xi,2).^2.*2);
    
    if lambda ~= 0
        % cost term 2, first order smooth cost, integrate with trapezoidal
        % rule, consistent with Boumal's paper. TODO change in paper.
        wv = [1 ones(1,N-2)];
        cost2 = sum(vecnorm(v,2).^2.*(2/tau).*wv);
    else
        cost2 = 0;
    end
    
    if miu ~= 0
        % cost term 3, second order smooth cost, integrate with trapezoidal
        % rule
        a = zeros(3,N-2);
        for i = 2:N-1
            a(:,i-1)=v(:,i)-v(:,i-1);
        end
        cost3 = sum(vecnorm(a,2).^2.*(2/tau^3));
    else
        cost3 = 0;
    end
    
    y = cost1 * 0.5 + cost2 * 0.5 * lambda + cost3 * 0.5 * miu;
end

function [x,cost] = qp_sol(xinit, miu, N, xid, indices, tau, options)
    % plan from local chart, then right jacobian will be identity
%     Jr = eye(3);
    
    % firstly, form cost function, T make no pysical meaning in this case
    % could provide different T to weight pysically
    Q = lut_angular_acceleration(tau);
    Qt = [];
    for i = 1:N-1
        Qt = blkdiag(Qt,Q);
    end
    
    Qd = zeros(size(Qt));
    Atmp1 = [1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0; ...
             0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0; ...
             0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1];
    Atmp2 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; ...
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ...
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    f = zeros(size(Qt,1),1);
    costterm = 0;
    for i = 1:size(xid,2)
        % data term loss
        xitmp = logSO3(expSO3(xid(:,i))'*expSO3(xinit(:,indices(i))));
        Jr = rightJinv(xitmp);
        if indices(i) ~= N
            Qd((indices(i)-1)*18+1:indices(i)*18, (indices(i)-1)*18+1:indices(i)*18) = Atmp2'*(Jr'*Jr)*Atmp2;
            f((indices(i)-1)*18+1:indices(i)*18,:) = 2.*(xitmp'*Jr*Atmp2)';
        else
            Qd((N-2)*18+1:(N-1)*18, (N-2)*18+1:(N-1)*18) = Qd((N-2)*18+1:(N-1)*18, (N-2)*18+1:(N-1)*18)+Atmp1'*(Jr'*Jr)*Atmp1;
            f((N-2)*18+1:(N-1)*18,:) = 2.*(xitmp'*Jr*Atmp1)';
        end
        costterm = costterm + xitmp'*xitmp;
    end
    ll = 1;
    Qfull = miu.*Qt + ll*Qd*2;
    f = f.*ll;
    costterm = costterm*ll;
    
    % secondly, continuity 
    Aeq = zeros(9*(N-2),size(Qt,2));
    beq = zeros(size(Aeq,1),1);
    for i = 1:N-2
        [Astart,Aend] = lut_pva(tau,tau);
        Aeq((i-1)*9+1:i*9, (i-1)*18+1:(i-1)*18+36) = [Aend -Astart];
        % r equal plus rd, rdd equal
        beq((i-1)*9+1:i*9,:) = [-xinit(:,i)+xinit(:,i+1);zeros(6,1)];
    end
    
    % thirdly, deviation from anchor points, start and end
    AA1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; ...
           1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0; ...
           0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0; ...
           0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1];
    eps = 0.001;
    Aineq = zeros(12*(N-1),size(Qt,2));
    bineq = zeros(12*(N-1),1);
    for i = 1:N-1
        Aineq((i-1)*12+1:i*12,(i-1)*18+1:i*18) = [AA1;-AA1];
        bineq((i-1)*12+1:i*12,:) = [eps*ones(3,1);xinit(:,i+1)-xinit(:,i)+eps;eps*ones(3,1);-xinit(:,i+1)+xinit(:,i)+eps];
    end
    
    % fourthly, overlap constraint, only valid for image stabilization
    Qfull = sparse(Qfull);
    Aineq = sparse(Aineq);
    Aeq = sparse(Aeq);
    x = quadprog(Qfull,f,Aineq,bineq,Aeq,beq);
    
    cost = 0.5*x'*Qfull*x + f'*x + costterm;
end

function [Astart,Aend] = lut_pva(T1,T2)
    Aend = [ ...
           1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0; ...
           0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0; ...
           0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1; ...
           5 0 0 4 0 0 3 0 0 2 0 0 1 0 0 0 0 1; ...
           0 5 0 0 4 0 0 3 0 0 2 0 0 1 0 0 0 0; ...
           0 0 5 0 0 4 0 0 3 0 0 2 0 0 1 0 0 0; ...
           20 0 0 12 0 0 6 0 0 2 0 0 0 0 0 0 0 0; ...
           0 20 0 0 12 0 0 6 0 0 2 0 0 0 0 0 0 0; ...
           0 0 20 0 0 12 0 0 6 0 0 2 0 0 0 0 0 0];
    Aend(4:6,:) = Aend(4:6,:) / T1;
    Aend(7:9,:) = Aend(7:9,:) / T1^2;
    Astart = [ ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; ...
           0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; ...
           0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0; ...
           0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0; ...
           0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0];
    Astart(4:6,:) = Astart(4:6,:) / T2;
    Astart(7:9,:) = Astart(7:9,:) / T2^2;
end

function Q = lut_angular_velocity(T)
    Q = 1/T.*[[ 25/9,    0,    0,  5/2,    0,    0, 15/7,    0,    0, 5/3,   0,   0, 1, 0, 0, 0, 0, 0]
                [    0, 25/9,    0,    0,  5/2,    0,    0, 15/7,    0,   0, 5/3,   0, 0, 1, 0, 0, 0, 0]
                [    0,    0, 25/9,    0,    0,  5/2,    0,    0, 15/7,   0,   0, 5/3, 0, 0, 1, 0, 0, 0]
                [  5/2,    0,    0, 16/7,    0,    0,    2,    0,    0, 8/5,   0,   0, 1, 0, 0, 0, 0, 0]
                [    0,  5/2,    0,    0, 16/7,    0,    0,    2,    0,   0, 8/5,   0, 0, 1, 0, 0, 0, 0]
                [    0,    0,  5/2,    0,    0, 16/7,    0,    0,    2,   0,   0, 8/5, 0, 0, 1, 0, 0, 0]
                [ 15/7,    0,    0,    2,    0,    0,  9/5,    0,    0, 3/2,   0,   0, 1, 0, 0, 0, 0, 0]
                [    0, 15/7,    0,    0,    2,    0,    0,  9/5,    0,   0, 3/2,   0, 0, 1, 0, 0, 0, 0]
                [    0,    0, 15/7,    0,    0,    2,    0,    0,  9/5,   0,   0, 3/2, 0, 0, 1, 0, 0, 0]
                [  5/3,    0,    0,  8/5,    0,    0,  3/2,    0,    0, 4/3,   0,   0, 1, 0, 0, 0, 0, 0]
                [    0,  5/3,    0,    0,  8/5,    0,    0,  3/2,    0,   0, 4/3,   0, 0, 1, 0, 0, 0, 0]
                [    0,    0,  5/3,    0,    0,  8/5,    0,    0,  3/2,   0,   0, 4/3, 0, 0, 1, 0, 0, 0]
                [    1,    0,    0,    1,    0,    0,    1,    0,    0,   1,   0,   0, 1, 0, 0, 0, 0, 0]
                [    0,    1,    0,    0,    1,    0,    0,    1,    0,   0,   1,   0, 0, 1, 0, 0, 0, 0]
                [    0,    0,    1,    0,    0,    1,    0,    0,    1,   0,   0,   1, 0, 0, 1, 0, 0, 0]
                [    0,    0,    0,    0,    0,    0,    0,    0,    0,   0,   0,   0, 0, 0, 0, 0, 0, 0]
                [    0,    0,    0,    0,    0,    0,    0,    0,    0,   0,   0,   0, 0, 0, 0, 0, 0, 0]
                [    0,    0,    0,    0,    0,    0,    0,    0,    0,   0,   0,   0, 0, 0, 0, 0, 0, 0]];
end

function Q = lut_angular_acceleration(T)
   Q = 1/T^3.*[[ 400/7,     0,     0,    40,     0,     0, 24,  0,  0, 10,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0, 400/7,     0,     0,    40,     0,  0, 24,  0,  0, 10,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0, 400/7,     0,     0,    40,  0,  0, 24,  0,  0, 10, 0, 0, 0, 0, 0, 0]
                [    40,     0,     0, 144/5,     0,     0, 18,  0,  0,  8,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0,    40,     0,     0, 144/5,     0,  0, 18,  0,  0,  8,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0,    40,     0,     0, 144/5,  0,  0, 18,  0,  0,  8, 0, 0, 0, 0, 0, 0]
                [    24,     0,     0,    18,     0,     0, 12,  0,  0,  6,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0,    24,     0,     0,    18,     0,  0, 12,  0,  0,  6,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0,    24,     0,     0,    18,  0,  0, 12,  0,  0,  6, 0, 0, 0, 0, 0, 0]
                [    10,     0,     0,     8,     0,     0,  6,  0,  0,  4,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0,    10,     0,     0,     8,     0,  0,  6,  0,  0,  4,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0,    10,     0,     0,     8,  0,  0,  6,  0,  0,  4, 0, 0, 0, 0, 0, 0]
                [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0]
                [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0]];
end