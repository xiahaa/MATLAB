function varargout = traj_opt_by_optimization(Rdata, Rreg, miu, indices, tau)
    % N here is only used to resampling
    N = size(Rreg,2)/3;
    xreg = group_log(Rreg);
    xid = group_log(Rdata);
    % use SQP for trajectroy regression
    oldcost = Inf;
    oldnorm = Inf;
    maxiter = 50;
    costs1 = [];
    costs2 = [];
    for i = 1:maxiter
        cost2 = ncost(Rdata,Rreg,tau,0,miu,indices);
        costs2 = [costs2 cost2];
        % qp 
        [xi,cost1] = qp_sol(xreg, miu, N, xid, indices,tau);
        costs1 = [costs1 cost1];
        % update
        Rreg = group_update(xreg,xi);
        % converge
        if abs(cost2 - oldcost) < 1e-5 || abs(oldnorm - norm(xi)) < 1e-5
            break;
        end
        oldcost = cost2;
        oldnorm = norm(xi);
        % expand
        xreg = group_log(Rreg);
    end
    figure(10)
    plot(costs1,'r-o');hold on;
    figure(11)
    plot(costs2,'r-o');hold on;
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
    R(:,:,1) = expSO3(xinit(:,1))*expSO3(f);
    for i = 1:N
        x = xi((i-1)*18+1:i*18,1);
        a = x(1:3);b = x(4:6);c=x(7:9);d=x(10:12);e=x(13:15);f=x(16:18);
        R(:,:,i+1) = expSO3(xinit(:,i))*expSO3(a+b+c+d+e+f);
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

function [x,cost] = qp_sol(xinit, miu, N, xid, indices,tau)
    % plan from local chart, then right jacobian will be identity
%     Jr = eye(3);
    
    % firstly, form cost function, T make no pysical meaning in this case
    % could provide different T to weight pysically
    Q = lut_angular_acceleration(1);
    Qt = [];
    for i = 1:N-1
        Qt = blkdiag(Qt,Q);
    end
    
    Atmp1 = [1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0; ...
             0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0; ...
             0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1];
    Atmp2 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; ...
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ...
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    
    % firstly, deviation from data points
    Aineq2 = zeros(6*length(indices),size(Qt,2));
    bineq2 = zeros(6*length(indices),1);
    epsd = 0.1;
    for i = 1:size(xid,2)
        if indices(i) ~= N
            Jr = rightJinv(xinit(:,indices(i)));
            Aineq2((i-1)*6+1:i*6,(indices(i)-1)*18+1:indices(i)*18) = [Jr*Atmp2;-Jr*Atmp2];
            bineq2((i-1)*6+1:i*6,:) = [xid(:,i)-xinit(:,indices(i))+epsd;-xid(:,i)+xinit(:,indices(i))+epsd];
        else
            Jr = rightJinv(xinit(:,indices(i)-1));
            Aineq2((i-1)*6+1:i*6,(indices(i)-2)*18+1:indices(i)*18-18) = [Jr*Atmp1;-Jr*Atmp1];
            bineq2((i-1)*6+1:i*6,:) = [xid(:,i)-xinit(:,indices(i))+epsd;-xid(:,i)+xinit(:,indices(i))+epsd];
        end     
    end
    
    Qfull = miu.*Qt;
    f = [];
    
    % secondly, continuity 
    Aeq = zeros(9*(N-2),size(Qt,2));
    beq = zeros(size(Aeq,1),1);
    for i = 1:N-2
        Jr1 = rightJinv(xinit(:,i));
        Jr2 = rightJinv(xinit(:,i+1));
        [Astart,Aend] = lut_pva(1,1,Jr1,Jr2);
        Aeq((i-1)*9+1:i*9, (i-1)*18+1:(i-1)*18+36) = [Aend (-Astart)];
        % r equal plus rd, rdd equal
        beq((i-1)*9+1:i*9,:) = [-xinit(:,i)+xinit(:,i+1);zeros(6,1)];
    end
    
    % thirdly, deviation from anchor points, start and end
    AA1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    AA2 = [ ...
           1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0; ...
           0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0; ...
           0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1];
    eps = 0.01;
    Aineq = zeros(12*(N-1),size(Qt,2));
    bineq = zeros(12*(N-1),1);
    for i = 1:N-1
        Jr = rightJinv(xinit(:,i));
        Aineq((i-1)*12+1:i*12,(i-1)*18+1:i*18) = [Jr*AA1;-Jr*AA1;Jr*AA2;-Jr*AA2];
        bineq((i-1)*12+1:i*12,:) = [eps*ones(3,1);eps*ones(3,1);eps*ones(3,1);eps*ones(3,1)];
    end
    
    % fourthly, overlap constraint, only valid for image stabilization
%     Qfull = sparse(Qfull);
%     Aineq = sparse(Aineq);
%     Aineq2 = sparse(Aineq2);
%     Aeq = sparse(Aeq);
    x = quadprog(Qfull,f,Aineq2,bineq2,Aeq,beq);
    
    cost = 0.5*x'*Qfull*x;
end

function [Astart,Aend] = lut_pva(T1,T2,Jr1,Jr2)
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
    Aend(1:3,:) = Jr1 * Aend(1:3,:);
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
    Astart(1:3,:) = Jr2 * Astart(1:3,:);
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