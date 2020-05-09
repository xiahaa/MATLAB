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
Nd = 40:20:200;
% Nd = 97;

% for i = 1:length(Nd)
%     [cost1, cost2, cost3] = regression_comparison(data, Nd(i));
%     costs1(i)=cost1;
%     costs2(i)=cost2;
%     costs3(i)=cost3;
% end
% save('res.mat','costs1','costs2','costs3');
load('data/res.mat');
figure
yyaxis left;plot(Nd,[costs1(:).cost],'b-o','LineWidth',2);hold on;grid on;
yyaxis left;plot(Nd,[costs2(:).cost],'r-d','LineWidth',2)
yyaxis left;
ylabel('Cost','FontSize',15,'FontName','Arial')
xlabel('Discrete Number','FontSize',15,'FontName','Arial');
ylabel('Cost','FontSize',15,'FontName','Arial','Color','k')
yyaxis right;plot(Nd,[costs1.time],'b--','LineWidth',2)
yyaxis right;plot(Nd,[costs2.time],'r--','LineWidth',2)
legend({'Cost: Bou','Cost: Ours','Time: Bou','Time: Ours'},'FontSize',15,'FontName','Arial')
ylabel('Runtime: (s)','FontSize',15,'FontName','Arial','Color','k')

function [cost1, cost2, cost3] = regression_comparison(data, Nd)
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
    tic
    [X11, info, optim_problem] = digress(problem, X0);
    cost1.time = toc;
    Rreg1 = reshape(X11,3,[]);
    xi = data_term_error(Rdata,Rreg1,indices);
    v = numerical_diff_v(Rreg1);
    cost1.cost = fcost(xi,v,tau,lambda,miu,0);
    
    % start optimization
    iter = 1;
    maxiter = N2*10;
    oldcost = -1e6;
    newcost = 1e6;
    tol1 = 1e-5;
    tr = 10;
    newcosts = zeros(1,maxiter);
    Rreg2 = Rreg;
    tic
%     while iter <= maxiter
%         % sequential update
%         newcost = 0;
% %         ids = randperm(N2,N2);
%         dxis=[];
%         for j = 1:N2
%             id = j;%ids(j);
%             xi = data_term_error(Rdata,Rreg2,indices,id);
%             v = numerical_diff_v(Rreg2,id);
%             dxi = seq_sol(xi, v, indices, tau, lambda, miu, N2, id);
%             % 
%             if norm(dxi) > tr
%                 dxi = dxi ./ norm(dxi) .* tr;
%             end
%             dxis(:,id)=dxi;
%             Rreg2(:,id*3-2:id*3) = Rreg2(:,id*3-2:id*3) * expSO3(dxi);
%         end
% %         dxis = dxis(:);
% %         Rreg2 = group_update(Rreg2, dxis, N2, ones(1,N2), 1);	
% 
%         xi = data_term_error(Rdata,Rreg2,indices);
%         v = numerical_diff_v(Rreg2);
%         newcost = fcost(xi,v,tau,lambda,miu);
%         newcosts(iter) = newcost;
%         if abs(newcost - oldcost) < tol1
%             break;
%         end
%         oldcost = newcost;
%         iter = iter + 1;
%     end
%     Rreg2 = coarse_to_fine_seq_sol(Rdata, Rreg, indices, tau, lambda, miu, Nd);
%     xi = data_term_error(Rdata,Rreg2,indices);
%     v = numerical_diff_v(Rreg2);
%     newcosts(1) = cost(xi,v,tau,lambda,miu);
    [Rreg2,newcosts] = non_optimization_on_so3(Rdata, Rreg, miu, lambda, indices, tau, 4, 0.7);
    xi = data_term_error(Rdata,Rreg2,indices);
    v = numerical_diff_v(Rreg2);
    newcosts(1) = fcost(xi,v,tau,lambda,miu,0);
    cost2.time = toc;
    cost2.cost = newcosts(end);
    
    for i = 1:length(indices)
        dataloss1(i) = sqrt(2)*norm(logSO3(Rdata(:,i*3-2:i*3)'*Rreg1(:,indices(i)*3-2:indices(i)*3)));
        dataloss2(i) = sqrt(2)*norm(logSO3(Rdata(:,i*3-2:i*3)'*Rreg2(:,indices(i)*3-2:indices(i)*3)));
    end
    
    % jia chao's method
    tic
    Rreg3 = traj_smoothing_via_jc(Rreg, indices, 10, 10);
%     Rreg3 = Rreg;
    cost3.time = toc;
    xi = data_term_error(Rdata,Rreg3,indices);
    v = numerical_diff_v(Rreg3);
    cost3.cost = fcost(xi,v,tau,lambda,miu);
    
    [speed0, acc0] = compute_profiles(problem, X0);
    [speed1, acc1] = compute_profiles(problem, reshape(Rreg1,3,3,[]));
    [speed2, acc2] = compute_profiles(problem, reshape(Rreg2,3,3,[]));
    [speed3, acc3] = compute_profiles(problem, reshape(Rreg3,3,3,[]));

    ts = ((1:N2)-1).*tau;
    figure(1);
    subplot(2, 1, 1);
    plot(1:N2,speed0,1:N2,speed1,1:N2,speed2,1:N2,speed3,'LineWidth',2);
    title('Speed Profile');
    xlabel('Time');
    ylabel('Speed');
    legend('Geodesic Interpolation', 'Bou', 'Ours', 'Jia', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);
    grid on;
    subplot(2, 1, 2);
    plot(1:N2,acc0,1:N2,acc1,1:N2,acc2,1:N2,acc3,'LineWidth',2);
    title('Acceleration Profile');
    xlabel('Time');
    ylabel('Acceleration');
%     legend('Initial', 'Trust-Region', 'Proposed', 'Newton', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);
    grid on;
    ylim([0,50]);
end

function new_group = group_update(old_group, d, N, coding, update_id)
    % update the sequence of SO(3) matrices based on the direction d
    new_group = old_group;
    for j = 1:N
        new_group(:,j*3-2:j*3) = old_group(:,j*3-2:j*3) * expSO3(d(j*3-2:j*3).*(coding(j)==update_id));
    end
end

function grad = grad_descent_armoji(Rdata, Rreg, indices, tau, lambda, miu, N)
    Rreg = reshape(Rreg,3,3,[]);
    Rdata = reshape(Rreg,3,3,[]);
    grad = zeros(3,N);
    
    so3reg = zeros(3,N);
    for i = 1:N
        so3reg(:,i)=logSO3(Rreg(:,:,i));
    end
    
    %% data
    for i = 1:length(indices)
        e = logSO3(Rdata(i)'*Rreg(:,:,indices(i)));
        Jr = rightJ(so3reg(:,indices(i)));
        grad(:,indices(i)) = grad(:,indices(i)) + 2*Jr'*e;
    end
    
    %% velocity
    c1 = 0;%lambda / tau;
    if lambda ~= 0
        i = 1;
        e = logSO3(Rreg(i)'*Rreg(:,:,i+1));
        Jl = leftJ(-so3reg(:,i));
        grad(:,i) = grad(:,i) + c1.*(-2*Jl'*e);
        
        i = N;
        e = logSO3(Rreg(i-1)'*Rreg(:,:,i));
        Jr = rightJ(so3reg(:,i));
        grad(:,i) = grad(:,i) + c1.*(2*Jr'*e);
        
        for i = 2:N-1
            e1 = logSO3(Rreg(i-1)'*Rreg(:,:,i));
            e2 = logSO3(Rreg(i)'*Rreg(:,:,i+1));
            Jr = rightJ(so3reg(:,i));
            Jl = leftJ(-so3reg(:,i));
            grad(:,i) = grad(:,i) + c1.*(-2*Jl'*e2+2*Jr'*e1);
        end
    end
    
    %% acceleration
    c2 = 1;%miu / (tau^3);
    if c2 ~= 0
        i = 1;
        e = logSO3(Rreg(i+1)'*Rreg(:,:,i))+logSO3(Rreg(i+1)'*Rreg(:,:,i+2));
        Jr = rightJ(so3reg(:,i));
        grad(:,i) = grad(:,i) + c2.*(2*Jr'*e);
        
        i = N;
        e = logSO3(Rreg(i-1)'*Rreg(:,:,i))+logSO3(Rreg(i-1)'*Rreg(:,:,i-2));
        Jr = rightJ(so3reg(:,i));
        grad(:,i) = grad(:,i) + c2.*(2*Jr'*e);
        
        i = 2;
        e1 = logSO3(Rreg(i)'*Rreg(:,:,i-1))+logSO3(Rreg(i)'*Rreg(:,:,i+1));
        Jl = leftJ(-so3reg(:,i));
        e2 = logSO3(Rreg(i+1)'*Rreg(:,:,i))+logSO3(Rreg(i+1)'*Rreg(:,:,i+2));
        Jr = rightJ(so3reg(:,i));
        grad(:,i) = grad(:,i) + c2.*(-4*Jl'*e1+2*Jr'*e2);
        
        i = N-1;
        e1 = logSO3(Rreg(i)'*Rreg(:,:,i-1))+logSO3(Rreg(i)'*Rreg(:,:,i+1));
        Jl = leftJ(-so3reg(:,i));
        e2 = logSO3(Rreg(i-1)'*Rreg(:,:,i))+logSO3(Rreg(i-1)'*Rreg(:,:,i-2));
        Jr = rightJ(so3reg(:,i));
        grad(:,i) = grad(:,i) + c2.*(-4*Jl'*e1+2*Jr'*e2);
        
        for i = 3:N-2
            e1 = logSO3(Rreg(i)'*Rreg(:,:,i-1))+logSO3(Rreg(i)'*Rreg(:,:,i+1));
            Jl = leftJ(-so3reg(:,i));
            
            e2 = logSO3(Rreg(i-1)'*Rreg(:,:,i))+logSO3(Rreg(i-1)'*Rreg(:,:,i-2));
            Jr = rightJ(so3reg(:,i));
            
            e3 = logSO3(Rreg(i+1)'*Rreg(:,:,i))+logSO3(Rreg(i+1)'*Rreg(:,:,i+2));
            
            grad(:,i) = grad(:,i) + c2.*(-4*Jl'*e1+2*Jr'*e2+2*Jr'*e3);
        end
    end
    grad = -grad(:);
end

