function [Rreg,newcosts] = non_optimization_on_so3(Rdata, Rreg, miu, lambda, indices, tau)
    Rdata = reshape(Rdata,3,3,[]);
    Rreg = reshape(Rreg,3,3,[]);
    
    x0 = zeros(18*(size(Rreg,3)-1),1);

%     fobj = @(x) (objfunx(x, tau, lambda, miu, indices, Rdata, Rreg));
%     fcons = @(x) (constraint(x, Rreg));
    %% nonlinear optimization using continuous representation
    options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','iter');
    for i = 1:2
        %% nonlinear optimization
%         x = fmincon(fobj, x0, [], [], [], []);% -1e-1.*ones(length(x0),1), 1e-1.*ones(length(x0),1), []
%         Rreg = fromso3(Rreg, x);
        
        Rreg = reshape(Rreg,3,3,[]);
        [A,f] = objfunx(x0, tau, lambda, miu, indices, Rdata, Rreg);
        [Aeq,beq] = constraint2(x0);
        [Aieq,bieq] = constraint1(x0, Rdata, indices);
        x = quadprog(2*A, f, Aieq, bieq, Aeq, beq, [], [], [], options);
        Rreg = fromso3(Rreg, x);

        Rreg = reshape(Rreg,3,[]);
        
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        newcost = cost(xi,v,tau,lambda,miu);
        newcosts(i) = newcost;
        
        if i > 1
            disp(norm(oldx-x));
            if norm(oldx-x)<1e-10, break; end
        end
        oldx = x;
    end
    Rreg = reshape(Rreg,3,[]);
end

function [Q,fg] = objfunx(x, tau, lambda, mu, indices, Rdata, Rreg)
    Q1 = [[ 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 1,   1,    1,    1,    1, 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 1, 4/3,  3/2,  8/5,  5/3, 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 1, 3/2,  9/5,    2, 15/7, 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 1, 8/5,    2, 16/7,  5/2, 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 1, 5/3, 15/7,  5/2, 25/9, 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 1,   1,    1,    1,    1, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 1, 4/3,  3/2,  8/5,  5/3, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 1, 3/2,  9/5,    2, 15/7, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 1, 8/5,    2, 16/7,  5/2, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 1, 5/3, 15/7,  5/2, 25/9, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0, 0, 1,   1,    1,    1,    1]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0, 0, 1, 4/3,  3/2,  8/5,  5/3]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0, 0, 1, 3/2,  9/5,    2, 15/7]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0, 0, 1, 8/5,    2, 16/7,  5/2]; ...
            [ 0, 0,   0,    0,    0,    0, 0, 0,   0,    0,    0,    0, 0, 1, 5/3, 15/7,  5/2, 25/9]];
    Q2 = [[ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  4,  6,     8,    10, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  6, 12,    18,    24, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  8, 18, 144/5,    40, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0, 10, 24,    40, 400/7, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  4,  6,     8,    10, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  6, 12,    18,    24, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  8, 18, 144/5,    40, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0, 10, 24,    40, 400/7, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  4,  6,     8,    10]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  6, 12,    18,    24]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0,  8, 18, 144/5,    40]; ...
            [ 0, 0,  0,  0,     0,     0, 0, 0,  0,  0,     0,     0, 0, 0, 10, 24,    40, 400/7]];

    Q1 = Q1 ./ tau * lambda;
    Q2 = Q2 ./ (tau^3) * mu;
    
    N = round(length(x)/18);
    Q1Cell = repmat({Q1}, 1, N);
    BigQ1 = blkdiag(Q1Cell{:});

    Q2Cell = repmat({Q2}, 1, N);
    BigQ2 = blkdiag(Q2Cell{:});
    
    % smoothing terms
%     f = x'*BigQ1*x + x'*BigQ2*x;
    BigQ = BigQ1 + BigQ2;
    
    %% data terms
%     fdata = 0;
%     A1 = [[1 0 0 0 0 0 0 0 0]; ...
%           [0 0 0 1 0 0 0 0 0]; ...
%           [0 0 0 0 0 0 1 0 0]];
%     
    Q3 = zeros(size(BigQ));
    fg = zeros(size(BigQ,1),1);
%     for i = 1:length(indices)
%         ii = indices(i);
% %         if ii <= N
% %         v = logSO3(Rdata(:,:,i)' * expSO3(A1 * x(ii*9-8:ii*9)));
%           v = logSO3(Rdata(:,:,i)');% + A1 * x(ii*9-8:ii*9);
% 
% %             v = logSO3(Rdata(:,:,i)' * Rreg(:,:,ii));
%             Q3(ii*9-8:ii*9,ii*9-8:ii*9) = A1'*A1;
%             fg(ii*9-8:ii*9,1) = 2*A1'*v;
% %         else
% %             v = logSO3(Rdata(:,:,i)' * Rreg(:,:,ii - 1));
% %             ii = ii - 1;
% %             v = logSO3(Rdata(:,:,i)' * expSO3(A2 * x(ii*18-17:ii*18)));
% %             Q3(ii*9-8:ii*9,ii*9-8:ii*9) = A2'*A2;
% %             fg(ii*9-8:ii*9,1) = 2*A2'*v;
% %         end
% %         fdata = fdata + 2 * (v'*v);
%     end
%     
% %     f = f + fdata;
    Q = BigQ + Q3;
    Q = (Q+Q')/2;
%     f = x'*Q*x + fg'*x;
end

function Rreg = fromso3(Rreg, x)
    P0 = [[1 0 0 0 0 0], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 0 0 0 0 0]];
     P1 = [[1 1 1 1 1 1], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 1 1 1 1 1], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 1 1 1 1 1]];
    N = round(length(x)/18);
    for i = 1:N
        Rreg(:,:,i) = expSO3(P0 * x(i*18-17:i*18));
    end
    Rreg(:,:,i+1) = expSO3(P1 * x(i*18-17:i*18));
end

function [Aeq,beq] = constraint2(x)
    P0 = [[1 0 0 0 0 0], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 0 0 0 0 0]];
    V0 = [[0 1 0 0 0 0], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 1 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [0 1 0 0 0 0]];
    A0 = [[0 0 2 0 0 0], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 2 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [0 0 2 0 0 0]];
      
    P1 = [[1 1 1 1 1 1], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 1 1 1 1 1], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 1 1 1 1 1]];
    V1 = [[0 1 2 3 4 5], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 1 2 3 4 5], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [0 1 2 3 4 5]];
    A1 = [[0 0 2 6 12 20], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 2 6 12 20], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [0 0 2 6 12 20]];
     
    N = round(length(x)/18);
    Aeq = zeros(9*(N-1), length(x));
    beq = zeros(9*(N-1), 1);
    for i = 1:N-1
        Aeq(i*9-8:i*9-6, i*18-17:i*18) = P1;
        Aeq(i*9-8:i*9-6, i*18+1:i*18+18) = -P0;
        
        Aeq(i*9-5:i*9-3, i*18-17:i*18) = V1;
        Aeq(i*9-5:i*9-3, i*18+1:i*18+18) = -V0;
        
        Aeq(i*9-2:i*9, i*18-17:i*18) = A1;
        Aeq(i*9-2:i*9, i*18+1:i*18+18) = -A0;
    end
end

function [Aeq,beq] = constraint1(x, Rdata, indices)
    P0 = [[1 0 0 0 0 0], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 0 0 0 0 0]];
      
    P1 = [[1 1 1 1 1 1], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 1 1 1 1 1], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 1 1 1 1 1]];
      
    N = round(length(x)/18);
    Aeq = zeros(6*length(indices), length(x));
    beq = zeros(6*length(indices), 1);
    for i = 1:length(indices)
        ii = indices(i);
        v1 = logSO3(Rdata(:,:,i)');
        
        if ii <= N
            Aeq(i*6-5:i*6-3, (ii*18-17:ii*18)) = P0;
            Aeq(i*6-2:i*6, (ii*18-17:ii*18)) = -P0;
        else
            ii = ii - 1
            Aeq(i*6-5:i*6-3, (ii*18-17:ii*18)) = P1;
            Aeq(i*6-2:i*6, (ii*18-17:ii*18)) = -P1;
        end
            beq(i*6-5:i*6-3) = -v1 + 0.2;
            beq(i*6-2:i*6) = v1 + 0.2;
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
