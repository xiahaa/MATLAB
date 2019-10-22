function [Rreg,newcosts] = pure_optimization_on_so3(Rdata, Rreg, miu, lambda, indices, tau)
    Rdata = reshape(Rdata,3,3,[]);
    Rreg = reshape(Rreg,3,3,[]);
    
    
%     for i = 1:size(Rdata,3)
%         so3data(:,i) = logSO3(Rdata(:,:,i));
%     end
%     
%     for i = 1:size(Rreg,3)
%         so3reg(:,i) = logSO3(Rreg(:,:,i));
%     end

    maxiter = 100;
    oldcost = -1e6;
    newcost = 1e6;
    
    N = round(size(Rreg,3));
    
%     fobj = @(x) (objfunx(x, tau, lambda, miu, indices, Rdata, Rreg));
%     fcons = @(x) (constraint(x, Rreg));
    
    %% nonlinear optimization using continuous representation
    options = optimoptions('quadprog',...
        'Algorithm','trust-region-reflective','Display','iter');
    for i = 1:10
%         xi = data_term_error(Rdata,Rreg,indices);
%         v = numerical_diff_v(Rreg);
%         newcost = cost(xi,v,tau,lambda,miu);
%         newcosts(i) = newcost;
%         if abs(newcost - oldcost) < 1e-6
%             break;
%         end
%         oldcost = newcost;
        
        %% nonlinear optimization
        x0 = zeros(9*(N-1),1);
%         x = fmincon(fobj, x0, [], [], [], [], -1e-1.*ones(length(x0),1), 1e-1.*ones(length(x0),1), fcons);
%         Rreg = reshape(Rreg,3,3,[]);
        [A,f] = objfunx(x0, tau, lambda, miu, indices, Rdata, Rreg);
        [Aeq,beq] = constraintapproax(x0, Rreg);
        x = quadprog(2*A, f, [], [], Aeq, beq, [], [], [], options);

        Rreg = fromso3(Rreg, x);
%         Rreg = reshape(Rreg,3,[]);
        
        if i > 1
            disp(norm(oldx-x));
            if norm(oldx-x)<1e-10, break; end
        end
        oldx = x;
        
        newcosts(i) = x'*A*x+f'*x;
        
    end
%     Rreg = fromso3(Rreg, x);
    Rreg = reshape(Rreg,3,[]);
end

function [Q,fg] = objfunx(x, tau, lambda, mu, indices, Rdata, Rreg)
    Q2 = [  [ 0, 0, 0, 0, 0, 0, 0, 0, 0]; ...
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0]; ...
            [ 0, 0, 4, 0, 0, 0, 0, 0, 0]; ...
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0]; ...
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0]; ...
            [ 0, 0, 0, 0, 0, 4, 0, 0, 0]; ...
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0]; ...
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0]; ...
            [ 0, 0, 0, 0, 0, 0, 0, 0, 4]];
    Q1 = [  [ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
            [ 0, 1,   1, 0, 0,   0, 0, 0,   0]; ...
            [ 0, 1, 4/3, 0, 0,   0, 0, 0,   0]; ...
            [ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
            [ 0, 0,   0, 0, 1,   1, 0, 0,   0]; ...
            [ 0, 0,   0, 0, 1, 4/3, 0, 0,   0]; ...
            [ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
            [ 0, 0,   0, 0, 0,   0, 0, 1,   1]; ...
            [ 0, 0,   0, 0, 0,   0, 0, 1, 4/3]];
    Q1 = Q1 ./ tau * lambda;
    Q2 = Q2 ./ (tau^3) * mu;
    
    N = round(length(x)/9);
    Q1Cell = repmat({Q1}, 1, N);
    BigQ1 = blkdiag(Q1Cell{:});

    Q2Cell = repmat({Q2}, 1, N);
    BigQ2 = blkdiag(Q2Cell{:});
    
    % smoothing terms
%     f = x'*BigQ1*x + x'*BigQ2*x;
    BigQ = BigQ1 + BigQ2;
    
    %% data terms
    fdata = 0;
    A1 = [1 0 0 0 0 0 0 0 0; ...
         0 0 0 1 0 0 0 0 0; ...
         0 0 0 0 0 0 1 0 0];
    A2 = [1 1 1 0 0 0 0 0 0; ...
          0 0 0 1 1 1 0 0 0; ...
          0 0 0 0 0 0 1 1 1];
    Q3 = zeros(size(BigQ));
    fg = zeros(size(Q1,1),1);
    for i = 1:length(indices)
        ii = indices(i);
        if ii <= N
            v = logSO3(Rdata(:,:,i)' * Rreg(:,:,ii));
%             v = logSO3(Rdata(:,:,i)' * Rreg(:,:,ii) * expSO3(A1 * x(ii*9-8:ii*9)));
            Q3(ii*9-8:ii*9,ii*9-8:ii*9) = A1'*A1;
            fg(ii*9-8:ii*9,1) = 2*A1'*v;
        else
            v = logSO3(Rdata(:,:,i)' * Rreg(:,:,ii - 1));
            ii = ii - 1;
%             v = logSO3(Rdata(:,:,i)' * Rreg(:,:,ii) * expSO3(A2 * x(ii*9-8:ii*9)));
            Q3(ii*9-8:ii*9,ii*9-8:ii*9) = A2'*A2;
            fg(ii*9-8:ii*9,1) = 2*A2'*v;
        end
%         fdata = fdata + 2 * v'*v;
    end
    
%     f = f + fdata;
    Q = BigQ + Q3;
%     f = x'*Q*x + fg'*x;
    
end 

function [Aeq,beq] = constraintapproax(x, Rreg)
    A1 = [1 0 0 0 0 0 0 0 0; ...
         0 0 0 1 0 0 0 0 0; ...
         0 0 0 0 0 0 1 0 0];
    A2 = [1 1 1 0 0 0 0 0 0; ...
          0 0 0 1 1 1 0 0 0; ...
          0 0 0 0 0 0 1 1 1];
      
    V1 = [0 1 0 0 0 0 0 0 0; ...
          0 0 0 0 1 0 0 1 0; ...
          0 0 1 0 0 1 0 0 1];
    V2 = [0 1 2 0 0 0 0 0 0; ...
          0 0 0 0 1 2 0 0 0; ...
          0 0 0 0 0 0 0 1 2];
      
    N = round(length(x)/9);
    Aeq = zeros(6*(N-1), length(x));
    beq = zeros(6*(N-1), 1);
    for i = 1:N-1
        ii = i + 1;
%         vi = logSO3();
        vi1 = logSO3(Rreg(:,:,i)'*Rreg(:,:,i+1));
        aa = (eye(3)) * A2;%  + 0.5*hat(vi)
        bb = (eye(3) + 0.5*hat(vi1)) * A1;%  + 0.5*hat(vi);
        Aeq(i*6-5:i*6-3, (i*9-8:i*9)) = aa;
        Aeq(i*6-5:i*6-3, (ii*9-8:ii*9)) = -bb;
        
        beq(i*6-5:i*6-3) = vi1;
%         ceq(i) = norm(*expSO3(A2 * x(i*9-8:i*9)))'*Rreg(i+1)*expSO3(A1 * x(ii*9-8:ii*9))));

        Aeq(i*6-2:i*6, (i*9-8:i*9)) = Rreg(:,:,i+1)*Rreg(:,:,i)'*V2;
        Aeq(i*6-2:i*6, (ii*9-8:ii*9)) = -V1;
        beq(i*6-2:i*6) = 0;
    end
end

function [c,ceq] = constraint(x, Rreg)
    c = [];
    A1 = [1 0 0 0 0 0 0 0 0; ...
         0 0 0 1 0 0 0 0 0; ...
         0 0 0 0 0 0 1 0 0];
    A2 = [1 1 1 0 0 0 0 0 0; ...
          0 0 0 1 1 1 0 0 0; ...
          0 0 0 0 0 0 1 1 1];
    N = round(length(x)/9);
    ceq = zeros(3*(N-1),1);
    for i = 1:N-1
        ii = i + 1;
        ceq(i*3-2:i*3) = logSO3((Rreg(:,:,i)*expSO3(A2 * x(i*9-8:i*9)))'*Rreg(:,:,i+1)*expSO3(A1 * x(ii*9-8:ii*9)));
    end
end

function Q1 = formQ1(t)
    Q1 = [[ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
    [ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
    [ 0, 0, 4*t, 0, 0,   0, 0, 0,   0]; ...
    [ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
    [ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
    [ 0, 0,   0, 0, 0, 4*t, 0, 0,   0]; ...
    [ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
    [ 0, 0,   0, 0, 0,   0, 0, 0,   0]; ...
    [ 0, 0,   0, 0, 0,   0, 0, 0, 4*t]];
end

function Q2 = formQ2(t)
    Q2 = [[ 0,   0,         0, 0,   0,         0, 0,   0,         0]; ...
            [ 0,   t,       t^2, 0,   0,         0, 0,   0,         0]; ...
            [ 0, t^2, (4*t^3)/3, 0,   0,         0, 0,   0,         0]; ...
            [ 0,   0,         0, 0,   0,         0, 0,   0,         0]; ...
            [ 0,   0,         0, 0,   t,       t^2, 0,   0,         0]; ...
            [ 0,   0,         0, 0, t^2, (4*t^3)/3, 0,   0,         0]; ...
            [ 0,   0,         0, 0,   0,         0, 0,   0,         0]; ...
            [ 0,   0,         0, 0,   0,         0, 0,   t,       t^2]; ...
            [ 0,   0,         0, 0,   0,         0, 0, t^2, (4*t^3)/3]];
end

function Rreg = fromso3(Rreg, x)
    A1 = [1 0 0 0 0 0 0 0 0; ...
         0 0 0 1 0 0 0 0 0; ...
         0 0 0 0 0 0 1 0 0];
    A2 = [1 1 1 0 0 0 0 0 0; ...
          0 0 0 1 1 1 0 0 0; ...
          0 0 0 0 0 0 1 1 1];
    for i = 1:size(Rreg,3)-1
        Rreg(:,:,i) = Rreg(:,:,i) * expSO3(A1 * x(i*9-8:i*9));
    end
    Rreg(:,:,end) = Rreg(:,:,end-1) * expSO3(A2 * x(i*9-8:i*9));
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
