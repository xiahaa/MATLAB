function [Rreg,newcosts] = non_optimization_on_so3(Rdata, Rreg, miu, lambda, indices, tau, solver)
    Rdata = reshape(Rdata,3,3,[]);
    Rreg = reshape(Rreg,3,3,[]);

    %% nonlinear optimization using continuous representation
%     fobj = @(x) (nonobjfunx(x, tau, lambda, miu, indices, Rdata, Rreg));
%     fcons = @(x) (constraint(x, Rreg));
%     x0 = genNaiveSol(Rreg);

%     fobj = @(x) (objfunx_endpoint_new(x, tau, lambda, miu, indices, Rdata, Rreg));
%     fcons = @(x) (constraint_end1_new(x, Rdata, Rreg, indices));

    %% use qp
    options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','iter');
    for i = 1:100
        if size(Rreg,3)==1
            Rreg = reshape(Rreg,3,3,[]);
        end
        
        %% qp with polynomial parameters
        if solver == 1
            x0 = zeros(18*(size(Rreg,3)-1),1);
            [A,f] = objfunx(x0, tau, lambda, miu, indices, Rdata, Rreg);
            [Aeq,beq] = constraint2(x0, Rreg);
            [Aieq,bieq] = constraint1(x0, Rdata, indices, Rreg);
            x = quadprog(2*A, f, Aieq, bieq, Aeq, beq, [], [], [], options);
            Rreg = fromso3(Rreg, x);
        elseif solver == 2
        %% qp with end-point representation
            x0 = zeros(9*(size(Rreg,3)),1);
            Rreg = reshape(Rreg,3,3,[]);
            [A,f] = objfunx_endpoint(x0, tau, lambda, miu, indices, Rdata, Rreg);
            [Aeq,beq] = constraint_end2(x0);
            [Aieq,bieq] = constraint_end1(x0, Rdata, indices);
            x = quadprog(2*A, f, Aieq, bieq, Aeq, beq, [], [], [], options);
            Rreg = fromso3_end(Rreg, x);
        %% end-point
        elseif solver == 3
            x0 = zeros(3*size(Rreg,3),1);
            [~,~,A] = objfunx_endpoint_new(x0, tau, lambda, miu, indices, Rdata, Rreg);
            [Aeq,beq] = constraint_end2(x0);
            [Aieq,bieq] = constraint_end1_new_lin(x0, Rdata, Rreg, indices);
            x = quadprog(2*A, [], Aieq, bieq, Aeq, beq, [], [], [], options);
            Rreg = fromso3_end_new(Rreg, x);            
        elseif solver == 4
            x0 = zeros(3*size(Rreg,3),1);
    %         for j = 1:size(Rreg,3)
    %             x0(j*3-2:j*3)=logSO3(Rreg(:,:,j));
    %         end
            [~,~,A] = objfunx_endpoint_new(x0, tau, lambda, miu, indices, Rdata, Rreg);
            [Aeq,beq] = constraint_end2(x0);
            fobj = @(x) (objfunx_endpoint_new(x0, tau, lambda, miu, indices, Rdata, Rreg));
            fhess = @(x,lambda) (hessinterior(x,lambda,A));
            fcons = @(x) (constraint_end1_new(x, Rdata, Rreg, indices));
            options = optimoptions(@fmincon,'Algorithm','interior-point',...
                    'Display','final','SpecifyObjectiveGradient',true,'HessianFcn',fhess,'OptimalityTolerance',1e-10,...
                    'StepTolerance',1e-10,'SpecifyConstraintGradient',true);
            [x,fval2,exitflag2,output2]=fmincon(fobj,x0,[],[],[],[],[],[],fcons,options);
            Rreg = fromso3_end_new(Rreg, x);
        end

        Rreg = reshape(Rreg,3,[]);
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        newcost = cost(xi,v,tau,lambda,miu);
        newcosts(i) = newcost;
        
        if i > 1
            disp(['norm diff:', num2str(norm(oldx-x)), ' cost diff:', num2str(abs(oldcost-newcost))]);
            if norm(oldx-x)<1e-6, break; end
            if abs(oldcost-newcost) < 1e-6, break; end
        end
        oldx = x;
        oldcost=newcost;
    end
    Rreg = reshape(Rreg,3,[]);
end


function x0 = genNaiveSol(Rreg)
% generate the naive initial value for nonlinear optimization using fmincon
    A = [1 0 0 0 0 0; ...
         0 1 0 0 0 0; ...
         0 0 1 0 0 0; ...
         1 1 1 1 1 1; ...
         0 1 2 3 4 5; ...
         0 0 2 6 12 20];
    if size(Rreg) ~= 3
        Rreg = reshape(Rreg,3,3,[]);
    end
    
    x0 = zeros(18*(size(Rreg,3)-1),1);
    for i = 1:size(Rreg,3)-1
        p0 = logSO3(Rreg(:,:,i));
        p1 = logSO3(Rreg(:,:,i+1));
        x0(i*18-17:i*18-12) = inv(A)*[p0(1);0;0;p1(1);0;0];
        x0(i*18-11:i*18-6) = inv(A)*[p0(2);0;0;p1(2);0;0];
        x0(i*18-5:i*18) = inv(A)*[p0(3);0;0;p1(3);0;0];
    end
end

function f = nonobjfunx(x, tau, lambda, mu, indices, Rdata, Rreg)
% the objective function of fmincon
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
    f = x'*BigQ1*x + x'*BigQ2*x;
%     BigQ = BigQ1 + BigQ2;
    
    %% data terms
    fdata = 0;
    P0 = [[1 0 0 0 0 0], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 0 0 0 0 0]];
    P1 = [[1 1 1 1 1 1], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 1 1 1 1 1], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 1 1 1 1 1]];
    
%     A1 = [[1 0 0 0 0 0 0 0 0]; ...
%           [0 0 0 1 0 0 0 0 0]; ...
%           [0 0 0 0 0 0 1 0 0]];
%     
%     Q3 = zeros(size(BigQ1));
%     fg = zeros(size(BigQ1,1),1);
%     for i = 1:length(indices)
%         ii = indices(i);
%         if ii <= N
%             v = logSO3(Rdata(:,:,i)' * expSO3(P0 * x(ii*18-17:ii*18)));
% 
% %           v = logSO3(Rdata(:,:,i)');% + A1 * x(ii*9-8:ii*9); 
% % %             v = logSO3(Rdata(:,:,i)' * Rreg(:,:,ii));
% %             Q3(ii*9-8:ii*9,ii*9-8:ii*9) = A1'*A1;
% %             fg(ii*9-8:ii*9,1) = 2*A1'*v;
%         else
%             ii = ii - 1;
%             v = logSO3(Rdata(:,:,i)' * expSO3(P1 * x(ii*18-17:ii*18)));
% %             Q3(ii*9-8:ii*9,ii*9-8:ii*9) = A2'*A2;
% %             fg(ii*9-8:ii*9,1) = 2*A2'*v;
%         end
%         fdata = fdata + 2 * (v'*v);
%     end
    
    f = f + fdata;
%     Q = BigQ + Q3;
%     Q = (Q+Q')/2;
%     f = x'*Q*x + fg'*x;
end

function [Q,fg] = objfunx(x, tau, lambda, mu, indices, Rdata, Rreg)
% objective function for qp with polynomial representation.
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
% recover SO3 
    P0 = [[1 0 0 0 0 0], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 0 0 0 0 0]];
     P1 = [[1 1 1 1 1 1], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 1 1 1 1 1], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 1 1 1 1 1]];
    N = round(length(x)/18);
    for i = 1:N
        Rreg(:,:,i) = Rreg(:,:,i)*expSO3(P0 * x(i*18-17:i*18));
    end
    Rreg(:,:,i+1) = Rreg(:,:,i)*expSO3(P1 * x(i*18-17:i*18));
end

function [Aeq,beq] = constraint2(x,Rreg)
% continuity constraints
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
        v1 = logSO3(Rreg(:,:,i+1)'*Rreg(:,:,i));
        beq(i*9-8:i*9-6) = -v1;
        Jl = leftJinv(v1);
        Jr = rightJinv(v1);
        Aeq(i*9-8:i*9-6, i*18-17:i*18) = Jr*P1;
        Aeq(i*9-8:i*9-6, i*18+1:i*18+18) = -Jl*P0;
        
        Aeq(i*9-5:i*9-3, i*18-17:i*18) = Rreg(:,:,i)*V1;
        Aeq(i*9-5:i*9-3, i*18+1:i*18+18) = -Rreg(:,:,i+1)*V0;
        
        Aeq(i*9-2:i*9, i*18-17:i*18) = Rreg(:,:,i)*A1;
        Aeq(i*9-2:i*9, i*18+1:i*18+18) = -Rreg(:,:,i+1)*A0;
    end
end

function [Aeq,beq] = constraint1(x, Rdata, indices, Rreg)
% endpoint bound constraints
    P0 = [[1 0 0 0 0 0], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 0 0 0 0 0]];
      
    P1 = [[1 1 1 1 1 1], [0 0 0 0 0 0], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [1 1 1 1 1 1], [0 0 0 0 0 0]; ...
          [0 0 0 0 0 0], [0 0 0 0 0 0], [1 1 1 1 1 1]];
      
    N = round(length(x)/18);
    Aeq = zeros(6*length(indices), length(x));
    beq = zeros(6*length(indices), 1);
    tol = 0.5;
    for i = 1:length(indices)
        ii = indices(i);
        
        if ii <= N
            v1 = logSO3(Rdata(:,:,i)'*Rreg(:,:,ii));
            Jr = rightJinv(v1);
            Aeq(i*6-5:i*6-3, (ii*18-17:ii*18)) = Jr*P0;
            Aeq(i*6-2:i*6, (ii*18-17:ii*18)) = -Jr*P0;
        else
            ii = ii - 1;
            v1 = logSO3(Rdata(:,:,i)'*Rreg(:,:,ii));
            Jr = rightJinv(v1);% use eye(3) gets almost the same reuslt
            Aeq(i*6-5:i*6-3, (ii*18-17:ii*18)) = Jr*P1;
            Aeq(i*6-2:i*6, (ii*18-17:ii*18)) = -Jr*P1;
        end
        beq(i*6-5:i*6-3) = -v1 + tol;
        beq(i*6-2:i*6) = v1 + tol;
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

function [Q,fg] = objfunx_endpoint(x, tau, lambda, mu, indices, Rdata, Rreg)
% objective of endpoint representation
    Q1 = [[  10/7,  3/14,   1/84, -10/7,   3/14,  -1/84,     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0]; ...
            [  3/14,  8/35,   1/60, -3/14,  -1/70,  1/210,     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0]; ...
            [  1/84,  1/60,  1/630, -1/84, -1/210, 1/1260,     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0]; ...
            [ -10/7, -3/14,  -1/84,  10/7,  -3/14,   1/84,     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0]; ...
            [  3/14, -1/70, -1/210, -3/14,   8/35,  -1/60,     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0]; ...
            [ -1/84, 1/210, 1/1260,  1/84,  -1/60,  1/630,     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0]; ...
            [     0,     0,      0,     0,      0,      0,  10/7,  3/14,   1/84, -10/7,   3/14,  -1/84,     0,     0,      0,     0,      0,      0]; ...
            [     0,     0,      0,     0,      0,      0,  3/14,  8/35,   1/60, -3/14,  -1/70,  1/210,     0,     0,      0,     0,      0,      0]; ...
            [     0,     0,      0,     0,      0,      0,  1/84,  1/60,  1/630, -1/84, -1/210, 1/1260,     0,     0,      0,     0,      0,      0]; ...
            [     0,     0,      0,     0,      0,      0, -10/7, -3/14,  -1/84,  10/7,  -3/14,   1/84,     0,     0,      0,     0,      0,      0]; ...
            [     0,     0,      0,     0,      0,      0,  3/14, -1/70, -1/210, -3/14,   8/35,  -1/60,     0,     0,      0,     0,      0,      0]; ...
            [     0,     0,      0,     0,      0,      0, -1/84, 1/210, 1/1260,  1/84,  -1/60,  1/630,     0,     0,      0,     0,      0,      0]; ...
            [     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0,  10/7,  3/14,   1/84, -10/7,   3/14,  -1/84]; ...
            [     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0,  3/14,  8/35,   1/60, -3/14,  -1/70,  1/210]; ...
            [     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0,  1/84,  1/60,  1/630, -1/84, -1/210, 1/1260]; ...
            [     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0, -10/7, -3/14,  -1/84,  10/7,  -3/14,   1/84]; ...
            [     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0,  3/14, -1/70, -1/210, -3/14,   8/35,  -1/60]; ...
            [     0,     0,      0,     0,      0,      0,     0,     0,      0,     0,      0,      0, -1/84, 1/210, 1/1260,  1/84,  -1/60,  1/630]];
    Q2 = [[  120/7,   60/7,   3/7, -120/7,   60/7,   -3/7,      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0]; ...
            [   60/7, 192/35, 11/35,  -60/7, 108/35,  -4/35,      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0]; ...
            [    3/7,  11/35,  3/35,   -3/7,   4/35,   1/70,      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0]; ...
            [ -120/7,  -60/7,  -3/7,  120/7,  -60/7,    3/7,      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0]; ...
            [   60/7, 108/35,  4/35,  -60/7, 192/35, -11/35,      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0]; ...
            [   -3/7,  -4/35,  1/70,    3/7, -11/35,   3/35,      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0]; ...
            [      0,      0,     0,      0,      0,      0,  120/7,   60/7,   3/7, -120/7,   60/7,   -3/7,      0,      0,     0,      0,      0,      0]; ...
            [      0,      0,     0,      0,      0,      0,   60/7, 192/35, 11/35,  -60/7, 108/35,  -4/35,      0,      0,     0,      0,      0,      0]; ...
            [      0,      0,     0,      0,      0,      0,    3/7,  11/35,  3/35,   -3/7,   4/35,   1/70,      0,      0,     0,      0,      0,      0]; ...
            [      0,      0,     0,      0,      0,      0, -120/7,  -60/7,  -3/7,  120/7,  -60/7,    3/7,      0,      0,     0,      0,      0,      0]; ...
            [      0,      0,     0,      0,      0,      0,   60/7, 108/35,  4/35,  -60/7, 192/35, -11/35,      0,      0,     0,      0,      0,      0]; ...
            [      0,      0,     0,      0,      0,      0,   -3/7,  -4/35,  1/70,    3/7, -11/35,   3/35,      0,      0,     0,      0,      0,      0]; ...
            [      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0,  120/7,   60/7,   3/7, -120/7,   60/7,   -3/7]; ...
            [      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0,   60/7, 192/35, 11/35,  -60/7, 108/35,  -4/35]; ...
            [      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0,    3/7,  11/35,  3/35,   -3/7,   4/35,   1/70]; ...
            [      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0, -120/7,  -60/7,  -3/7,  120/7,  -60/7,    3/7]; ...
            [      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0,   60/7, 108/35,  4/35,  -60/7, 192/35, -11/35]; ...
            [      0,      0,     0,      0,      0,      0,      0,      0,     0,      0,      0,      0,   -3/7,  -4/35,  1/70,    3/7, -11/35,   3/35]];

    Q1 = Q1 ./ tau * lambda;
    Q2 = Q2 ./ (tau^3) * mu;
    
    N = round(length(x)/9);
    BigQ = zeros(9*N,9*N);
    
    C = [eye(3) zeros(3,15); ...
         zeros(3,9) eye(3) zeros(3,6); ...
         zeros(3) eye(3) zeros(3,12); ...
         zeros(3,12) eye(3) zeros(3); ...
         zeros(3,6) eye(3) zeros(3,9); ...
         zeros(3,15), eye(3)];
    
    CtQ1C = C'*Q1*C;
    CtQ2C = C'*Q2*C;
     
    for i = 1:N-1
        BigQ(i*9-8:i*9+9,i*9-8:i*9+9) = BigQ(i*9-8:i*9+9,i*9-8:i*9+9) + CtQ1C + CtQ2C;
    end
    
    %% data terms
    Q3 = zeros(size(BigQ));
    fg = zeros(size(BigQ,1),1);

    Q = BigQ + Q3;
    Q = (Q+Q')/2;
end


function [Aeq,beq] = constraint_end2(x)
% if use endpoint, there is no equality constraints
    Aeq = [];
    beq = [];
end

function [Aeq,beq] = constraint_end1(x, Rdata, indices)
% endpoint bound constraints
    P0 = [[1 0 0], [0 0 0], [0 0 0]; ...
          [0 0 0], [1 0 0], [0 0 0]; ...
          [0 0 0], [0 0 0], [1 0 0]];
      
    N = round(length(x)/9);
    Aeq = zeros(6*length(indices), length(x));
    beq = zeros(6*length(indices), 1);
    tol = 0.3;
    for i = 1:length(indices)
        ii = indices(i);
        v1 = logSO3(Rdata(:,:,i)');
        
        Aeq(i*6-5:i*6-3, (ii*9-8:ii*9)) = P0;
        Aeq(i*6-2:i*6, (ii*9-8:ii*9)) = -P0;
        
        beq(i*6-5:i*6-3) = -v1 + tol;
        beq(i*6-2:i*6) = v1 + tol;
    end
end

function Rreg = fromso3_end(Rreg, x)
% recover SO3 from endpoint representation.
    P0 = [[1 0 0], [0 0 0], [0 0 0]; ...
          [0 0 0], [1 0 0], [0 0 0]; ...
          [0 0 0], [0 0 0], [1 0 0]];
    N = round(length(x)/9);
    for i = 1:N
        Rreg(:,:,i) = expSO3(P0 * x(i*9-8:i*9));
    end
end

function [f,gradf,Qout] = objfunx_endpoint_new(x, tau, lambda, mu, indices, Rdata, Rreg)
% objective of endpoint representation

    Q1 = [eye(3) -eye(3);-eye(3) eye(3)];
    Q2 = [eye(3) -2*eye(3) eye(3);-2*eye(3) 4*eye(3) -2*eye(3);eye(3) -2*eye(3) eye(3)];

    Q1 = Q1 ./ tau * lambda;
    Q2 = Q2 ./ (tau^3) * mu;
    
% %     Rregnew = fromso3_end_new(Rreg, x);
%     Rregnew=Rreg;

    N = round(length(x)/3);
    BigQ = zeros(3*N,3*N);
     
    fg = zeros(size(BigQ,1),1);
    
    % for v
    for i = 1:N-1
%         v1 = logSO3(Rregnew(:,:,i)'*Rregnew(:,:,i+1));
%         Jl = -eye(3);%leftJinv(v1);
%         Jr = eye(3);%rightJinv(v1);
%         Q1 = [Jl';Jr'];
%         Q1 = (Q1*Q1') ./ tau * lambda;
        BigQ(i*3-2:i*3+3,i*3-2:i*3+3) = BigQ(i*3-2:i*3+3,i*3-2:i*3+3) + Q1;
%         fg(i*3-2:i*3+3,1) = [2.*Jl'*v1;2*Jr'*v1];
    end
    
    for i = 2:N-1
% %         v1 = logSO3(Rregnew(:,:,i)'*Rregnew(:,:,i-1));
%         Jl1 = -eye(3);%leftJinv(v1);
%         Jr1 = eye(3);%rightJinv(v1);
% %         
% %         v2 = logSO3(Rregnew(:,:,i)'*Rregnew(:,:,i+1));
%         Jl2 = -eye(3);%leftJinv(v2);
%         Jr2 = eye(3);%rightJinv(v2);
% %         
%         Q2 = [Jr1';(Jl1+Jl2)';(Jr2)'];
%         Q2 = (Q2*Q2') ./ (tau^3) * mu;
%         
        BigQ(i*3-5:i*3+3,i*3-5:i*3+3) = BigQ(i*3-5:i*3+3,i*3-5:i*3+3) + Q2;
%         fg(i*3-5:i*3+3,1) = [2.*Jr1'*(v1+v2);2.*(Jl1+Jl2)'*(v1+v2);2.*(Jr2)'*(v1+v2)];
    end
    
    %% data terms
    Q3 = zeros(size(BigQ));
    

    Q = BigQ + Q3;
    Q = (Q+Q')/2;
    
%     f = x'*Q*x;
    f = x'*Q*x + fg'*x;
    if nargout > 1
        gradf=(Q+Q')*x+fg;
    end
    if nargout > 2
        Qout = Q;
    end
end

function [Aieq, bieq] = constraint_end1_new_lin(x, Rdata, Rreg, indices)
% endpoint bound constraints
    P0 = [[1 0 0]; ...
          [0 1 0]; ...
          [0 0 1]];
      
    N = round(length(x)/3);
    Aieq = zeros(6*length(indices), length(x));
    bieq = zeros(6*length(indices), 1);
    tol = 0.1.*ones(1,N);%[0.2 0.5 0.5 0.2];
    for i = 1:length(indices)
        ii = indices(i);
        v1 = logSO3(Rdata(:,:,i)');%*Rreg(:,:,ii)
        Jr = eye(3);%rightJinv(v1);
        Aieq(i*6-5:i*6-3, (ii*3-2:ii*3)) = (Jr+0.5.*hat(v1))*P0;
        Aieq(i*6-2:i*6, (ii*3-2:ii*3)) = -(Jr+0.5.*hat(v1))*P0;
        bieq(i*6-5:i*6-3) = -v1 + tol(i);
        bieq(i*6-2:i*6) = v1 + tol(i);
    end
end

function [c, ceq, gradc, gradceq] = constraint_end1_new(x, Rdata, Rreg, indices)
% endpoint bound constraints
    P0 = [[1 0 0]; ...
          [0 1 0]; ...
          [0 0 1]];
      
    N = round(length(x)/3);
%     Aieq = zeros(6*length(indices), length(x));
%     bieq = zeros(6*length(indices), 1);
%     tol = 0.1.*ones(1,N);%[0.2 0.5 0.5 0.2];
%     for i = 1:length(indices)
%         ii = indices(i);
%         v1 = logSO3(Rdata(:,:,i)');%*Rreg(:,:,ii)
%         Jr = eye(3);%rightJinv(v1);
%         Aieq(i*6-5:i*6-3, (ii*3-2:ii*3)) = (Jr+0.5.*hat(v1))*P0;
%         Aieq(i*6-2:i*6, (ii*3-2:ii*3)) = -(Jr+0.5.*hat(v1))*P0;
%         
%         bieq(i*6-5:i*6-3) = -v1 + tol(i);
%         bieq(i*6-2:i*6) = v1 + tol(i);
%     end
    tol = 0.01.*ones(1,N);%[0.2 0.5 0.5 0.2];
    c = zeros(6*length(indices),1);
    for i = 1:length(indices)
        ii = indices(i);
        tmp = logSO3(Rdata(:,:,i)'*expSO3(P0*x(ii*3-2:ii*3)));
        c(i*6-5:i*6-3) = tmp - tol(i);
        c(i*6-2:i*6) = -tmp - tol(i);
    end
    ceq = [];
    if nargout > 2
        gradceq = [];
        der = zeros(6*length(indices),length(x));
        for i = 1:length(indices)
            ii = indices(i);
            v1 = P0*x(ii*3-2:ii*3);
            v2 = logSO3(Rdata(:,:,i)'*expSO3(v1));
            Jr1 = rightJ(v1);
            Jr2 = rightJinv(v2);
            Jr = Jr2*Jr1;
            der(i*6-5:i*6-3, (ii*3-2:ii*3)) = Jr;
            der(i*6-2:i*6, (ii*3-2:ii*3)) = -Jr;
        end
        gradc = der';
    end
    
end

function Rreg = fromso3_end_new(Rreg, x)
% recover SO3 from endpoint representation.
    P0 = [[1 0 0]; ...
          [0 1 0]; ...
          [0 0 1]];
    N = round(length(x)/3);
    for i = 1:N
        Rreg(:,:,i) = expSO3(P0 * x(i*3-2:i*3));
    end
end

function h = hessinterior(x,lambda,A)
    h = A+A';
end
