function dxi = seq_sol(xi, v, indices, tau, lambda, miu, N, id, varargin)
    if nargin > 8
    	Rreg = varargin{1};
    end
    if nargin > 9
    	options = varargin{2};
    end

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
    if c2~=0
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
%
%     dxi = quadprog(2.*LHS,2*RHS',Aineq,bineq,[],[],[],[],[],options);
end
