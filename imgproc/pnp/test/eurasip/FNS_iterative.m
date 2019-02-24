function M = FNS_iterative(SO3s,M_1st,varargin)
% My implementation of iterative optimizing the mean using manfiold
% optimization technique. It works directly with the definition of the mean
% and applies optimization directly on SO3.
    M = M_1st;
    iter = 1;
    maxIter = 100;
    if nargin == 2
        weight = ones(1,size(SO3s,3))./size(SO3s,3);
    else
        weight = varargin{1};
    end
    fori = evalfunc(SO3s, M, weight);
    disp(fori);
    while iter < maxIter
        Minv = (M)';
        LHS = zeros(3,3);
        RHS = zeros(3,1);
        for k=1:size(SO3s,3)
            v1 = rot2vec(Minv*SO3s(:,:,k));
            invJ = vec2jacInv(v1);
            RHS = RHS + weight(k).*v1;
            LHS = LHS + weight(k).*invJ;
        end
        xhat = LHS\RHS;
        
        fori = evalfunc(SO3s, M, weight);
        % add a line search
        c = 0.7; tau = 0.9;
        m = -LHS * xhat;
        t = -c*m;
        alpha = 1;
        liter = 1;
        maxliter = 50;
        exit = 0;
        while liter <= maxliter
            xhat = xhat.*alpha;
            Mnew = M*vec2rot(xhat);
            fnew = evalfunc(SO3s, Mnew, weight);
            if abs(sum(fori - fnew)) < alpha*abs(sum(t))
                exit = 1;
                break;
            else
                alpha = alpha * tau;
            end
            liter = liter + 1;
        end
        
        if liter > maxliter || exit == 1
            % failure
            break;
        end
        
        M = M*vec2rot(xhat);
        if norm(xhat) < 1e-3
            break;
        end
        iter = iter + 1;
    end
    fnew = evalfunc(SO3s, M, weight);
    disp(fnew);
end

function f = evalfunc(SO3s, M, weight)
    Minv = (M)';
    f = zeros(3,1);
    for k=1:size(SO3s,3)
        v1 = rot2vec(Minv*SO3s(:,:,k));
        f = f + weight(k).*v1;
    end
end

%% plan, i think I need to add a line search in order to make sure the obtained result is definitely a descending direction