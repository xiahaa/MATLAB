function M = FNS_iterative2(SO3s,M_1st,varargin)
% My implementation of iterative optimizing the mean using manfiold
% optimization technique. It works directly with the definition of the mean
% and applies optimization directly on SO3.
    M = M_1st;
    iter = 1;
    maxIter = 50;
    if nargin == 2
        weight = ones(1,size(SO3s,3))./size(SO3s,3);
    else
        weight = varargin{1};
    end
    LHS = zeros(3,3);
    RHS = zeros(3,1);
    for k=1:size(SO3s,3)
        v1 = rot2vec(SO3s(:,:,k));
        RHS = weight(k).*v1;
        LHS = LHS + weight(k).*(eye(3)-0.5.*skewm(v1));
    end
    xhat = LHS\RHS;
    Mnew = vec2rot(xhat);
    
    fori = evalfunc(SO3s, M, weight);
%     disp(fori);
    fnew = evalfunc(SO3s, Mnew, weight);
%     disp(fnew);
    
    if sum(abs(fnew)) < sum(abs(fori))
        M = Mnew;
    end
    
%     while iter < maxIter
%         Minv = (M)';
%         LHS = zeros(3,3);
%         RHS = zeros(3,1);
%         for k=1:size(SO3s,3)
%             v1 = rot2vec(Minv*SO3s(:,:,k));
%             RHS = RHS + weight(k).*v1;
%             LHS = LHS + weight(k).*(eye(3)-0.5.*skewm(v1));
%         end
%         xhat = LHS\RHS;
%         
%         fori = evalfunc(SO3s, M, weight);
%         % add a line search
%         c = 0.7; tau = 0.9;
%         m = -LHS * xhat;
%         t = -c*m;
%         alpha = 1;
%         liter = 1;
%         maxliter = 50;
%         while liter <= maxliter
%             xhat = xhat.*alpha;
%             Mnew = M*vec2rot(xhat);
%             fnew = evalfunc(SO3s, Mnew, weight);
%             gain = norm(fnew-fori);
%             if gain < 0 || gain < alpha*norm(t)
%                 break;
%             else
%                 alpha = alpha * tau;
%             end
%             liter = liter + 1;
%         end
%         
%         if liter > maxliter || gain < 0
%             % failure
%             break;
%         end
%         
%         M = M*vec2rot(xhat);
%         if norm(xhat) < 1e-6
%             break;
%         end
%         iter = iter + 1;
%     end
%     fnew = evalfunc(SO3s, M, weight);
%     disp(fnew);
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