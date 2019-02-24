function MX = mean_Taylor_2nd_adv_recursive3( X, varargin )
%% This function calculates the 2nd order approximation of the mean of a
% bunch of matrices based on the Taylor expansion of the matrix logarithm
% and the definition of mean of a probability density function.

% Input: X is a cell of 4 by 4*n rigid transformation matrices
% Output: M_T1 is the mean of the 1st order approximation of Taylor
% expansion

% Output: MX is the 2nd order approximation of the Taylor expansion
% coder.extrinsic('mean_Taylor_1st_mex');

    %%
    n = size(X,2)/4;
    M_t1 = mean_Taylor_1st( X );

    MX = M_t1;

    iter = 1;
    maxIter = 20;
    if nargin == 2
        weight = varargin{1};
    else
        weight = ones(1,n)./n;
    end
    
%     se3 = zeros(6,n);
%     for k=1:n
%         se3(:,k) = tran2vec(X(:,(k-1)*4+1:k*4));
%     end
    
    
    while iter < maxIter
%         Minv = [MX(1:3,1:3)' -MX(1:3,1:3)'*MX(1:3,4);[0 0 0 1]];
%         AdT = tranAd(Minv);
        LHS = zeros(6,6);
        RHS = zeros(6,1);
        SE3s=MX\X;
%         vs = AdT*se3;
%         vs = vs/AdT;
%         RHS = sum(repmat(weight,6,1).*vs, 1);
        for k=1:n
            SE3 = SE3s(:,(k-1)*4+1:k*4);
            v1 = tran2vec(SE3);
%             v1 = vs(:,k);
            invJ = vec2jacInv(v1);
            RHS = RHS + weight(k).*v1;
            LHS = LHS + weight(k).*invJ;
        end
        xhat = LHS\RHS;
        
%         fori = evalfunc(SO3s, M, weight);
%         % add a line search
%         c = 0.7; tau = 0.9;
%         m = -LHS * xhat;
%         t = -c*m;
%         alpha = 1;
%         liter = 1;
%         maxliter = 50;
%         exit = 0;
%         while liter <= maxliter
%             xhat = xhat.*alpha;
%             Mnew = M*vec2rot(xhat);
%             fnew = evalfunc(SO3s, Mnew, weight);
%             if abs(sum(fori - fnew)) < alpha*abs(sum(t))
%                 exit = 1;
%                 break;
%             else
%                 alpha = alpha * tau;
%             end
%             liter = liter + 1;
%         end
%         if liter > maxliter || exit == 1
            % failure
%             break;
%         end
        MX = MX*vec2tran(xhat);
        if norm(xhat) < 1e-3
            break;
        end
        iter = iter + 1;
    end
%     MX = orthog(MX);
    disp(['iter: ', num2str(iter)]);
end
