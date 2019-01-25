function [R,t] = pnp_long(p,q,K)
%% implementation of PnP algorithm proposed in 
%% Linear N-Point Camera Pose Determination
%% Long Quan and Zhongdan Lan, IEEE PAMI, 1999.
    n = size(p,2);
    
    if n < 4
        warning("only work for 4 points and more!!!!");
        R = [];t = [];
    end
    
    if size(q,1) == 2
        q = [q;ones(1,n)];
    end
    qn = K\q;
    qn = qn./sqrt(qn(1,:).^2+qn(2,:).^2+qn(3,:).^2);
    
    depth = zeros(1,n);
    for i = 1:n
        selset = 1:n;
        selset(i) = [];
        seli = nchoosek(selset,2);
        A = zeros(size(seli,1),5);
        for j = 1:size(seli,1)
            id = [i,seli(j,:)];
            A(j,:) = form_polynomial(p(:,id),qn(:,id));
        end
        %% sol
        sol = solve_x(A);
        %% recover x
        x = sol(2)/sol(1);
        if x > 0
            x = sqrt(x);
        else
            %% error
            R = [];
            t = [];
            return;
        end        
        depth(i) = x;
    end
   
    q3d = qn.*depth;
    [R,t] = svd_3d23d(p, q3d);
end

function [Ropt,topt] = svd_3d23d(ptsrc, ptdst)
%% SVD 

ptsrcmean = mean(ptsrc,2);
ptdstmean = mean(ptdst,2);

ptsrcrefine = ptsrc - repmat(ptsrcmean, 1, size(ptsrc,2));
ptdstrefine = ptdst - repmat(ptdstmean, 1, size(ptsrc,2));

Y = ptdstrefine';
X = ptsrcrefine;
S = X*Y;
[U,~,V] = svd(S);

D = V*U';
if det(D) < 0
    Ropt = V*[1 0 0;0 1 0;0 0 -1]*U';
else
    Ropt = V*U';
end
topt = ptdstmean - Ropt * ptsrcmean;


end

function sol = solve_x(A)
    if size(A,1) > size(A,2)
        %% unique solution
        nvid = 5;
    else
        nvid = 4:5;
    end
    
    [~,~,V] = svd(A);
    
    if numel(nvid) == 1
        sol = V(:,nvid);
        return;
    end
    
    %% null space
    v4 = V(:,4);
    v5 = V(:,5);
    
    %% 
    tuple = [[4,2,3,3]; ...
             [4,1,3,2]; ...
             [4,0,3,1]; ...
             [4,0,2,2]; ...
             [3,1,2,2]; ...
             [3,0,2,1]; ...
             [2,0,1,1]];
    tuple = tuple + 1;
    A = zeros(7,3);
    A(:,1) = v4(tuple(:,1)).*v4(tuple(:,2))-v4(tuple(:,3)).*v4(tuple(:,4));
    A(:,2) = v4(tuple(:,1)).*v5(tuple(:,2))+v5(tuple(:,1)).*v4(tuple(:,2)) - (v4(tuple(:,3)).*v5(tuple(:,4))+v5(tuple(:,3)).*v4(tuple(:,4)));
    A(:,3) = v5(tuple(:,1)).*v5(tuple(:,2))-v5(tuple(:,3)).*v5(tuple(:,4));
    
    [~,~,V2] = svd(A);
    y3 = V2(:,end);
    
    %% average
    s1 = (y3(1)/y3(2) + y3(2) / y3(3)) * 0.5;
    lambda = 1 / (v4(1) + v5(1)/s1);
    rho = lambda / s1;
    
    sol = v4.*lambda + v5.*rho;
    
end