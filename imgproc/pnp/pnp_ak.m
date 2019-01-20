function [R,t] = pnp_ak(P, q, K)
% The implementation of PnP algorithm proposed in
% Linear Pose Estimation from Points or Lines
% Adnan Ansar and Kostas Daniilidis, 
% IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, 
% VOL. 25, NO. 5, MAY 2003
% Author: xiahaa@space.dtu.dk
    n = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];
    end
    qn = K\q;
    qn = qn ./ sqrt(qn(1,:).^2+qn(2,:).^2+qn(3,:).^2);
    
    qij = qn'*qn;
    % dij
    m1 = fix(n*(n-1)*0.5); m2 = fix(n*(n+1)*0.5)+1;
    M = zeros(m1,m2);
    k = 1;
    
    lut = zeros(1,n-1);
    for i = 1:n-1
        Pi = P(:,i);
        mm = i+1:1:n;
        Pjs = P(:,mm);
        nmm = numel(mm);
        d1 = repmat(Pi,1,nmm) - Pjs;
        dij2 = diag(d1'*d1);
        ii1 = k:1:(k+nmm-1);
        M(ii1,end) = -dij2;
        M(ii1,ii1) = -1.*diag(qij(i, mm));
        M(ii1,m1+i) = qij(i,i);
        M(ii1,m1+mm) = diag(diag(qij(mm,mm)));
        lut(i) = k;
        k = k + nmm;
    end
    [~,~,V] = svd(M);
    
    %% null space dimention should be: n+1
    kerM = V(:, n+1:end);
    
    mk1 = fix((n+1)*(n+2)*0.5);
    mk2 = fix(n*n*(n-1)*0.5);
    K = zeros(mk2,mk1);
    
    cnt = 0;
    for i = 1:1:n
        for j = 1:1:n
            for k = j+1:1:n
                cnt = cnt + 1;
            end
        end
    end
    
end