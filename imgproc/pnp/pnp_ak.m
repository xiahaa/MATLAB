function [R,t] = pnp_ak(P, q, K, Q)%
% The implementation of PnP algorithm proposed in
% Linear Pose Estimation from Points or Lines
% Adnan Ansar and Kostas Daniilidis, 
% IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, 
% VOL. 25, NO. 5, MAY 2003
% Author: xiahaa@space.dtu.dk
    
    %%
%     P = [1.1698   -0.9294   -1.0804    2.2232; ...
%         -1.6257   -0.7562   -0.9304   -2.0611; ...
%         5.0000    5.0000    5.0000    5.0000];
%     q = [0.3843    0.2999    0.3327    0.4265; ...
%         0.4157    0.0465    0.0266    0.6009; ...
%         1.0000    1.0000    1.0000    1.0000];
% %     
%     Q = [2.3056    1.7991    1.9962    2.5588; ...
%         2.4941    0.2791    0.1593    3.6055; ...
%         6.0000    6.0000    6.0000    6.0000];
%     K = eye(3);

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
        M(ii1,ii1) = -2.*diag(qij(i, mm));
        M(ii1,m1+i) = qij(i,i);
        M(ii1,m1+mm) = diag(diag(qij(mm,mm)));
        lut(i) = k;
        k = k + nmm;
    end
    [~,S,V1] = svd(M);
    
%     %% just for debug
%     dq = sqrt(Q(1,:).^2+Q(2,:).^2+Q(3,:).^2);
%     xx = zeros(m2,1);
%     k = 1;
%     for i = 1:n-1
%         for j = i+1:n
%             xx(k) = dq(i)*dq(j);
%             k = k+1;
%         end
%     end
%     xx(k:k+n-1) = dq.^2;
%     xx(end)=1;
% %     disp(['err1:',num2str(M*xx)]);
%     M*xx
    
    %% null space dimention should be: n+1
    kerM = V1(:, end-(n):end);
    
    mk1 = fix((n+1)*(n+2)*0.5);
    mk2 = fix(n*n*(n-1)*0.5);
    K = zeros(mk2,mk1);
    
    cnt = 0;
    N = n+1;
    
%     %% just for debug
%     lll = kerM\xx;
%     xxx = zeros(mk1,1);
%     k = 1;
%     xxx(k:k+N-1) = (lll.^2)';
%     k = k + N;
%     for i = 1:N-1
%         for j = i+1:N
%             xxx(k) = lll(i)*lll(j);
%             k = k+1;
%         end
%     end
    
    %% use iijk and ijik
    for i = 1:n
        vii = kerM(m1+i,:);
        for j=1:n-1
            if i > j
                idij = lut(j) + i - j - 1;
            elseif i < j
                idij = lut(i) + j - i - 1;
            else 
                %% i == j
                idij = m1+j;
            end
            vij = kerM(idij,:);
            for k=j:n-1
                cnt = cnt+1;
                if i > k
                    idik = lut(k) + i - k - 1;
                elseif i < k
                    idik = lut(i) + k - i - 1;
                else
                    %% i == k
                    idik = m1+k;
                end                
                vik = kerM(idik,:);

                if j == k
                    vjk = kerM(m1+j,:);
                else
                    idjk = lut(j) + k - j - 1;
                    vjk = kerM(idjk,:);
                end
                
                %% laa lab
                vv1 = vii.*vjk-vij.*vik;
                K(cnt,1:N) = vv1;
                
                idab = 1;
                for a = 1:N
                    for b = a+1:N
                        vv2 = (vii(a).*vjk(b) + vii(b).*vjk(a))-(vij(a).*vik(b) + vij(b).*vik(a));
                        K(cnt,N+idab) = vv2;
                        idab = idab+1;
                    end
                end
            end
        end
    end
    
    [~,~,V2]=svd(K);
    sol = V2(:,end);
    l11 = sol(1);
    l = zeros(2, N);
    l(1,1) = sqrt(l11);
    l(2,1) = -l(1,1);
    
    for i = 2:N
        lij = sol(N+i-1);
        l(1,i) = lij/l(1,1);
        l(2,i) = lij/l(2,1);
    end
    
    x1 = zeros(m2,1);
    x2 = zeros(m2,1);
    for i = 1:N
        x1 = x1+kerM(:,i).*l(1,i);
        x2 = x2+kerM(:,i).*l(2,i);
    end
    
    %% if valid
    xsol1 = x1(end-n:end);
    xsol2 = x2(end-n:end);
    
    xii1 = xsol1(1:end-1)./xsol1(end);
    xii2 = xsol2(1:end-1)./xsol2(end);
    
    ti = [];
    if isempty(find(xii1<=0,1))
        % valid;
        ti = sqrt(xii1);
    elseif isempty(find(xii2<=0,1))
        % valid;
        ti = sqrt(xii2);
    end    
    
    %% potential an error when xsol(end) == 0
    if ~isempty(ti) && isempty(find(isnan(ti),1))
        Qr = qn.*ti';
        [R,t] = svd_3d23d(P, Qr);
    else
        R = [];
        t = [];
    end
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