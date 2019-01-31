function [R, t] = epnp_re(p,q,K)
% Implementation of EPnP algorithm.
    n = size(p,2);
    %% 1st step, define the virtual point
    cw = [1 0 0;0 1 0;0 0 1;0 0 0]';
    Cw = [1 0 0 0;0 1 0 0;0 0 1 0;1 1 1 1];
    ph = [p;ones(1,n)];
    alpha = Cw\ph;
    
    %% create M
    fu = K(1,1);fv = K(2,2);
    cu = K(1,3);cv = K(2,3);
    
    M = zeros(2*n,12);
    for i = 1:n
        cu_ui = cu-q(1,i);
        M(i*2-1, :) = [alpha(1,i)*fu 0 alpha(1,i)*(cu_ui) ...
                       alpha(2,i)*fu 0 alpha(2,i)*(cu_ui) ...
                       alpha(3,i)*fu 0 alpha(3,i)*(cu_ui) ...
                       alpha(4,i)*fu 0 alpha(4,i)*(cu_ui)];
        cv_vi = cv-q(2,i);
        M(i*2, :) = [0 alpha(1,i)*fv alpha(1,i)*(cv_vi) ...
                     0 alpha(2,i)*fv alpha(2,i)*(cv_vi) ...
                     0 alpha(3,i)*fv alpha(3,i)*(cv_vi) ...
                     0 alpha(4,i)*fv alpha(4,i)*(cv_vi)];
    end
    
    %% deomposition, 
    [ev,ee] = eig(M'*M);
    V = ev(:,4:-1:1);
    
    %% l1*V1+l2*V2+l3*V3+l4*V4 = solution
    v1 = V(:,end);%% >= 6 points
    v2 = V(:,end-1);%% 5 points
    v3 = V(:,end-2);
    v4 = V(:,end-3);%% 4 points
    
    cwdist = calc_cwdist(cw);
    
    %% solution 1
    [x1, ~] = case1(v1, cwdist);
    qc1 = pc_in_cam(x1, p, alpha);
    [R1(:,:,1),t1(:,1)] = svd_3d23d(p, qc1);
    
    %% solution 2
    [x2, ~] = case2(v1, v2, cwdist);
    qc2 = pc_in_cam(x2, p, alpha);
    
    %% solution 3
    [x3, ~] = case3(v1, v2, v3, cwdist);
    qc3 = pc_in_cam(x3, p, alpha);

    %% solution 4
    [x4, ~] = case4(v1, v2, v3, v4, cwdist);
    qc4 = pc_in_cam(x4, p, alpha);
    
    %% 3D-3D pose estimation
    R1 = zeros(3,3,4);
    t1 = zeros(3,4);
    
    [R1(:,:,2),t1(:,2)] = svd_3d23d(p, qc2);
    [R1(:,:,3),t1(:,3)] = svd_3d23d(p, qc3);
    [R1(:,:,4),t1(:,4)] = svd_3d23d(p, qc4);
    
    %% find the best solution
    errs = zeros(1,4);
    for i = 1:4
        errs(i)=calc_reproj_err(p,q,R1(:,:,i),t1(:,i),K);
    end
    [~,minid] = min(errs);
    R = R1(:,:,minid);
    t = t1(:,minid);
end

function err = calc_reproj_err(p,q,R,t,K)
    pc = R*p+repmat(t,1,size(p,2));
    pc = pc./pc(3,:);
    pq = K*pc;
    err = pq(1:2,:) - q(1:2,:);%2xn
    err = sum(diag(err'*err));
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

function [pc, x] = pc_in_cam(x, p, alpha)
    %% firstly, obtain beta
    v1 = x(1:3);
    v2 = x(4:6);
    v3 = x(7:9);
    v4 = x(10:12);
    
    pc = v1.*alpha(1,:)+v2.*alpha(2,:)+v3.*alpha(3,:)+v4.*alpha(4,:);
    pcm = mean(pc,2);
    cdist2c = pc-repmat(pcm,1,size(p,2));
    cdist2c = sqrt(diag(cdist2c'*cdist2c));
    
    pm = mean(p,2);
    wdist2c = p-repmat(pm,1,size(p,2));
    wdist2c = sqrt(diag(wdist2c'*wdist2c));
    
    beta = (wdist2c)'*(cdist2c)/(cdist2c'*cdist2c);
    x = beta .* x;
    pc = pc.*beta;
    
    if ~isempty(find(pc(3,:)<0,1))
        pc = pc.*-1;
        x = x.*(-1);
    end
end

function dists = calc_cwdist(cw)
    ii = [1 1 1 2 2 3];
    jj = [2 3 4 3 4 4];
    b = cw(:,ii)-cw(:,jj);
    dists = b(1,:).^2+b(2,:).^2+b(3,:).^2;
end

function [x, beta] = case1(v, cwdist)
    v1 = zeros(3,4);
    v1(:,1) = v(1:3);
    v1(:,2) = v(4:6);
    v1(:,3) = v(7:9);
    v1(:,4) = v(10:12);
    %% Ax = b -> x = Ab/(AA) 
    iis = [1 1 1 2 2 3];
    jjs = [2 3 4 3 4 4];

    A = zeros(6,1);
    for i = 1:6
        ii = iis(i);
        jj = jjs(i);
        A(i,1) = v1(:,ii)'*v1(:,ii)+v1(:,jj)'*v1(:,jj)-2.*v1(:,ii)'*v1(:,jj);
    end
    b = cwdist';
    
    xx = A\b;
    
    beta(1) = sqrt(xx(1));
    x = beta(1).*v;
end

function [x, beta] = case2(v, vv, cwdist)
    v1 = zeros(3,4);
    v1(:,1) = v(1:3);
    v1(:,2) = v(4:6);
    v1(:,3) = v(7:9);
    v1(:,4) = v(10:12);
    
    v2 = zeros(3,4);
    v2(:,1) = vv(1:3);
    v2(:,2) = vv(4:6);
    v2(:,3) = vv(7:9);
    v2(:,4) = vv(10:12);
    
    A = zeros(6,3);
    
    iis = [1 1 1 2 2 3];
    jjs = [2 3 4 3 4 4];
    
    %% linearization
    for i = 1:6
        ii = iis(i);
        jj = jjs(i);
        A(i,1) = v1(:,ii)'*v1(:,ii)+v1(:,jj)'*v1(:,jj)-2.*v1(:,ii)'*v1(:,jj);
        A(i,2) = 2.*v1(:,ii)'*v2(:,ii)+2.*v1(:,jj)'*v2(:,jj)-4.*v1(:,ii)'*v2(:,jj);
        A(i,3) = v2(:,ii)'*v2(:,ii)+v2(:,jj)'*v2(:,jj)-2.*v2(:,ii)'*v2(:,jj);
    end
    
    b = cwdist';
    
    xx = A\b;
    
    beta(1) = sqrt(xx(1));
    beta(2) = xx(2)/beta(1);
    
    x = beta(1).*v+beta(2).*vv;
end

function [x, beta] = case3(v, vv, vvv, cwdist)
    v1 = zeros(3,4);
    v1(:,1) = v(1:3);
    v1(:,2) = v(4:6);
    v1(:,3) = v(7:9);
    v1(:,4) = v(10:12);
    
    v2 = zeros(3,4);
    v2(:,1) = vv(1:3);
    v2(:,2) = vv(4:6);
    v2(:,3) = vv(7:9);
    v2(:,4) = vv(10:12);
    
    v3 = zeros(3,4);
    v3(:,1) = vvv(1:3);
    v3(:,2) = vvv(4:6);
    v3(:,3) = vvv(7:9);
    v3(:,4) = vvv(10:12);
    
    A = zeros(6,6);
    
    iis = [1 1 1 2 2 3];
    jjs = [2 3 4 3 4 4];
    
    %% linearization
    for i = 1:6
        ii = iis(i);
        jj = jjs(i);
        A(i,1) = v1(:,ii)'*v1(:,ii)+v1(:,jj)'*v1(:,jj)-2.*v1(:,ii)'*v1(:,jj);
        A(i,2) = 2.*v1(:,ii)'*v2(:,ii)+2.*v1(:,jj)'*v2(:,jj)-4.*v1(:,ii)'*v2(:,jj);
        A(i,3) = 2.*v1(:,ii)'*v3(:,ii)+2.*v1(:,jj)'*v3(:,jj)-4.*v1(:,ii)'*v3(:,jj);
        A(i,4) = v2(:,ii)'*v2(:,ii)+v2(:,jj)'*v2(:,jj)-2.*v2(:,ii)'*v2(:,jj);
        A(i,5) = 2.*v2(:,ii)'*v3(:,ii)+2.*v2(:,jj)'*v3(:,jj)-4.*v2(:,ii)'*v3(:,jj);
        A(i,6) = v3(:,ii)'*v3(:,ii)+v3(:,jj)'*v3(:,jj)-2.*v3(:,ii)'*v3(:,jj);
    end
    
    b = cwdist';
    
    xx = A\b;
    
    beta(1) = sqrt(xx(1));
    beta(2) = xx(2)/beta(1);
    
    x = beta(1).*v+beta(2).*vv;
end

function [x, beta] = case4(v, vv, vvv, vvvv, cwdist)
    v1 = zeros(3,4);
    v1(:,1) = v(1:3);
    v1(:,2) = v(4:6);
    v1(:,3) = v(7:9);
    v1(:,4) = v(10:12);
    
    v2 = zeros(3,4);
    v2(:,1) = vv(1:3);
    v2(:,2) = vv(4:6);
    v2(:,3) = vv(7:9);
    v2(:,4) = vv(10:12);
    
    v3 = zeros(3,4);
    v3(:,1) = vvv(1:3);
    v3(:,2) = vvv(4:6);
    v3(:,3) = vvv(7:9);
    v3(:,4) = vvv(10:12);
    
    v4 = zeros(3,4);
    v4(:,1) = vvvv(1:3);
    v4(:,2) = vvvv(4:6);
    v4(:,3) = vvvv(7:9);
    v4(:,4) = vvvv(10:12);
    
    A = zeros(6,10);
    
    iis = [1 1 1 2 2 3];
    jjs = [2 3 4 3 4 4];
    
    %% linearization
    for i = 1:6
        ii = iis(i);
        jj = jjs(i);
        A(i,1) =    v1(:,ii)'*v1(:,ii)+v1(:,jj)'*v1(:,jj)-2.*v1(:,ii)'*v1(:,jj);
        A(i,2) = 2.*v1(:,ii)'*v2(:,ii)+2.*v1(:,jj)'*v2(:,jj)-4.*v1(:,ii)'*v2(:,jj);
        A(i,3) = 2.*v1(:,ii)'*v3(:,ii)+2.*v1(:,jj)'*v3(:,jj)-4.*v1(:,ii)'*v3(:,jj);
        A(i,4) = 2.*v1(:,ii)'*v4(:,ii)+2.*v1(:,jj)'*v4(:,jj)-4.*v1(:,ii)'*v4(:,jj);

        A(i,5) =    v2(:,ii)'*v2(:,ii)+v2(:,jj)'*v2(:,jj)-2.*v2(:,ii)'*v2(:,jj);
        A(i,6) = 2.*v2(:,ii)'*v3(:,ii)+2.*v2(:,jj)'*v3(:,jj)-4.*v2(:,ii)'*v3(:,jj);
        A(i,7) = 2.*v2(:,ii)'*v4(:,ii)+2.*v2(:,jj)'*v4(:,jj)-4.*v2(:,ii)'*v4(:,jj);

        A(i,8) =    v3(:,ii)'*v3(:,ii)+v3(:,jj)'*v3(:,jj)-2.*v3(:,ii)'*v3(:,jj);
        A(i,9) = 2.*v3(:,ii)'*v4(:,ii)+2.*v3(:,jj)'*v4(:,jj)-4.*v3(:,ii)'*v4(:,jj);

        A(i,10) =    v4(:,ii)'*v4(:,ii)+v4(:,jj)'*v4(:,jj)-2.*v4(:,ii)'*v4(:,jj);    
    end
    
    %% relinearization
    V = null(A);
    l1 = V(:,1);l2 = V(:,2);l3 = V(:,3);l4 = V(:,4);
    
    AA = zeros(10,10);
    tuple = [[1 5 2 2]; ...
             [1 6 2 3]; ...
             [1 7 2 4]; ...
             [1 8 3 3]; ...
             [1 9 3 4]; ...
             [1 10 4 4]; ...
             [5 8 6 6]; ...
             [5 9 6 7]; ...
             [5 10 7 7]; ...
             [8 10 9 9]];

    for i = 1:10
        AA(i,1) = l1(tuple(i,1))*l1(tuple(i,2))-l1(tuple(i,3))*l1(tuple(i,4));% 11
        AA(i,2) = l1(tuple(i,1))*l2(tuple(i,2))+l1(tuple(i,2))*l2(tuple(i,1)) - (l1(tuple(i,3))*l2(tuple(i,4))+l1(tuple(i,4))*l2(tuple(i,3)));%12
        AA(i,3) = l1(tuple(i,1))*l3(tuple(i,2))+l1(tuple(i,2))*l3(tuple(i,1)) - (l1(tuple(i,3))*l3(tuple(i,4))+l1(tuple(i,4))*l3(tuple(i,3)));%13
        AA(i,4) = l1(tuple(i,1))*l4(tuple(i,2))+l1(tuple(i,2))*l4(tuple(i,1)) - (l1(tuple(i,3))*l4(tuple(i,4))+l1(tuple(i,4))*l4(tuple(i,3)));%14

        AA(i,5) = l2(tuple(i,1))*l2(tuple(i,2))-l2(tuple(i,3))*l2(tuple(i,4));% 22
        AA(i,6) = l2(tuple(i,1))*l3(tuple(i,2))+l2(tuple(i,2))*l3(tuple(i,1)) - (l2(tuple(i,3))*l3(tuple(i,4))+l2(tuple(i,4))*l3(tuple(i,3)));%23
        AA(i,7) = l2(tuple(i,1))*l4(tuple(i,2))+l2(tuple(i,2))*l4(tuple(i,1)) - (l2(tuple(i,3))*l4(tuple(i,4))+l2(tuple(i,4))*l4(tuple(i,3)));%24

        AA(i,8) = l3(tuple(i,1))*l3(tuple(i,2))-l3(tuple(i,3))*l3(tuple(i,4));% 33
        AA(i,9) = l3(tuple(i,1))*l4(tuple(i,2))+l3(tuple(i,2))*l4(tuple(i,1)) - (l3(tuple(i,3))*l4(tuple(i,4))+l3(tuple(i,4))*l4(tuple(i,3)));%34

        AA(i,10) = l4(tuple(i,1))*l4(tuple(i,2))-l4(tuple(i,3))*l4(tuple(i,4));% 44
    end
    
    [~,~,V1] = svd(AA);
    sol1 = V1(:,end);
    lambda1 = sqrt(sol1(1));    
    lambda2 = sol1(2)/lambda1;
    lambda3 = sol1(3)/lambda1;
    lambda4 = sol1(4)/lambda1;
    
    b = cwdist';
    x0 = A\b;
    xx = x0+lambda1.*l1 + lambda2.*l2 + lambda3.*l3 + lambda4.*l4;
    
    beta(1) = sqrt(xx(1));%%11
    beta(2) = xx(2)/beta(1);%%12
    beta(3) = xx(3)/beta(1);%%13
    beta(4) = xx(3)/beta(1);%%14
    
    x = beta(1).*v+beta(2).*vv+beta(3).*vvv+beta(4).*vvvv;
end