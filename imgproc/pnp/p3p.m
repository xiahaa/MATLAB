function [R, t] = p3p(P, q, K)
    %% 
    n = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];% to homogeneous
    end

    %% normalize
    qn = K\q;

    %% 3 tetrahedron and 3 angle
    ab = norm(P(:,1)-P(:,2));%% euclidean transformation, distance is perceived
    ac = norm(P(:,1)-P(:,3));
    bc = norm(P(:,2)-P(:,3));

    va = qn(:,1)./norm(qn(:,1));
    vb = qn(:,2)./norm(qn(:,2));
    vc = qn(:,3)./norm(qn(:,3));
    
    cosab = dot(va,vb);
    cosac = dot(va,vc);
    cosbc = dot(vb,vc);

    %% input prepare
    p = 2*cosbc;
    q = 2*cosac;
    r = 2*cosab;
    a = bc/ab; a = a*a;
    b = ac/ab; b = b*b;

    I0 = p^2+q^2+r^2-p*q*r-1;

    if abs(I0) < 1e-8
        R = [];
        t = [];
        warning('Degenerate situation~!!!!');
        return;
    end

%     a0 = -2*b+b^2+a^2+1-b*r^2*a+2*b*a-2*a;
%     a1 = -2*b*q*a-2*a^2*q+b*r^2*q*a-2*q+2*b*q+4*a*q+p*b*r+b*r*p*a-b^2*r*p;
%     a2 = q^2+b^2*r^2-b*p^2-q*p*b*r+b^2*p^2-b*r^2*a+2-2*b^2-a*b*r*p*q+2*a^2-4*a-2*q^2*a+q^2*a^2;
%     a3 = -b^2*r*p+b*r*p*a-2*a^2*q+q*p^2*b+2*b*q*a+4*a*q+p*b*r-2*b*q-2*q;
%     a4 = 1-2*a+2*b+b^2-b*p^2+a^2-2*b*a;
% 
%     b0 = b*(p^2*a-p^2+b*p^2+p*q*r-q*a*r*p+a*r^2-r^2-b*r^2)^2;

    %% fast
    a_2 = a*a;
    b2 = b*b;
    q2 = q*q;
    p2 = p*p;
    r2 = r*r;
    ba = a*b;
    baq = ba*q;
    bp = b*p;
    aq = a*q;
    rp = r*p;

    a0 = -2*b+b2+a_2+1-ba*r2+2*ba-2*a;
    a1 = -2*baq-2*a_2*q+baq*r2-2*q+2*b*q+4*aq+b*rp+ba*rp-b2*rp;
    a2 = q2+b2*r2-b*p2-q*bp*r+b2*p2-ba*r2+2-2*b2-baq*rp+2*a_2-4*a-2*q2*a+q2*a_2;
    a3 = -b2*rp+ba*rp-2*a_2*q+q*p2*b+2*baq+4*aq+bp*r-2*b*q-2*q;
    a4 = 1-2*a+2*b+b2-b*p2+a_2-2*ba;
    b0 = b*(p2*a-p2+b*p2+rp*q-aq*rp+a*r2-r2-b*r2)^2;

    if abs(a0) < 1e-8 || abs(b0) < 1e-8
        R = [];
        t = [];
        warning('Degenerate situation~!!!!');
        return;
    end

    rs = roots([a0 a1 a2 a3 a4]);
    valid = abs(imag(rs)) < 1e-8;

    rrs = real(rs(valid));
    
    turple = [];
    for i = 1:numel(rrs)
        x = rrs(i);
        
        if x <= 0
            continue;
        end
        
        if ~isempty(turple)
            duplicate_root = turple(:,1) - x;
            if ~isempty(find(abs(duplicate_root)<1e-6,1))
                continue;
            end
        end
        
        b1 = calcb1(x,a,b,p,q,r);
        if b1 > 0
            y = b1 / b0;
        else
            continue;
        end
        v = x^2+y^2-x*y*r;
        if v <= 0
            continue;
        end
        turple = [turple;[x y v]];
    end
    
    R = zeros(3,3,size(turple,1));
    t = zeros(3,1,size(turple,1));
    for i = 1:size(turple,1)
        x = turple(i,1);y = turple(i,2);v = turple(i,3);
        
        Z = ab / sqrt(v);
        X = x*Z;
        Y = y*Z;
        
        Q = [va.*X vb.*Y vc.*Z];
        [Ropt,topt] = svd_3d23d(P(:,1:3), Q);
        
        R(:,:,i) = Ropt;
        t(:,:,i) = topt;
    end
end

function b1 = calcb1(x,a,b,p,q,r)
    %% todo simplify this
%     b1 = ((1-a-b)*x^2+(q*a-q)*x+1-a+b)*((a^2*r^3+2*b*r^3*a-b*r^5*a-2*a*r^3+r^3+b^2*r^3-2*r^3*b)*x^3 + ...
%           (p*r^2+p*a^2*r^2-2*b*r^3*q*a+2*r^3*b*q-2*r^3*q-2*p*a*r^2-2*p*r^2*b+r^4*p*b+4*a*r^3*q+b*q*a*r^5-2*r^3*a^2*q + ...
%            2*r^2*p*b*a+b^2*r^2*p-r^4*p*b^2)*x^2+(r^3*q^2+r^5*b^2+r*p^2*b^2-4*a*r^3-2*a*r^3*q^2+r^3*q^2*a^2 + ...
%            2*a^2*r^3-2*b^2*r^3-2*p^2*b*r+4*p*a*r^2*q+2*a*p^2*r*b-2*a*r^2*q*b*p-2*p^2*a*r+r*p^2-b*r^5*a+2*p*r^2*b*q + ...
%            r*p^2*a^2-2*p*q*r^2+2*r^3-2*r^2*p*a^2*q-r^4*q*b*p)*x+4*a*r^3*q+p*r^2*q^2+2*p^3*b*a-4*p*a*r^2 - ...
%            2*r^3*b*q-2*p^2*q*r-2*b^2*r^2*p+r^4*p*b+2*p*a^2*r^2-2*r^3*a^2*q-2*p^3*a+p^3*a^2+2*p*r^2+p^3+2*b*r^3*q*a + ... 
%            2*q*p^2*b*r+4*q*a*r*p^2-2*p*a*r^2*q^2-2*p^2*a^2*r*q+p*a^2*r^2*q^2-2*r^3*q-2*p^3*b+p^3*b^2-2*p^2*b*r*q*a);
    r2=r*r;
    r3=r2*r;
    r4=r3*r;    
    r5=r4*r;

    b2=b*b;
    p2=p*p;
    a_2 = a*a;
    q2 = q*q;

    ab = a*b;
    bq = b*q;
    pa = p*a;
    bp = b*p;

    b1 = ((1-a-b)*x^2+(q*a-q)*x+1-a+b)*((a_2*r3+2*r3*ab-r5*ab-2*a*r3+r3+b2*r3-2*r3*b)*x^3 + ...
          (p*r2+p*a_2*r2-2*r3*q*ab+2*r3*bq-2*r3*q-2*pa*r2-2*r2*bp+r4*bp+4*a*r3*q+ab*q*r5-2*r3*a_2*q + ...
           2*r2*p*ab+b2*r2*p-r4*p*b2)*x^2+(r3*q2+r5*b2+r*p2*b2-4*a*r3-2*a*r3*q2+r3*q2*a_2 + ...
           2*a_2*r3-2*b2*r3-2*p2*b*r+4*pa*r2*q+2*p2*r*ab-2*r2*q*ab*p-2*p2*a*r+r*p2-b*r5*a+2*p*r2*bq + ...
           r*p2*a_2-2*p*q*r2+2*r3-2*r2*p*a_2*q-r4*bq*p)*x+4*a*r3*q+p*r2*q2+2*p^3*ab-4*pa*r2 - ...
           2*r3*bq-2*p2*q*r-2*b2*r2*p+r4*p*b+2*p*a_2*r2-2*r3*a_2*q-2*p^3*a+p^3*a_2+2*p*r2+p^3+2*r3*q*ab + ... 
           2*q*p2*b*r+4*q*a*r*p2-2*pa*r2*q2-2*p2*a_2*r*q+p*a_2*r2*q2-2*r3*q-2*p^3*b+p^3*b2-2*p2*ab*r*q);
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