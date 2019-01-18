function [R, t] = p3p_fisch(P, q, K)
    %% 
    n = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];% to homogeneous
    end
    
    %% normalize
    qn = K\q;
    
    %% 3 tetrahedron and 3 angle
    Rab = norm(P(:,1)-P(:,2));%% euclidean transformation, distance is perceived
    Rac = norm(P(:,1)-P(:,3));
    Rbc = norm(P(:,2)-P(:,3));
        
    wab = dot(qn(:,1),qn(:,2))./(norm(qn(:,1))*norm(qn(:,2)));
    wac = dot(qn(:,1),qn(:,3))./(norm(qn(:,1))*norm(qn(:,3)));
    wbc = dot(qn(:,2),qn(:,3))./(norm(qn(:,2))*norm(qn(:,3)));
    
%     %% example
%     Rab = 2*sqrt(3);
%     Rac = 2*sqrt(3);
%     Rbc = 2*sqrt(3);
%     wab = acos(20 / 32);
%     wac = acos(20 / 32);
%     wbc = acos(20 / 32);

    %% K1, K2
    K1 = Rbc^2/(Rac^2);
    K2 = Rbc^2/(Rab^2);
    
    c1 = cos(wbc);
    c2 = c1^2;
    c3 = cos(wab);
    c4 = cos(wac);
    c5 = c4^2;
    c6 = K1*K2;
    
    %% polynomial coefficients
    G4 = (c6-K1-K2)^2-4*c6*c2;
    G3 = 4*(c6-K1-K2)*K2*(1-K1)*c3+4*K1*c1*((c6+K2-K1)*c4+2*K2*c3*c1);
    G2 = (2*K2*(1 - K1)*c3)^2 + 2*(c6 + K1 - K2)*(c6 - K1 - K2) + 4*K1*((K1 - K2)* (c2) + (1 - K2)*K1*(c5) - 2*K2*(1 + K1)*c3*c4*c1);
    G1 = 4*(K1*K2 + K1 - K2)*K2*(1 - K1)*c3 + 4*K1*((c6 - K1+ K2)* c4* c1 + 2*K1*K2* c3* (c5));
    G0 =(K1*K2 + K1 - K2)^2 - 4*(K1^2)*K2 * (c5);
    
    Poly = [G4, G3, G2, G1, G0];
    
    root = roots(Poly);
    
    nvalid = (abs(imag(root))<1e-6);
    
    x = real(root(nvalid));
    pos = x > 0;
    x = x(pos);
    
    turple = [];
    for i = 1:numel(x)
        if ~isempty(turple)
            duplicate_root = turple(:,1) - x(i);
            if ~isempty(find(abs(duplicate_root)<1e-6,1))
                continue;
            end
        end
        
        a = Rab ./ sqrt(x(i)^2-2*x(i)*c3+1);
        b = a*x(i);
        
        m1 = 1-K1;
        p1 = 2*(K1*c4-x(i)*c1);
        q1 = x(i)^2 -  K1;
        m2 = 1;
        p2 = 2*(-x(i)*c1);
        q2 = x(i)^2*(1-K2)+2*x(i)*K2*c3-K2;
        
        den = m1*q2 - m2*q1;
        if abs(den) < 1e-3
            cc1 = c5+((Rac)^2-a^2)/(a^2);
            y(1) = c4 + sqrt(cc1);
            y(2) = c4 - sqrt(cc1);
            
            for j = 1:numel(y)
                if y(j) > 0 && isreal(y(j)) == 1
                    c = y(j)*a;
                    cond1 = b^2+c^2-2*b*c*c1 - Rbc^2;
                    if abs(cond1) < 1e-6
                        turple = [turple;[x(i) a b c y(j)]];
                    end
                end
            end
        else
            y = (p2*q1-p1*q2)/den;
            c = y*a;
            turple = [turple;[x(i) a b c y]];
        end
    end
    
    %%
    R = zeros(3,3,size(turple,1));
    t = zeros(3,1,size(turple,1));
    
    for i = 1:size(turple,1)
        a = turple(i,2);
        b = turple(i,3);
        c = turple(i,4);
%         v1 = qn(:,1);
%         v2 = qn(:,2);
%         v3 = qn(:,3);
        
        [Ropt,topt] = postprocess(P(:,1), P(:,2), P(:,3), ...
                                  a, b, c, Rab, Rac, qn(:,1), qn(:,2), qn(:,3));
        
%         v1 = v1 ./ norm(v1);
%         v2 = v2 ./ norm(v2);
%         v3 = v3 ./ norm(v3);
%         Q = [v1*a v2*b v3*c];
%         [Ropt,topt] = svd_3d23d(P, Q);
        R(:,:,i) = Ropt;
        t(:,:,i) = topt;
    end
    
end

function cosa = coseqn(A,B,C)
    cosa = (A^2+B^2-C^2)/(2*A*B);
end

function varargout = postprocess(A, B, C, a, b, c, ab, ac, ia, ib, ic)
    coslab = coseqn(a,ab,b);
    vecAB = B-A;
    v1 = vecAB./ab;
    Q = A + v1 * coslab * a;
    %% plane 1
    n1 = v1;
    d1 = -n1'*Q;
    plane1 = [n1;d1];
    
    coslac = coseqn(a,ac,c);
    vecAC = C-A;
    v2 = vecAC./ac;
    QQ = A + v2 * coslac * a;
    %% plane 1
    n2 = v2;
    d2 = -n2'*QQ;
    plane2 = [n2;d2];
    
    %% plane 3
    n3 = cross(v2, v1);
    n3 = n3 ./ norm(n3);
    d3 = -n3'*A;
    plane3 = [n3;d3];
    
    %% intersect point
    A1 = [plane1(1:3)';plane2(1:3)';plane3(1:3)'];
    b1 = [-plane1(4);-plane2(4);-plane3(4)];
    R = A1\b1;
    
    %% AR
    ar = norm(R-A);
    rl = sqrt(a^2 - ar^2);
    
    L = R + n3.*rl;
    
    VX = B - A;VX = VX./norm(VX);
    VY = C - A;VY = VY./norm(VY);
    VZ = cross(VX, VY);VZ = VZ./norm(VZ);
    VY = cross(VZ, VX);VY = VY./norm(VY);
    RW = [VX VY VZ];
    
    %% 3 rays
    LA = A - L;vla = LA./a;
    LB = B - L;vlb = LB./b;
    LC = C - L;vlc = LC./c;
    
    %% image
    nia = sqrt(norm(ia)^2);nib = sqrt(norm(ib)^2);nic = sqrt(norm(ic)^2);
    
    iag = L+vla.*nia;
    ibg = L+vlb.*nib;
    icg = L+vlc.*nic;
    
    %% norm image plane
%     nimage = cross(ibg-iag,icg-iag);
%     nimage = nimage./norm(nimage);
    
    %% principle point
    vx = ibg - iag;vx = vx ./ norm(vx);
    vy = icg - iag;vy = vy ./ norm(vy);
    vz = cross(vx,vy);vz = vz ./ norm(vz);
    vy = cross(vz,vx);vy = vy ./ norm(vy);
    Rimg = [vx vy vz];
%     pp = L-nimage;
%     vpp2iag = iag - pp;
%     vpp2ibg = ibg - pp;
%     vpp2icg = icg - pp;
%     
% %     vi1 = vpp2iag./norm(vpp2iag);
% %     vi2 = ia ./ norm(ia);
% %     angle = acos(vi2'*vi1);
% %     
% %     vx = [cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1] * [1;0;0];
% %     vz = -niamge;
% %     vy = cross(vz, vx);
%     AA = [vpp2iag';vpp2ibg';vpp2icg';];
%     bb = [ia(1);ib(1);ic(1)];
%     vx = AA\bb;
%     vy = cross(vz, vx);
    
    R = Rimg*RW';
    t = R'*L;

    varargout{1} = R;
    varargout{2} = t;
end

function [Ropt,topt] = svd_3d23d(ptsrc, ptdst)
    %% SVD 
    ptsrcmean = mean(ptsrc')';
    ptdstmean = mean(ptdst')';

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