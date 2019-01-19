function [R,t] = p3p_kneip(P, q, K)
%% TODO, I should add some information about this code.
    %% 
    n = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];% to homogeneous
    end

    %% normalize
    qn = K\q;

    %% tau frame
    f1 = qn(:,1);f2 = qn(:,2);f3 = qn(:,3);
    f1 = f1 ./ norm(f1);f2 = f2./norm(f2);f3 = f3./norm(f3);
    t1 = f1;t2 = f2;
    t3 = cross(t1,t2);t3 = t3./norm(t3);
    t2 = cross(t3,t1);t2 = t2./norm(t2);
    T = [t1';t2';t3'];
    f3tau = T*f3;
    P1 = P(:,1);P2 = P(:,2);P3 = P(:,3);
    if ( f3tau(3,1) > 0 )
        f1 = qn(:,2);
        f2 = qn(:,1);
        f3 = qn(:,3);

        f1 = f1 ./ norm(f1);f2 = f2./norm(f2);f3 = f3./norm(f3);
        t1 = f1;t2 = f2;
        t3 = cross(t1,t2);t3 = t3./norm(t3);
        t2 = cross(t3,t1);t2 = t2./norm(t2);

        T = [t1';t2';t3'];

        f3tau = T*f3;
        
        P1 = P(:,2);
        P2 = P(:,1);
        P3 = P(:,3);
        
    end
    
    %% eta frame
    nx = P2-P1;ny = P3-P1;
    nx = nx ./ norm(nx);ny = ny./norm(ny);
    nz = cross(nx,ny);nz = nz./norm(nz);
    if abs(norm(nz)) < 1e-8
        R=[];t=[];return;
    end
    ny = cross(nz,nx);ny = ny./norm(ny);
    N = [nx ny nz]';
    P3eta = N *(P3-P1);
    
    p1 = P3eta(1);p2 = P3eta(2);
    d12 = norm(P2-P1);
    cosbeta = f1'*f2;
    cosbeta2 = cosbeta * cosbeta;
    if abs(cosbeta2-1)<1e-8
        R=[];t=[];return;
    end
    b = sign(cosbeta)*sqrt(1/(1-cosbeta2)-1);
    psi1 = f3tau(1)/f3tau(3);psi2 = f3tau(2)/f3tau(3);
    
    %% precompute
    psi12 = psi1*psi1;
    psi22 = psi2*psi2;
    psi1_2 = psi1*psi2;
    p22 = p2*p2; p23 = p22*p2; p24 = p23*p2;
    d12b = d12*b;
    p12 = p1*p1;
    d122 = d12*d12;
    b2 = b*b;

    a4 = -psi22*p24 - psi12*p24 - p24;

    a3 = 2*p23*d12b + 2*psi22*p23*d12b - 2*psi1_2*p23*d12;

    a2 = -psi22*p12*p22 - psi22*p22*d122*b2 - psi22*p22*d122 + psi22*p24 ...
        +psi12*p24 + 2*p1*p22*d12 + 2*psi1_2*p1*p22*d12b ...
        -psi12*p12*p22 + 2*psi22*p1*p22*d12 - p22*d122*b2 - 2*p12*p22;

    a1 = 2*p12*p2*d12b + 2*psi1_2*p23*d12 - 2*psi22*p23*d12b - 2*p1*p2*d122*b;

    a0 = -2*psi1_2*p1*p22*d12b + psi22*p22*d122 + 2*p1^3*d12 - p12*d122 + psi22*p12*p22 ...
        - p1^4 - 2*psi22*p1*p22*d12 + psi12*p12*p22 + psi22*p22*d122*b2;
    
    
    rs = roots([a4 a3 a2 a1 a0]);
    valid = abs(imag(rs)) < 0.01;
    rss = real(rs(valid));
    
    psi1_psi2 = psi1/psi2;
    
    k = 1;
    R = [];
    t = [];
    for i = 1:numel(rss)
        num = psi1_psi2*p1-d12b;
        den = d12-p1;
        if abs(rss(i)) <= 1
            %% cotalpha
            num = num + rss(i)*p2;
            den = den + psi1_psi2*rss(i)*p2;
            cotalpha = num/den;
            
            %% sinalpha, cosalpha, sintheta
            salpha = sqrt(1 / (1+cotalpha*cotalpha));
            calpha = salpha * cotalpha;
            
            stheta = sqrt(1-rss(i)*rss(i));
            if f3tau(3) > 0
                stheta = -1*stheta;
            end
            
            %% Ceta, Q
            cons1 = salpha*b+calpha;
            Ceta = [d12*calpha*cons1; ...
                    d12*salpha*rss(i)*cons1; ...
                    d12*salpha*stheta*cons1];
            Q = [-calpha -salpha*rss(i) -salpha*stheta; ...
                  salpha -calpha*rss(i) -calpha*stheta; ...
                       0       -stheta           rss(i)];
                   
           %% R t
           Cw = P1 + N'*Ceta;
           Rw = N'*Q'*T;
           
           R(:,:,k) = Rw';
           t(:,:,k) = -Rw'*Cw; 
           k = k + 1;
            
        end
    end

    
end
