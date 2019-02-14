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
        
    wab = acos(dot(qn(:,1),qn(:,2))./(norm(qn(:,1))*norm(qn(:,2))));
    wac = acos(dot(qn(:,1),qn(:,3))./(norm(qn(:,1))*norm(qn(:,3))));
    wbc = acos(dot(qn(:,2),qn(:,3))./(norm(qn(:,2))*norm(qn(:,3))));
    
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
    
    tuple = [];
    for i = 1:numel(x)
        if ~isempty(tuple)
            duplicate_root = tuple(:,1) - x(i);
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
                        tuple = [tuple;[x(i) a b c y(j)]];
                    end
                end
            end
        else
            y = (p2*q1-p1*q2)/den;
            c = y*a;
            tuple = [tuple;[x(i) a b c y]];
        end
    end
    
    %%
    R = zeros(3,3,size(tuple,1));
    t = zeros(3,1,size(tuple,1));
    
    for i = 1:size(tuple,1)
        a = tuple(i,2);
        b = tuple(i,3);
        c = tuple(i,4);
        v1 = qn(:,1);
        v2 = qn(:,2);
        v3 = qn(:,3);
        
        v1 = v1 ./ norm(v1);
        v2 = v2 ./ norm(v2);
        v3 = v3 ./ norm(v3);
        Q = [v1*a v2*b v3*c];
        [Ropt,topt] = svd_3d23d(P(:,1:3), Q);
        R(:,:,i) = Ropt;
        t(:,:,i) = topt;
    end
    
end

function cosa = coseqn(A,B,C)
    cosa = (A^2+B^2-C^2)/(2*A*B);
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