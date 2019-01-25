function [R,t] = p3p_pst(p,q,K)

    n = size(p,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];
    end
    qn = K\q;
    qn = qn./sqrt(qn(1,:).^2+qn(2,:).^2+qn(3,:).^2);
    
    [B,A,D1,l,C1] = poly4order(p,qn);
    
    t1 = roots(B);
%     pdot = [4*B(1) 3*B(2) 2*B(3) B(4)];
%     pdotdot = [12*B(1) 6*B(2) 2*B(3)];
%     textreme = roots(pdot);
%     
%     for i = 1:numel(textreme)
%         if abs(imag(textreme(i))<1e-6)
%             %% real
%             xx = real(textreme(i));
%             value1 = pdotdot * [xx^2 xx 1]';
%             value2 = B * [xx^4 xx^3 xx^2 xx 1]';
%             if (value1 > 0 && value2 > 0) || (value1 < 0 && value2 < 0)
%                 t1 = [t1;xx];
%             end
%         end
%     end
    
    validt1 = [];
    m = 1;
    for i = 1:numel(t1)
        if abs(imag(t1(i))<1e-6)
            %% real
            if isempty(validt1) || ~isempty(find(abs(validt1-real(t1(i)))<1e-6,1))
                validt1 = [validt1 real(t1(i))];
            end
        else
            continue;
        end
        t = real(t1(i));
        if abs(A(4)+A(5)*t) > 1e-6
            t2 = -(A(3)*t+A(6)*t*t+A(7))/(A(4)+A(5)*t);
        else
            c1 = sqrt(A(1)*t*t+A(2));
            t = [t t];
            t2 = [c1 c2];
        end
        
        for j = 1:numel(t)
            lambda = D1 / sqrt(t(j)*t(j)+C1*C1);
            d0 = lambda;
            d1 = lambda * (l(1)+t(j));
            d2 = lambda * (l(2)+t2(j));
            q3d = qn.*[d0 d1 d2];
            [Rs,ts] = svd_3d23d(p, q3d);
            R(:,:,m) = Rs;
            t(1:3,1,m) = ts;
            m = m + 1;
        end
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