function Hn = normalizeHomography(H)
    %% 1
%     tic
    [~,S,~]=svd(H);
    Hn = H./S(2,2);
%     toc
    
    %% 2, analytical solution, use Cardano formula for 3rd order polynomial root function.
%     tic
%     M = H'*H;
%     m11 = M(1,1);m12 = M(1,2);m13 = M(1,3);
%     m22 = M(2,2);m23 = M(2,3);m33 = M(3,3);
%     
%     a2 = -(m11+m22+m33);
%     a1 = m11*m22+m11*m33+m22*m33-(m12*m12+m13*m13+m23*m23);
%     a0 = m12*m12*m33+m13*m13*m22+m23*m23*m11-m11*m22*m33-2*m12*m13*m23;
%     
%     Q = (3*a1-a2*a2)/9;
%     R = (9*a2*a1-27*a0-2*a2*a2*a2)/54;
%     
%     C = sqrt(Q^3+R^2);
%     S1 = (R+C)^(1/3);
%     T = (R-C)^(1/3);
%     
%     l2 = -a2/3 - 0.5*(S1+T) - sqrt(3)/2*(S1-T)*complex(0,1);
%     l2 = real(l2);
%     l = sqrt(l2);
%     toc
    
end