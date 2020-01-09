% the benefit of this observor is that no direct measuremnt of H is
% necessary
function Hest = SL3filter5(Hest, q, p, A, dt)
% Manerikar N, Hua M D, Hamel T. 
% Homography Observer Design on Special Linear Group SL (3) with Application to Optical Flow Estimation[C]
% 2018 European Control Conference (ECC). IEEE, 2018: 1-5.

    e = Hest * q;
    e = e./vecnorm(e);
    Y = e-p;
    Y = Y(:);
    [B1,B2,B3,B4,B5,B6,B7,B8] = sl3basis();
    B1p = B1 * p;B3p = B3 * p;B5p = B5 * p;B7p = B7 * p;
    B2p = B2 * p;B4p = B4 * p;B6p = B6 * p;B8p = B8 * p;
    
    C = zeros(3*size(p,2),8);
    for i = 1:size(p,2)
        pix = eye(3) - p(:,i)*p(:,i)';
        C(3*i-2:3*i,:) = [pix*B1p(:,i) pix*B2p(:,i) pix*B3p(:,i) pix*B4p(:,i) pix*B5p(:,i) pix*B6p(:,i) pix*B7p(:,i) pix*B8p(:,i)];
    end
    k = 0.5;
    delta = -k*C'*Y;
    Delta = delta(1)*B1+delta(2)*B2+delta(3)*B3+delta(4)*B4+delta(5)*B5+delta(6)*B6+delta(7)*B7+delta(8)*B8;

    Hdot = A+inv(Hest)*Delta*Hest;
    Hest = Hest*(eye(3)+(Hdot)*dt);
end

function varargout = sl3basis()
    e1 = [1;0;0];e2 = [0;1;0];e3 = [0;0;1];
    B1 = e1*e2';
    B2 = e2*e1';
    B3 = e2*e3';
    B4 = e3*e2';
    B5 = e3*e1';
    B6 = e1*e3';
    B7 = e1*e1'-1/3*eye(3);
    B8 = e2*e2'-1/3*eye(3);
    varargout{1} = B1;
    varargout{2} = B2;
    varargout{3} = B3;
    varargout{4} = B4;
    varargout{5} = B5;
    varargout{6} = B6;
    varargout{7} = B7;
    varargout{8} = B8;
end