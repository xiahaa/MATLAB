function varargout = homo_decom_malis(varargin)
%% decompose homography matrix using olivier algo.
    H = varargin{1};
    format long;
    
    S = H'*H-eye(3);
    
    s11 = S(1,1);s12 = S(1,2);s13 = S(1,3);s22 = S(2,2);s23 = S(2,3);s33 = S(3,3);
%     s21 = S(2,1);
%     s31 = S(3,1);s32 = S(3,2);
    
    Ms11 = s23*s23 - s22*s33;
    Ms22 = s13*s13 - s11*s33;
    Ms33 = s12*s12 - s11*s22;
    
    Ms12 = s23*s13 - s12*s33;
    Ms13 = s22*s13 - s12*s23;
    Ms23 = s12*s13 - s11*s23;
    
    vscalar = 2*sqrt(1+trace(S)-Ms11-Ms22-Ms33);
    tenorm = sqrt(2+ trace(S)-vscalar);

    %% intermediate varibales
    if abs(s22) > 1e-6
        epsilon13 = signc(Ms13);
        nea = [s12+sqrt(Ms33);s22;s23-epsilon13*sqrt(Ms11)];
        neb = [s12-sqrt(Ms33);s22;s23+epsilon13*sqrt(Ms11)];
        
        na = nea./norm(nea);
        nb = neb./norm(neb);
        
        epsilons22 = signc(s22);
        rho = sqrt(tenorm*tenorm+2*vscalar);
        tas = tenorm*0.5*(epsilons22*rho*nb-tenorm*na);
        tbs = tenorm*0.5*(epsilons22*rho*na-tenorm*nb);
        
    elseif abs(s11) > 1e-6
        epsilon23 = signc(Ms23);
        nea = [s11;s12+sqrt(Ms33);s13+epsilon23*sqrt(Ms22)];
        neb = [s11;s12-sqrt(Ms33);s13-epsilon23*sqrt(Ms22)];
        
        na = nea./norm(nea);
        nb = neb./norm(neb);
        
        epsilons11 = signc(s11);
        rho = sqrt(tenorm*tenorm+2*vscalar);
        tas = tenorm*0.5*(epsilons11*rho*nb-tenorm*na);
        tbs = tenorm*0.5*(epsilons11*rho*na-tenorm*nb);
        
    elseif abs(s33) > 1e-6
        epsilon12 = signc(Ms12);
        nea = [s13+epsilon12*sqrt(Ms22);s23+sqrt(Ms11);s33];
        neb = [s13-epsilon12*sqrt(Ms22);s23-sqrt(Ms11);s33];
        
        na = nea./norm(nea);
        nb = neb./norm(neb);
        
        epsilons33 = signc(s33);
        rho = sqrt(tenorm*tenorm+2*vscalar);
        tas = tenorm*0.5*(epsilons33*rho*nb-tenorm*na);
        tbs = tenorm*0.5*(epsilons33*rho*na-tenorm*nb);
    else
        error('no solution');
    end
    
    t = [tas tbs];
    n = [na nb];
    R(:,:,1) = H*(eye(3)-2/vscalar*tas*na');
    R(:,:,2) = H*(eye(3)-2/vscalar*tbs*nb');
    
        
%         za1 = s12+sqrt(Ms33)/s22; zb1 = s12-sqrt(Ms33)/s22;
%         za3 = s23-epsilon13*sqrt(Ms11)/s22;zb3 = s23+epsilon13*sqrt(Ms11)/s22;
%         aa = 1 + za1*za1+za3*za3;
%         ab = 1 + zb1*zb1+zb3*zb3;
%         b = 2+trace(S);    
    
    varargout{1} = R;
    varargout{2} = t;
    varargout{3} = n;
end

function s = signc(a)
    if a > 0 || abs(a) < 1e-6
        s = 1;
    else
        s = -1;
    end
end