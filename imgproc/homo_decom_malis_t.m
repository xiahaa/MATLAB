function varargout = homo_decom_malis_t(varargin)
%% decompose homography matrix using olivier algo.
    H = varargin{1};
    format long;
    
    hr1 = H(1,:)'; hr2 = H(2,:)'; hr3 = H(3,:)';
    
    Sr = H*H'-eye(3);    
    sr11 = Sr(1,1);sr12 = Sr(1,2);sr13 = Sr(1,3);sr22 = Sr(2,2);sr23 = Sr(2,3);sr33 = Sr(3,3);
    
    Msr11 = sr23*sr23 - sr22*sr33;
    Msr22 = sr13*sr13 - sr11*sr33;
    Msr33 = sr12*sr12 - sr11*sr22;
    
    Msr12 = sr23*sr13 - sr12*sr33;
    Msr13 = sr22*sr13 - sr12*sr23;
    Msr23 = sr12*sr13 - sr11*sr23;
    
    vscalar = 2*sqrt(1+trace(Sr)-Msr11-Msr22-Msr33);
    tenorm = sqrt(2+ trace(Sr)-vscalar);
    rho = sqrt(tenorm*tenorm+2*vscalar);

    if abs(sr22) > 1e-6
%         Ms13 = -hr1'*(eye(3)+skewm(hr2)*skewm(hr2))*hr3;
        epsilonr13 = signc(Msr13);
        
        c1 = hr1'*hr2;
        c2 = norm(hr1);
        c3 = norm(hr2);
        c4 = norm(hr3);
        c5 = hr2'*hr3;
        
        tea = [c1+sqrt(c1*c1-(c2*c2-1)*(c3*c3-1)); ...
               c3*c3-1; ...
               c5-epsilonr13*sqrt(c5*c5)-(c3*c3-1)*(c4*c4-1)];

        teb = [c1-sqrt(c1*c1-(c2*c2-1)*(c3*c3-1)); ...
               c3*c3-1; ...
               c5+epsilonr13*sqrt(c5*c5)-(c3*c3-1)*(c4*c4-1)];
        
        ta = tenorm*tea./norm(tea);
        tb = tenorm*teb./norm(teb);
        
        epsilonsr22 = signc(sr22);
        na = 0.5*(epsilonsr22 * rho / tenorm * tb - ta);
        nb = 0.5*(epsilonsr22 * rho / tenorm * ta - tb);
    elseif abs(sr11) > 1e-6
        epsilonr23 = signc(Msr23);
        
        c1 = hr1'*hr2;
        c2 = norm(hr1);
        c3 = norm(hr2);
        c4 = norm(hr3);
        c5 = hr1'*hr3;
        
        tea = [c2*c2-1; ...
               c1+sqrt(c1*c1-(c2*c2-1)*(c3*c3-1)); ...
               c5+epsilonr23*sqrt(c5*c5)-(c2*c2-1)*(c4*c4-1)];

        teb = [c2*c2-1; ...
               c1-sqrt(c1*c1-(c2*c2-1)*(c3*c3-1)); ...
               c5-epsilonr23*sqrt(c5*c5)-(c2*c2-1)*(c4*c4-1)];
        
        ta = tenorm*tea./norm(tea);
        tb = tenorm*teb./norm(teb);
        
        epsilonsr11 = signc(sr11);
        na = 0.5*(epsilonsr11 * rho / tenorm * tb - ta);
        nb = 0.5*(epsilonsr11 * rho / tenorm * ta - tb);
        
    elseif abs(sr33) > 1e-6
        epsilonr12 = signc(Msr12);
        
        c1 = hr1'*hr3;
        c2 = norm(hr1);
        c3 = norm(hr2);
        c4 = norm(hr3);
        c5 = hr2'*hr3;
        
        tea = [c1+epsilonr12*sqrt(c1*c1-(c2*c2-1)*(c4*c4-1)); ...
               c5+sqrt(c5*c5)-(c3*c3-1)*(c4*c4-1); ...
               c4*c4-1];
               
        teb = [c1-epsilonr12*sqrt(c1*c1-(c2*c2-1)*(c4*c4-1)); ...
               c5-sqrt(c5*c5)-(c3*c3-1)*(c4*c4-1); ...
               c4*c4-1];
           
        ta = tenorm*tea./norm(tea);
        tb = tenorm*teb./norm(teb);
           
        epsilonsr33 = signc(sr33);
        na = 0.5*(epsilonsr33 * rho / tenorm * tb - ta);
        nb = 0.5*(epsilonsr33 * rho / tenorm * ta - tb);
    else
        error('no definition!!!');
    end
    
    t = [ta tb];
    n = [na nb];
    R(:,:,1) = (eye(3)-2/vscalar*ta*na')*H;
    R(:,:,2) = (eye(3)-2/vscalar*tb*nb')*H;
    
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