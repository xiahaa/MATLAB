function varargout = homo_decom_svd_mayi(varargin)
%% decompose homography matrix using olivier algo.
    H = varargin{1};
    format long;
    %% SVD 
    [U,S,V] = svd(H);
    
    d1 = S(1,1);d2 = S(2,2);d3 = S(3,3);
    v1 = V(:,1);v2 = V(:,2);v3 = V(:,3);
    
    if abs(d1-d3) > 1e-6
        c1 = sqrt(d1*d1-d3*d3);
        u1 = (v1*sqrt(1-d3*d3)+v3*(d1*d1-1))/c1;
        u2 = (v1*sqrt(1-d3*d3)-v3*(d1*d1-1))/c1;

        n1 = cross(v2,u1);n1=n1./norm(n1);
        
        n2 = cross(v2,u2);n2=n2./norm(n2);

        U1 = [v2 u1 n1];
        U2 = [v2 u2 n2];
        W1 = [H*v2 H*u1 cross(H*v2, H*u1)];
        W2 = [H*v2 H*u2 cross(H*v2, H*u2)];

        R1 = W1*U1';
        R2 = W2*U2';

        t1 = (H-R1)*n1;
        t2 = (H-R2)*n2;
    else
        error('Pure rotation');
    end
    
    R(:,:,1) = R1;R(:,:,2) = R1;R(:,:,3) = R2;R(:,:,4) = R2;
    t(:,1) = t1;t(:,2) = -t1;t(:,3) = t2;t(:,4) = -t2;
    n(:,1) = n1;n(:,2) = -n1;n(:,3) = n2;n(:,4) = -n2;
    
    %% ambiguity removing
    m1 = varargin{2};
    valid1 = zeros(1,size(t,2));
    for i = 1:size(t,2)
        if n(:,i)'*m1(:,1) > 0
            valid1(i) = 1;
        end
    end
    valid1 = valid1 == 1;
    Rf = R(:,:,valid1);
    tf = t(:,valid1);
    nf = n(:,valid1);
    
    
    varargout{1} = Rf;
    varargout{2} = tf;
    varargout{3} = nf;
end