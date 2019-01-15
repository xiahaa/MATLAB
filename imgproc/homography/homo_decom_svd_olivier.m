function varargout = homo_decom_svd_olivier(varargin)
%% decompose homography matrix using olivier algo.
    H = varargin{1};
    format long;
    %% SVD 
    [U,S,V] = svd(H);
    %% 
    d1 = S(1,1);d2 = S(2,2);d3 = S(3,3);
    s = det(U)*det(V);%% a sign

    if abs(d1-d3) > 1e-6
        if abs(d1-d2) > 1e-6 && abs(d2-d3) > 1e-6
            c1 = sqrt((d1*d1-d2*d2)/(d1*d1-d3*d3));
            c2 = sqrt((d2*d2-d3*d3)/(d1*d1-d3*d3));
            xsol = zeros(3,4);
            xsol(:,1) = [c1;0;c2];
            xsol(:,2) = [c1;0;-c2];
            xsol(:,3) = [-c1;0;c2];
            xsol(:,4) = [-c1;0;-c2];
            %% case 1
            [Rd, td, nd] = case1(xsol, d1, d2, d3);
            dsols = zeros(1,size(td,2));
            dsols(1:4) = s * d2;
            dsols(5:8) = -s * d2;
        else
            if abs(d1-d2) < 1e-6
                %% case 2
                xsol = zeros(3,2);
                xsol(:,1) = [0;0;1];
                xsol(:,2) = [0;0;-1];
                [Rd, td, nd] = case2(xsol, d1, d2, d3, 1);
            else
                %% case 2
                xsol = zeros(3,2);
                xsol(:,1) = [1;0;0];
                xsol(:,2) = [-1;0;0];
                [Rd, td, nd] = case2(xsol, d1, d2, d3, 2);
            end
            dsols = zeros(1,size(td,2));
            dsols(1:2) = s * d2;
            dsols(3:4) = -s * d2;
        end
        [R,t,n] = final_sol(Rd, td, nd, U, V, s);
    else
        warning('Pure rotation');
        %% the 3rd case
        % case dd > 0
        td(:,1) = zeros(3,1);
        Rd(:,:,1) = eye(3);
%         R = U*Rd*V';
%         t = td;
        nd = zeros(3,1);%% undefined
        [R,t,n] = final_sol(Rd, td, nd, U, V, s);
        dsols = zeros(1,size(td,2));
        dsols(1) = s * d2;
        dsols(2) = -s * d2;
    end

    
    %% ambiguity removing
    m1 = varargin{2};
    num1 = H(3,1)*m1(1,1) + H(3,2)*m1(2,1) + H(3,3);
    valid1 = zeros(1,size(t,2));
    for i = 1:size(t,2)
        if (num1)/dsols(i) > 0
            if n(:,i)'*m1(:,1)/dsols(i) > 0
                valid1(i) = 1;
            end
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

function varargout = final_sol(Rd, td, nd, U, V, s)
    R = zeros(3,3,size(td,2));
    t = zeros(3,size(td,2));
    n = zeros(3,size(td,2));
    %% R,t,n
    for i = 1:size(td,2)
        R(:,:,i) = s*U*Rd(:,:,i)*V';
        t(:,i) = U*td(:,i);
        n(:,i) = V*nd(:,i);
    end
    varargout{1} = R;
    varargout{2} = t;
    varargout{3} = n;
end

function varargout = case1(xsol, d1, d2, d3)
    %% case dd > 0
    c1 = sqrt((d1*d1-d2*d2)*(d2*d2-d3*d3));
    c2 = (d1+d3)*d2;
    c3 = d2*d2+d1*d3;
    
    sin1 = c1/c2;sin2 = -c1/c2;
    cos1 = c3/c2;
    
    Rd = zeros(3,3,8);
    Rd(:,:,1) = [cos1 0 -sin1;0 1 0;sin1 0 cos1];
    Rd(:,:,2) = [cos1 0 -sin2;0 1 0;sin2 0 cos1];
    Rd(:,:,3) = [cos1 0 -sin2;0 1 0;sin2 0 cos1];
    Rd(:,:,4) = [cos1 0 -sin1;0 1 0;sin1 0 cos1];
    
    td = zeros(3,8);
    c4 = d1 - d3;
    td(:,1) = c4*[xsol(1,1);0;-xsol(3,1)];
    td(:,2) = c4*[xsol(1,2);0;-xsol(3,2)];
    td(:,3) = c4*[xsol(1,3);0;-xsol(3,3)];
    td(:,4) = c4*[xsol(1,4);0;-xsol(3,4)];
    
    %% case dd < 0
    c2 = (d1-d3)*d2;
    c3 = -d2*d2+d1*d3;
    
    sin1 = c1/c2;sin2 = -c1/c2;
    cos1 = c3/c2;
    
    Rd(:,:,5) = [cos1 0 sin1;0 -1 0;sin1 0 -cos1];
    Rd(:,:,6) = [cos1 0 sin2;0 -1 0;sin2 0 -cos1];
    Rd(:,:,7) = [cos1 0 sin2;0 -1 0;sin2 0 -cos1];
    Rd(:,:,8) = [cos1 0 sin1;0 -1 0;sin1 0 -cos1];
    
    c4 = d1 + d3;
    td(:,5) = c4*[xsol(1,1);0;xsol(3,1)];
    td(:,6) = c4*[xsol(1,2);0;xsol(3,2)];
    td(:,7) = c4*[xsol(1,3);0;xsol(3,3)];
    td(:,8) = c4*[xsol(1,4);0;xsol(3,4)];
    
    nd = zeros(3,8);
    nd(:,1:4) = xsol;
    nd(:,5:8) = xsol;
    
    varargout{1} = Rd;
    varargout{2} = td;
    varargout{3} = nd;
end

function varargout = case2(xsol, d1, d2, d3, type)
    %% case dd > 0
    Rd = zeros(3,3,4);
    Rd(:,:,1) = eye(3);
    Rd(:,:,2) = eye(3);
  
    td = zeros(3,4);
    c4 = d3 - d1;
    td(:,1) = c4*[xsol(1,1);0;xsol(3,1)];
    td(:,2) = c4*[xsol(1,2);0;xsol(3,2)];
    
    %% case dd < 0
    if type == 1
        Rd(:,:,3) = diag([-1,-1,1]);
        Rd(:,:,4) = diag([-1,-1,1]);
    else
        Rd(:,:,3) = diag([1,-1,-1]);
        Rd(:,:,4) = diag([1,-1,-1]);
    end
    
    c4 = d1 + d3;
    td(:,5) = c4*[xsol(1,1);0;xsol(3,1)];
    td(:,6) = c4*[xsol(1,2);0;xsol(3,2)];
    
    nd = zeros(3,4);
    nd(:,1:2) = xsol;
    nd(:,3:4) = xsol;
    
    varargout{1} = Rd;
    varargout{2} = td;
    varargout{3} = nd;
end