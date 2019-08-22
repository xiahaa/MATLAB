function dsicrete_trajectory_regression_on_manifold
    
    N1 = 4;
    Rdata = zeros(3,3*N1);
    N2 = 97;
    Rreg = zeros(3,3*N2);
    
    % fill in data to Rdata
    
    % initialize with piecewise geodesic path using park's method
    
    % start optimization
    iter = 1;
    maxiter = 1e6;
    
    oldcost = -1e6;
    newcost = 1e6;
    
    tol1 = 1e-3;
    
    while iter < maxiter
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        newcost = cost(xi,v,tau,lambda,miu);
        
        if abs(newcost - oldcost) < tol1
            break;
        end
        oldcost = newcost;
        
        % compute gradient, batch optimization
        % here, I will compute graduate for all so3s and then update as a
        % batch. I think perhaps sequencial update could be another option
        % since somehow this is a nonconvex optimization, no global minimum
        % is guaranteed.
        [LHS, RHS] = batch_sol(xi, v, indices, tau, lambda, miu, N2);
        
        dxis = -LHS\RHS;
        
        % update
        for j = 1:N2
            dxi = dxis(j*3-2:j*3);
            Rreg(:,j*3-2:j*3) = Rreg(:,j*3-2:j*3) * expSO3(dxi);
        end
        
        % TODO, do we need to check the norm of the gradient and exit if
        % the norm of gradient is lower than a threshold.
        
        iter = iter + 1;
    end
    
    
end

function [LHS, RHS] = batch_sol(xi, v, indices, tau, lambda, miu, N)
    lhs = zeros(3,3,N);
    rhs = zeros(3,N);
    
    % fist deal with term 1
    for i = 1:length(indices)
        Jr = rightJinv(xi(:,i));
        lhs(3,3,indices(i)) = lhs(3,3,indices(i)) + Jr'*Jr;
        rhs(3,indices(i)) = rhs(3,indices(i)) + Jr'*xi(:,i);
    end
    
    % second term
    % endpoints 
    c1 = lambda / tau;
    
    Jr = rightJinv(-v(:,1));
    lhs(3,3,1) = lhs(3,3,1) + Jr'*Jr.*c1;
    rhs(3,1) = rhs(3,1) + Jr'*(-v(:,1)).*c1;
    
    Jr = rightJinv(v(:,end));
    lhs(3,3,end) = lhs(3,3,end) + Jr'*Jr.*c1;
    rhs(3,end) = rhs(3,end) + Jr'*(v(:,end)).*c1;
    
    for i = 2:N-1
        Jr1 = rightJinv(v(:,i-1));
        Jr2 = rightJinv(-v(:,i));
        A1 = Jr1'*Jr1;
        b1 = Jr1'*v(:,i-1);
        A2 = Jr2'*Jr2;
        b2 = Jr2'*(-v(:,i));
        lhs(3,3,i) = lhs(3,3,i) + (A1+A2).*c1;
        rhs(3,i) = rhs(3,i) + (b1+b2).*c1;
    end
    
    % third term
    c2 = miu / tau^3;
    % end points
    Jr = rightJinv(-v(:,1));
    lhs(3,3,1) = lhs(3,3,1) + Jr'*Jr.*c2;
    rhs(3,1) = rhs(3,1) + Jr'*(-v(:,1)+v(:,2)).*c2;
    
    Jr = rightJinv(v(:,end));
    lhs(3,3,end) = lhs(3,3,end) + Jr'*Jr.*c2;
    rhs(3,end) = rhs(3,end) + Jr'*(-v(:,end-1)+v(:,end)).*c2;
    
    for i = 2:N-1
        Jr1 = rightJinv(-v(:,i-1));
        Jr2 = rightJinv(-v(:,i));
        Jr = Jr1+Jr2;
        Jrtb = Jr'*(-v(:,i) - v(:,i-1));
        lhs(3,3,i) = lhs(3,3,i) + Jr'*Jr.*c2;
        rhs(3,i) = rhs(3,i) + Jrtb.*c2;
    end
    
    LHS = spdiags(lhs);
    RHS = rhs(:);
end

function xi = data_term_error(Rdata,Rreg,indices)
    xi = zeros(3,length(indices));
    for i = length(indices)
        ii = indices(i);
        xi(:,i) = logSO3(Rdata(:,(i*3-2):i*3)'*Rreg(:,(ii*3-2):ii*3));
    end
end

function v = numerical_diff_v(Rreg)
    N = round(size(Rreg)/3);
    v = zeros(3,N-1);
    for i = 1:N-1
        v(:,i) = logSO3(Rreg(:,(i*3-2):i*3)'*Rreg(:,(i*3+1):i*3+3));
    end
end

function y = cost(xi,v,tau,lambda,miu)
    % cost term 1, data cost
    cost1 = sum(vecnorm(xi,2).^2.*2);
    
    % cost term 2, first order smooth cost, integrate with trapezoidal
    % rule, consistent with Boumal's paper. TODO change in paper.
    N = round(size(Rreg)/3);
    wv = [1 ones(1,N-2)];
    cost2 = sum(vecnorm(v,2).^2.*(2/tau).*wv);
    
    % cost term 3, second order smooth cost, integrate with trapezoidal
    % rule
    a = zeros(3,N-2);
    for i = 2:N-1
        a(:,i-1)=v(:,i)-v(:,i-1);
    end
    cost3 = sum(vecnorm(a,2).^2.*(2/tau^3));
    
    y = cost1 * 0.5 + cost2 * 0.5 * lambda + cost3 * 0.5 * miu;
end

function vechat = hat(vec)
    if size(vec,1) == 3 
        vechat = [  0,     -vec(3),  vec(2);
                vec(3),   0    , -vec(1);
               -vec(2),  vec(1),   0    ];  
    elseif size(vec,1) == 6
        vechat = [ hat( vec(4:6,1) ) vec(1:3,1); zeros(1,4) ];    
    end   
end

function R = expSO3(r)
    angle = norm(r);
    tol = 1e-12;
    if angle < tol
        R = eye(3);
        xM = eye(3);
        cmPhi = hat(phi);
        N = 10;% finite series
        for n = 1:N
            xM = xM * (cmPhi / n);
            R = R + xM;
        end
        [U,~,V] = svd(R);
        R = V*diag([1,1,det(V*U')])*U';% projection to SO3
    else
        so3 = hat(r);
        R = eye(3) + sin(angle)/angle*so3 + (1-cos(angle))/angle^2*so3^2;
    end
end

function r = logSO3(R)
    r = [0;0;0];
    angle = acos((trace(R)-1)*0.5);
    if angle > 1e-10
        r = angle / (2*sin(angle))*(R-R');
    else
        so3 = zeros(3,3);
        for i = 1:2
            so3 = so3 + (-1)^(i-1)/i.*(R-eye(3))^(i);
        end
        so3 = 0.5.*(so3-so3');
        r = [-so3(2,3);so3(1,3);-so3(1,2)];
    end
end


function J = leftJ(r)
    tolerance = 1e-12;
    angle = norm(r);
    if angle < tolerance
        % If the angle is small, fall back on the series representation
        N = 10;
        J = eye(3);
        pxn = eye(3);
        px = hat(r);
        for n = 1:N
            pxn = pxn*px/(n + 1);    
            J = J + pxn;
        end
    else
        axis = r/angle;

        cph = (1 - cos(angle))/angle;
        sph = sin(angle)/angle;

        J = sph * eye(3) + (1 - sph) * axis * axis' + cph * hat(axis);
    end       
end

function J = rightJ(r)
    tolerance = 1e-12;
    angle = norm(r);
    if angle < tolerance
        % If the angle is small, fall back on the series representation
        N = 10;
        J = eye(3);
        pxn = eye(3);
        px = -hat(r);
        for n = 1:N
            pxn = pxn*px/(n + 1);    
            J = J + pxn;
        end
    else
        axis = r/angle;

        cph = (1 - cos(angle))/angle;
        sph = sin(angle)/angle;

        J = sph * eye(3) + (1 - sph) * axis * axis' - cph * hat(axis);
    end    
end

function Jinv = leftJinv(r)
    tolerance = 1e-12;
    phi = r;
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation        
        Jinv = eye(3);
        pxn = eye(3);
        px = hat(vec);
        for n = 1:N
            pxn = pxn * px/n;
            Jinv = Jinv + lut_bernoullinumber(n) * pxn;
        end
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        Jinv =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               - ph_2 * hat(axis);
    end   
end

function Jinv = rightJinv(r)
    tolerance = 1e-12;
    phi = r;
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation        
        Jinv = eye(3);
        pxn = eye(3);
        px = -hat(vec);
        for n = 1:N
            pxn = pxn * px/n;
            Jinv = Jinv + lut_bernoullinumber(n) * pxn;
        end
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        Jinv =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               + ph_2 * hat(axis);
    end   
end

function b = lut_bernoullinumber(n)
    lut = [1,-0.500000000000000,0.166666666666667,0,-0.0333333333333333,0, ...
           0.0238095238095238,0,-0.0333333333333334,0,0.0757575757575749,0, ...
           -0.253113553113554,0,1.16666666666667,0,-7.09215686274515,0,54.9711779448621,0,-529.124242424242];
    b = lut(n+1);
end