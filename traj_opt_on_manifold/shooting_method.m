% this implements the shooting method for optimal trajectory generation on
% SO3.
function shooting_method()
%     cost = 1e6;
%     minuvk = [0;0;0];
%     
%     for u = -10:0.1:10
%         for v = -10:0.1:10
%             for k = -10:0.1:10
%                 Rf = numrical_integration(R0,w0,a0,[u,v,k]');
%                 if norm(logSO3(R1'*Rf)) < cost
%                     cost = norm(logSO3(R1'*Rf));
%                     minuvk = [u,v,k]';
%                 end
%             end
%         end
%     end
%     Rf = numrical_integration(R0,w0,a0,minuvk)
%     R1
%     minuvk

%     tspan = 0:0.1:10;
%     x0 = 1;
%     each = 2;
%     [t,x]=ode45(@fun,tspan,[0,x0],[],each);
%     kk = 1;
%     for u = -3:0.1:3
%         for v = -3:0.1:3
%             for k = -3:0.1:3
%                 cost(kk)=solver(R0,w0,a0,[u,v,k]',R1);
%                 kk = kk+1;
%                 disp(kk)
%             end
%         end
%     end

    [wdd0,fval] = fminsearch(@solver,[0,0,0]');
    disp(wdd0);
    disp(fval);
    
    r0 = [0.2,0.1,0.1]';
    r1 = [0.6,0.4,0.4]';
    R0 = expSO3(r0);
    R1 = expSO3(r1);
    w0 = [5,0.1,0.1]';
    wd0 = [0.5,0.1,0.1]';
    [t,x]=ode45(@fun,[0 1],[vec(R0);w0;wd0;wdd0],[]);
    s = length(t);
    reshape(x(s,1:9),3,3)
    R1
end

%%% test
function F=solver(wdd0)
    % example
    r0 = [0.2,0.1,0.1]';
    r1 = [0.6,0.4,0.4]';
    R0 = expSO3(r0);
    R1 = expSO3(r1);
    w0 = [5,0.1,0.1]';
    wd0 = [0.5,0.1,0.1]';

    tspan = [0 1];
    [t,x]=ode45(@fun,tspan,[vec(R0);w0;wd0;wdd0],[]);
    s = length(t);
    F = norm(reshape(x(s,1:9),3,3)-R1,'fro');
    
%     figure(1);
%     plot(cost);
end

%%%
function dy = fun(t,y)
    dy = zeros(3*3+9,1);
    dyb = zeros(9,1);
    
    R = reshape(y(1:9),3,3);
    % y in global sense
    yb = R*reshape(y(10:end),3,3);
    yb = yb(:);
    
    dy(1:9) = vec(R*hat(yb(1:3)));
    dyb(1:3) = yb(4:6);%dw
    dyb(4:6) = yb(7:9);%ddw
    dyb(7:9) = -cross(yb(1:3),yb(7:9));
    
    dy(10:end) = vec(R'*reshape(dyb,3,3));
end


function R1 = numrical_integration(R0,w0,wd0,wdd0)
    t = 0;
    dt = 0.02;
    w = w0;
    wd = wd0;
    wdd = wdd0;
    R = R0;
    while t <= 1
        % solve Euler-Lagrange
        wddd = -cross(w,wdd);
        
        Rp = R;
        R = R*expm(R*hat(w).*dt);
        w = w + wd*dt;
        wd = wd + wdd*dt;
        wdd = wdd + wddd * dt;
        
        % project to new frame
        w = R'*Rp*w;
        wd = R'*Rp*wd;
        wdd = R'*Rp*wdd;
        
        t = t + dt;
    end
    R1 = R;
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