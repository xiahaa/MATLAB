% this implements the shooting method for optimal trajectory generation on
% SO3.
function shooting_method()
%     [wdd0,fval] = fminsearch(@solver,[0,0,0]');
%     disp(wdd0);
%     disp(fval);
%     
    r0 = [0.2,0.1,0.1]';
    r1 = [0.6,0.4,0.4]';
%     
    w0 = [1.5,0.1,0.1]';
    wd0 = [0.5,0.1,0.1]';
%     [t,x]=ode45(@fun,[0 1],[vec(R0);w0;wd0;wdd0],[]);
%     s = length(t);
%     reshape(x(s,1:9),3,3)
%     R1

    [R1opt,w1opt,~,costopt,tspan] = traj_gen_by_shooting(r0,r1,w0,wd0);
    
%     costopt
    
    R0 = expSO3(r0);
    R1 = expSO3(r1);
    cost_cubic = traj_gen_by_cubic(R0,R1,w0,R1*w1opt',tspan);
    
%     cost_cubic
    
    [cost_quintic,R1quintic] = traj_gen_by_quintic(R0,R1,w0,R1*w1opt',tspan);
    
%     cost_quintic
    
    disp(['costopt: ',num2str(costopt),', cost_cubic: ',num2str(cost_cubic),', cost_quintic: ',num2str(cost_quintic)]);
    disp(['shooting: ',num2str(norm(logSO3(R1'*R1opt)))]);
    disp(['optimization: ',num2str(norm(logSO3(R1'*R1quintic)))]);
end

function varargout = traj_gen_by_quintic(R0,R1,w0,w1,tspan)
    dr1 = logSO3(R0'*R1);
    % polynomial via optimization, lift order from cubic to quintic
%     1 0 0 1 0 0 1 0 0 1 0 0 1 0 0; ...
%        0 1 0 0 1 0 0 1 0 0 1 0 0 1 0; ...
%        0 0 1 0 0 1 0 0 1 0 0 1 0 0 1; ...
    
    Aeq = [ ...
           0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; ...
           5 0 0 4 0 0 3 0 0 2 0 0 1 0 0; ...
           0 5 0 0 4 0 0 3 0 0 2 0 0 1 0; ...
           0 0 5 0 0 4 0 0 3 0 0 2 0 0 1];
   beq = [w0;calcAr(dr1)\w1];
   Q1 = [[ 25/9,    0,    0,  5/2,    0,    0, 15/7,    0,    0, 5/3,   0,   0, 1, 0, 0]
        [    0, 25/9,    0,    0,  5/2,    0,    0, 15/7,    0,   0, 5/3,   0, 0, 1, 0]
        [    0,    0, 25/9,    0,    0,  5/2,    0,    0, 15/7,   0,   0, 5/3, 0, 0, 1]
        [  5/2,    0,    0, 16/7,    0,    0,    2,    0,    0, 8/5,   0,   0, 1, 0, 0]
        [    0,  5/2,    0,    0, 16/7,    0,    0,    2,    0,   0, 8/5,   0, 0, 1, 0]
        [    0,    0,  5/2,    0,    0, 16/7,    0,    0,    2,   0,   0, 8/5, 0, 0, 1]
        [ 15/7,    0,    0,    2,    0,    0,  9/5,    0,    0, 3/2,   0,   0, 1, 0, 0]
        [    0, 15/7,    0,    0,    2,    0,    0,  9/5,    0,   0, 3/2,   0, 0, 1, 0]
        [    0,    0, 15/7,    0,    0,    2,    0,    0,  9/5,   0,   0, 3/2, 0, 0, 1]
        [  5/3,    0,    0,  8/5,    0,    0,  3/2,    0,    0, 4/3,   0,   0, 1, 0, 0]
        [    0,  5/3,    0,    0,  8/5,    0,    0,  3/2,    0,   0, 4/3,   0, 0, 1, 0]
        [    0,    0,  5/3,    0,    0,  8/5,    0,    0,  3/2,   0,   0, 4/3, 0, 0, 1]
        [    1,    0,    0,    1,    0,    0,    1,    0,    0,   1,   0,   0, 1, 0, 0]
        [    0,    1,    0,    0,    1,    0,    0,    1,    0,   0,   1,   0, 0, 1, 0]
        [    0,    0,    1,    0,    0,    1,    0,    0,    1,   0,   0,   1, 0, 0, 1]];
   Q2 = [[ 400/7,     0,     0,    40,     0,     0, 24,  0,  0, 10,  0,  0, 0, 0, 0]
        [     0, 400/7,     0,     0,    40,     0,  0, 24,  0,  0, 10,  0, 0, 0, 0]
        [     0,     0, 400/7,     0,     0,    40,  0,  0, 24,  0,  0, 10, 0, 0, 0]
        [    40,     0,     0, 144/5,     0,     0, 18,  0,  0,  8,  0,  0, 0, 0, 0]
        [     0,    40,     0,     0, 144/5,     0,  0, 18,  0,  0,  8,  0, 0, 0, 0]
        [     0,     0,    40,     0,     0, 144/5,  0,  0, 18,  0,  0,  8, 0, 0, 0]
        [    24,     0,     0,    18,     0,     0, 12,  0,  0,  6,  0,  0, 0, 0, 0]
        [     0,    24,     0,     0,    18,     0,  0, 12,  0,  0,  6,  0, 0, 0, 0]
        [     0,     0,    24,     0,     0,    18,  0,  0, 12,  0,  0,  6, 0, 0, 0]
        [    10,     0,     0,     8,     0,     0,  6,  0,  0,  4,  0,  0, 0, 0, 0]
        [     0,    10,     0,     0,     8,     0,  0,  6,  0,  0,  4,  0, 0, 0, 0]
        [     0,     0,    10,     0,     0,     8,  0,  0,  6,  0,  0,  4, 0, 0, 0]
        [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0]
        [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0]
        [     0,     0,     0,     0,     0,     0,  0,  0,  0,  0,  0,  0, 0, 0, 0]];
    
    AA1 = [1 0 0 1 0 0 1 0 0 1 0 0 1 0 0; ...
           0 1 0 0 1 0 0 1 0 0 1 0 0 1 0; ...
           0 0 1 0 0 1 0 0 1 0 0 1 0 0 1];
    AA2 = 0.5.*hat(dr1)*AA1;
    Aineq = AA1-AA2;
    Aineq = [Aineq;-Aineq];
    eps = 0.1;
    bineq = [eps+dr1;eps-dr1];
    
    x = quadprog(Q2,[],Aineq,bineq,Aeq,beq);
    a = x(1:3);b = x(4:6);c=x(7:9);d=x(10:12);e=x(13:15);
    
     % numerical integration
%      ws = zeros(length(tspan),3);
%      for i = 1:length(tspan)
%         t = tspan(i);
%         ri = a.*t^5+b.*t^4+c.*t^3+d.*t^2+e.*t;
%         Ri = R0*expSO3(ri);
%         wi = calcAr(ri)*(a.*(5*t^4)+b.*(4*t^3)+c.*(3*t^2)+d.*(2*t)+e);
%         ws(i,:) = Ri'*wi;
%      end
%      dt = diff(tspan);
%      wds = diff(ws,1,1)./dt;
%      cost1 = sum(vecnorm(wds,2,2).^2.*dt);
%      
     % analytical sol
     cost2 = costFunc2(a,b,c,d,e);
     
     varargout{1} = cost2;
     varargout{2} = R0*expSO3(a+b+c+d+e);
end

function varargout = traj_gen_by_cubic(R0,R1,w0,w1,tspan)
    dr1 = logSO3(R0'*R1);
    % polynomial, utilize boundary values to derive a closed-form solution.
    A = [1 0 0 1 0 0 1 0 0; ...
         0 1 0 0 1 0 0 1 0; ...
         0 0 1 0 0 1 0 0 1; ...
         0 0 0 0 0 0 1 0 0; ...
         0 0 0 0 0 0 0 1 0; ...
         0 0 0 0 0 0 0 0 1; ...
         3 0 0 2 0 0 1 0 0; ...
         0 3 0 0 2 0 0 1 0; ...
         0 0 3 0 0 2 0 0 1];
     p = A\[dr1;w0;calcAr(dr1)\w1];
     %
     a = p(1:3);b = p(4:6);c = p(7:9);
          
     % numerical integration
%      ws = zeros(length(tspan),3);
%      for i = 1:length(tspan)
%         t = tspan(i);
%         ri = a.*t^3+b.*t^2+c.*t;
%         Ri = R0*expSO3(ri);
%         wi = calcAr(ri)*(a.*(3*t^2)+b.*(2*t)+c);
%         ws(i,:) = Ri'*wi;
%      end
%      dt = diff(tspan);
%      wds = diff(ws,1,1)./dt;
%      cost1 = sum(vecnorm(wds,2,2).^2.*dt);
     
     % analytical sol
     cost2 = costFunc1(a,b,c);
     
     varargout{1} = cost2;
end

function cost = costFunc1(a,b,c)
    cost = integral(@intf,0,1);
    function y = intf(t)
        rt = a.*t.^3 + b.*t.^2 + c.*t;
        drt = 3.*a.*t.^2 + 2.*b.*t + c;
        ddrt = 6.*a.*t + 2.*b;

        normr = vecnorm(rt);

        alphat = ddrt - dot(rt,drt)./normr.^4.*(2*cos(normr)+normr.*sin(normr)-2).*(cross(rt,drt)) - (1-cos(normr))./normr.^2.*(cross(rt,ddrt)) ...
                 + dot(rt,drt)./normr.^5.*(3*sin(normr)-normr.*cos(normr)-2*normr).*(cross(rt, cross(rt,drt))) ...
                 + (normr-sin(normr))./normr.^3.*(cross(drt, cross(rt,drt))+cross(rt, cross(rt,ddrt)));
             
        y = dot(alphat,alphat);
    end
end

function cost = costFunc2(a,b,c,d,e)
    cost = integral(@intf,0,1);
    function y = intf(t)
        rt = a.*t.^5 + b.*t.^4 + c.*t.^3+d.*t.^2+e.*t;
        drt = 5*a.*t.^4 + 4*b.*t.^3 + 3*c.*t.^2+2*d.*t+e;
        ddrt = 20*a.*t.^3 + 12*b.*t.^2 + 6*c.*t+2*d;

        normr = vecnorm(rt);

        alphat = ddrt - dot(rt,drt)./normr.^4.*(2*cos(normr)+normr.*sin(normr)-2).*(cross(rt,drt)) - (1-cos(normr))./normr.^2.*(cross(rt,ddrt)) ...
                 + dot(rt,drt)./normr.^5.*(3*sin(normr)-normr.*cos(normr)-2*normr).*(cross(rt, cross(rt,drt))) ...
                 + (normr-sin(normr))./normr.^3.*(cross(drt, cross(rt,drt))+cross(rt, cross(rt,ddrt)));
             
        y = dot(alphat,alphat);
    end
end

function Ar = calcAr(r)
    theta = norm(r);
    rskew = skew(r);
    if theta < 1e-6
        Ar = eye(3);
    else
        Ar = eye(3) - (1-cos(theta))/theta^2*rskew + (theta-sin(theta))/theta^3*rskew*rskew;
    end
end

function varargout = traj_gen_by_shooting(r0,r1,w0,wd0)
    R0 = expSO3(r0);
    R1 = expSO3(r1);
    
    function F=solver_shooting(wdd0)
        tspan = [0 1];
        [t,x]=ode45(@fun,tspan,[vec(R0);w0;wd0;wdd0],[]);
        s = length(t);
        F = norm(reshape(x(s,1:9),3,3)-R1,'fro');
    end
    % find wdd0 by shooting
    [wdd0,fval] = fminsearch(@solver_shooting,[0,0,0]');
    % re-integrate using ode45
    tspan = 0:0.001:1;
    [t,x]=ode45(@fun,tspan,[vec(R0);w0;wd0;wdd0],[]);
    % some results
    s = length(t);
    R1opt = reshape(x(s,1:9),3,3);
    disp(norm(logSO3(R1'*R1opt)));
    wds = x(2:end,13:15);
    dt = diff(t);
    cost = sum(vecnorm(wds,2,2).^2.*dt);
    varargout{1} = R1opt;% R1
    varargout{2} = x(s,10:12);% w1
    varargout{3} = x(s,13:15);% wd1
    varargout{4} = cost;
    varargout{5} = t;
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
        cmPhi = hat(r);
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
        so3 = angle / (2*sin(angle))*(R-R');
        r = [-so3(2,3);so3(1,3);-so3(1,2)];
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