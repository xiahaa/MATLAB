function Rinterp = interpolation_on_SO3_park(Rraw, traw, w0, a0, tinterp)
    if size(Rraw,3) == 1
        Rraw = reshape(Rraw,3,3,[]);
    end
    N0 = size(Rraw,3);
    N1 = length(tinterp);
    Rinterp = zeros(3,3,N1);
    % 1st, preprocessing
    ri = zeros(3,N0-1);
    Ar = zeros(3,3,N0-1);
    for i = 2:N0
        ri(:,i-1) = logSO3(Rraw(:,:,i-1)'*Rraw(:,:,i));
        Ar(:,:,i-1) = calcAr(ri(:,i-1));
    end
    
    c = zeros(3,N0-1);
    b = zeros(3,N0-1);
    a = zeros(3,N0-1);
    
    c(:,1) = w0;
    b(:,1) = a0/2;
    a(:,1) = ri(:,1) - b1 - c1;
    
    for i = 2:size(ri,2)
        s = ri(:,i);
        t = 3*a(:,i-1)+2*b(:,i-1)+c(:,i-1);
        u = 6*a(:,i-1)+2*b(:,i-1);
        c(:,i) = Ar(:,:,i-1)*c(:,i-1);
        ns = norm(s);
        b(:,i) = 0.5*(u-s'*t/ns^4*(2*cos(ns)+ns*sin(ns)-2)*(cross(s,t)) - (1-cos(ns))/ns^2*(cross(s,u)) ...
               + s'*t/ns^5*(3*sin(ns)-ns*cos(ns)-2*ns)*(cross(s,cross(s,t))) ...
               + (ns-sin(ns))/ns^3*(cross(t,cross(s,t))+cross(s,cross(s,u))));
        a(:,i) = s - b(:,i) - c(:,i);
    end
    
    % do interpolation
    for i = 1:N1
        id = find(traw>tinterp(i),1);
        if id == 1
            Rinterp(:,:,i) = Rraw(:,:,id);
        elseif isempty(id)
            Rinterp(:,:,i) = nan(3,3);
        else
            ti = traw(id-1);
            tj = traw(id);
            tau = (tinterp(i) - ti) / (tj-ti);
            Rinterp(:,:,i) = Rraw(:,:,id-1)*expSO3(a(:,id-1).*tau^3+b(:,id-1).*tau^2+c(:,id-1).*tau);% id-1 since a is 1 less than Rraw
        end
    end
end

function Ar = calcAr(r)
    theta = norm(r);
    rskew = skew(r);
    Ar = eye(3) - (1-cos(theta))/theta^2*rskew + (theta-sin(theta))/theta^3*rskew*rskew;
end