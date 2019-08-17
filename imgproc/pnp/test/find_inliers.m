function varargout = find_inliers(v,n)
% trial to check if my idea works

    N = size(v,2);
    
    A = n';
 
    % to ll
    for i = 1:N
        ll(1,i) = find_theta(v(:,i));
        ll(2,i) = find_phi(v(:,i));
    end
    
    id = 1;
    
    % non heuristics
    min_cost1 = -1e6;
    min_cost2 = 1e6;
    best_err = [];
    min_d = [];
    
    threshold = 0.01;
    
%     rv = R*v;
%     rv = rv ./ sqrt(rv(1,:).^2 + rv(2,:).^2 + rv(3,:).^2 );
%     for i = 1:N
%         llr(1,i) = find_theta(rv(:,i));
%         llr(2,i) = find_phi(rv(:,i));
%     end
%     tmp1 = [cos(ll(2,:));sin(ll(2,:))];
%     tmp2 = [cos(llr(2,:));sin(llr(2,:))];
%     anglediff = acos(dot(tmp1,tmp2));
    
    phis = linspace(-pi,pi,360);
    for i = 1:length(phis)
        [ll1,valid] = proj_on_ring2(phis(i),n(:,id));
        [ll2,valid] = proj_on_ring2(phis(i),n(:,id+1));
        
        %
        vfake = [cos(ll1(2,:)).*cos(ll1(1,:));sin(ll1(2,:)).*cos(ll1(1,:));sin(ll1(1,:))];
        rotaxis = cross(v(:,id),vfake);rotaxis = rotaxis./norm(rotaxis);
        rotangle = acos(dot(vfake,v(:,id)));
        R1 = vec2rot(rotaxis.*rotangle);
        
        vfake = [cos(ll2(2,:)).*cos(ll2(1,:));sin(ll2(2,:)).*cos(ll2(1,:));sin(ll2(1,:))];
        rotaxis = cross(v(:,id+1),vfake);rotaxis = rotaxis./norm(rotaxis);
        rotangle = acos(dot(vfake,R1*v(:,id+1)));
        R2= vec2rot(rotaxis.*rotangle);
        
        R = R2*R1;
        
        corxyz1 = R*v;
        corxyz1 = corxyz1 ./ vecnorm(corxyz1,2);
        % dx dy
%         d1 = ll1 - ll(:,id);
%         d2 = ll2 - ll(:,id);
%         corll1 = ll + d1;
%         corll2 = ll + d2;
        % to xyz
%         corxyz1 = [cos(corll1(2,:)).*cos(corll1(1,:));sin(corll1(2,:)).*cos(corll1(1,:));sin(corll1(1,:))];
%         corxyz2 = [cos(corll2(2,:)).*cos(corll2(1,:));sin(corll2(2,:)).*cos(corll2(1,:));sin(corll2(1,:))];
        err1 = diag(abs(A * corxyz1));
%         err2 = diag(abs(A * corxyz2));
        indices = err1 < threshold;
        cost = sum(err1(indices));
        if sum(indices) > min_cost1 || (sum(indices)==min_cost1 && cost < min_cost2)
            min_cost1 = sum(indices);
            min_cost2 = cost;
            min_d = (err1 < threshold);
            best_err = err1;
        end
%         cost = sum(err2 < threshold);
%         if cost < min_cost
%             min_cost = cost;
%             min_d = d2;
%         end
        plot(err1);
        pause(0.1);
    end
    
%     corll = ll + min_d;
%     corxyz = [cos(corll(2,:)).*cos(corll(1,:));sin(corll(2,:)).*cos(corll(1,:));sin(corll(1,:))];
%     err = diag(abs(A * corxyz));
    
%     inlier_id = err < threshold;
    inlier_id = min_d;
    plot(best_err);
    
    
    % fail, cannot only use tilt angle to justify inliers and outliers
%     betas = zeros(1,N);
%     valids = zeros(1,N);
%     for i = 1:N
% %         tilt = find_tilt(n(:,i));
%         phi = find_phi(v(:,i));
%         
%         [p1,valid] = proj_on_ring(phi,n(:,i));
%         
%         if valid == true
%             valids(i) = 1;
%             beta1 = acos(dot(v(:,i),p1));
% %             beta2 = acos(dot(v(:,i),p2));
% %             beta = min(beta1,beta2);
%             betas(1,i) = beta1;
% %             betas(2,i) = beta2;
%         end
%     end
%     
%     % remove invalid
%     betas(:,valids == 0) = [];
%     betas = betas(:);
%     hist(betas,180);
%     varargout = [];
end

function theta = find_theta(v)
% compute relative latitude
    if norm(v(1:2)) > 1e-6
        theta = asin(v(3));
    else
        theta = sign(v(3))*pi/2;
    end
end

function phi = find_phi(v)
% compute relative longitude
    phi = atan2(v(2),v(1));
end

function tilt = find_tilt(n)
% compute tilt angle
    tilt = acos(dot([0;0;1],n));
end

function [p1,valid] = proj_on_ring(phi,n)
    if abs(n(3)) < 1e-6
        p1 = [];
%         p2 = [];
        valid = false;
    else
        c1 = -(cos(phi)*n(1)+sin(phi)*n(2));
        den = (c1)^2 + n(3)^2;
        num = n(3)^2;
        co1 = sqrt(num/den);
%         co2 = -sqrt(num/den);
        so1 = co1 * c1 / n(3);
%         so2 = co2 * c1 / n(3);
        
        p1 = [cos(phi)*co1;sin(phi)*co1;so1];
%         p2 = [cos(phi)*co2;sin(phi)*co2;so2];
        
        valid = true;
    end
end

function [ll1,valid] = proj_on_ring2(phi,n)
    if abs(n(3)) < 1e-6
        ll1 = [];
        valid = false;
    else
        c1 = -(cos(phi)*n(1)+sin(phi)*n(2));
        den = (c1)^2 + n(3)^2;
        num = n(3)^2;
        co1 = sqrt(num/den);
%         co2 = -sqrt(num/den);
        so1 = co1 * c1 / n(3);
        ll1 = [atan2(so1,co1);phi];
%         ll2 = [acos(co2);phi];
        valid = true;
    end
end