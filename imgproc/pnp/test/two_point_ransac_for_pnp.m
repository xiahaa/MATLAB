function [retinliers, cost, itercnt] = two_point_ransac_for_pnp(p,q)
% implementation of two point ransac for pnp
    minset = 2;
    maxiter = 1e6;
    iter = 0;

    bestcost = 0;
    inliers = [];
    % total number
    n = size(p,2);
    pd = 0.99;         % Desired probability of choosing at least one sample
                      % free from outliers (probably should be a parameter)
    idn = 1:n;
    idthreshold = round((n - 2)*0.25);
    
    while iter < maxiter
        
        % random pick two
        ss = [n-1,n];%randperm(n, minset);
        % for others form fourth order polynomial
        A = zeros(n-minset,5);
        k = 1;
        for i = 1:n
            if i == ss(1) || i == ss(2) continue; end
            ind = [ss i];
            A(k,:) = form_polynomial(p(:,ind),q(:,ind));
            k = k + 1;
        end
        
% 
%         v1= q(:,ss(1));
%         v2= q(:,ss(2));
%         cg1= v1.'*v2;
%         sg1= sqrt(1-cg1^2);
%         D1= norm(p(:,ss(1))-p(:,ss(2)));
%         k = 1;
%         for i= 1:n
%             if i == ss(1) || i == ss(2) continue; end
%             
%             vi= q(:,i);
%             cg2= v1.'*vi;
%             cg3= v2.'*vi;
%             sg2= sqrt(1-cg2^2);
%             D2= norm(p(:,ss(1))-p(:,i));
%             D3= norm(p(:,i)-p(:,ss(2)));
%             % get the coefficients of the P3P equation from each subset.
%             A(k,:)= getp3p(cg1,cg2,cg3,sg1,sg2,D1,D2,D3);
%             k = k + 1;
%         end

        
%         oldxi = inf;
        idraw = 1:n-2;
%         inlierloop = 1;
        idinlier = idraw;
        id2 = idn ~= ss(1) & idn ~= ss(2);r = zeros(n,0);
        
%         while inlierloop <= 50
            Ain = A(idinlier,:);
            % find the minimum singular value as well the rightmost singular
            % vector
            [~,~,V] = svd(Ain); 
            vmin = V(:,end);

            % residual
            rid2 = abs(A * vmin);  
            r(id2) = rid2;
            [r_sort] = sort(rid2,'ascend');
            inlier_threshold = max(r_sort(idthreshold), 1e-6);
            if inlier_threshold > oldxi
                break;
            else
                oldxi = inlier_threshold;
                idinlier = idraw(rid2<=inlier_threshold);
            end
%             inlierloop = inlierloop + 1;
%         end
        
        inlierid = r<oldxi;
        rinlier = r(inlierid);
        ninliers = numel(rinlier);
%         avg_err_per_inliers = sum(rinlier)/ninliers;
        if ninliers > bestcost
            bestcost = ninliers;
            inliers = inlierid;
%             ninliers = numel(rinlier); 
            % Update estimate of N, the number of trials to ensure we pick,
            % with probability p, a data set with no outliers.
            fracinliers =  ninliers/n;
            pNoOutliers = 1 -  fracinliers^minset;
            pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
            pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
            maxiter = log(1-pd)/log(pNoOutliers);
        end
        iter = iter + 1;
    end
    
    retinliers = inliers;
    cost = bestcost;
    itercnt = iter;
end

function B = getp3p(l1,l2,A5,C1,C2,D1,D2,D3)

A1= (D2/D1)^2;
A2= A1*C1^2-C2^2;
A3= l2*A5-l1;
A4= l1*A5-l2;
A6= (D3^2-D1^2-D2^2)/(2*D1^2);
A7= 1-l1^2-l2^2+l1*l2*A5+A6*C1^2;

B= [A6^2-A1*A5^2, 2*(A3*A6-A1*A4*A5), A3^2+2*A6*A7-A1*A4^2-A2*A5^2,...
    2*(A3*A7-A2*A4*A5), A7^2-A2*A4^2];

end