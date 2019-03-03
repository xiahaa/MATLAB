function SO3s = formeSO3(p)
    np = size(p,2);
    % mean
    center = mean(p,2);
    % diff
    dist = p - repmat(center,1,np);
    dist = diag(dist'*dist);
    % sort
    [val,sortid] = sort(dist,'ascend');
    % select from 25%-75%
    lowerbd = round(0.3*np);
    upperbd = round(0.6*np);
    
    if (upperbd-lowerbd) > 10
        lowerbd = round(0.5*np)- 5;
        upperbd = 10 + lowerbd - 1;
    end
    
    % filtering
    filterid = sortid(lowerbd:upperbd);
    % cn2
    setp = nchoosek(1:length(filterid),2);
    % pre allocate
    n_so3 = size(setp,1);
    SO3s = zeros(3,3,n_so3);
    k = 1;
    for i=1:n_so3
        ii = filterid(setp(i,1));
        jj = filterid(setp(i,2));
%         kk = setp(i,3);
        v1 = -center+p(:,ii); 
        v2 = -center+p(:,jj);
        
        if norm(v1) < norm(v2)
            tmp = v1;v1 = v2;v2 = tmp;
        end
         
        if norm(cross(v1,v2)) < 1e-3 || norm(v1) < 1e-3 || norm(v2) < 1e-3 
            continue;
        else
            %% case 1
            vn1 = v1./norm(v1);vn2 = v2./norm(v2);
            vn3 = cross(vn1,vn2);vn3 = vn3./norm(vn3);
            vn2 = cross(vn3,vn1);vn2 = vn2./norm(vn2);
            SO3s(:,:,k) = [vn1 vn2 vn3];
            k = k + 1;
            %% case 2
%             vn1 = v2./norm(v2);vn2 = v1./norm(v1);
%             vn3 = cross(vn1,vn2);vn3 = vn3./norm(vn3);
%             vn2 = cross(vn3,vn1);vn2 = vn2./norm(vn2);
%             SO3s(:,:,k) = [vn1 vn2 vn3];
%             k = k + 1;
            
        end
    end
end