function SO3s = formeSO3(p, setp)
    n_so3 = size(setp,1);
    SO3s = zeros(3,3,n_so3);
    center = mean(p,2);
    k = 1;
    for i=1:n_so3
        ii = setp(i,1);
        jj = setp(i,2);
%         kk = setp(i,3);
        v1 = center-p(:,ii); 
        v2 = center-p(:,jj);
        
%         if norm(v1) < 1e-3 || norm(v2) < 1e-3 
%             continue;
%         else
            v1 = v1./norm(v1);v2 = v2./norm(v2);
            v3 = cross(v1,v2);v3 = v3./norm(v3);
            v2 = cross(v3,v1);v2 = v2./norm(v2);
            SO3s(:,:,k) = [v1 v2 v3];
            k = k + 1;
%         end
    end
end