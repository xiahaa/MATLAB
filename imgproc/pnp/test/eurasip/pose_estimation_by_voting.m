function varargout = pose_estimation_by_voting(p, q, varargin) 
    nq = size(q,2);
    np = size(p,2);
    
    %% use distance to select consistent subsets  
    % find consistent subsets
    subq = zeros(1,nq);
    subp = zeros(1,np);
    threshold = 0.15;
    tmp = 1:np;
    for i = 1:nq
        
        
        
        after_threshold = tmp(cross_error(i,:) < threshold);
        if ~isempty(after_threshold)
            subq(i) = 1;
            subp(after_threshold) = 1;
        end
    end
    
    % only choose consistent data
    q = q(:,subq == 1);
    p = p(:,subp == 1);
    
    nq = size(q,2);
    np = size(p,2);
    
    % form SO3
    SO3_q = formeSO3self(q);
    SO3_p = formeSO3self(p);

    nq1 = size(SO3_q,3);
    np1 = size(SO3_p,3);
    
%     Mq1 = mean_1st_order(SO3_q);
%     Mp1 = mean_1st_order(SO3_p);
%     
%     Mq3 = FNS_iterative(SO3_q,Mq1);
%     Mp3 = FNS_iterative(SO3_p,Mp1);
% 
%     R1 = Mq3*Mp3';
% %     R2 = Mq2*inv(Mp2)
%     R2 = Mq1*Mp1';
    
%     
    % voting, brute-force
    SO3_tb = zeros(3,3,nq1*np1);
    so3_tb = zeros(3,nq1*np1);
    vote = zeros(1,nq1*np1);
    pairs = cell(nq1*np1);
    k = 1;
    for i = 1:nq1
        R1 = SO3_q(:,:,i);
        for j = 1:np1
            R2 = SO3_p(:,:,j);
            R = R1*R2';
            so3r = rot2vec(R);
            if k == 1
                so3_tb(:,k) = so3r;
                SO3_tb(:,:,k) = R;
                vote(k) = 1;
%                 pairs(k) = [pairs(k);[i,j]];
                k = k + 1;
            else
                err = so3_tb(:,1:k-1) - repmat(so3r,1,k-1);
                err = sqrt(diag(err'*err));
                [minval,minid] = min(err);
                if minval < 0.2
                    vote(minid) = vote(minid) + 1;
%                     pairs(minid) = [pairs(minid);[i,j]];
                else
                    so3_tb(:,k) = so3r;
                    SO3_tb(:,:,k) = R;
                    vote(k) = 1;
%                     pairs(k) = [pairs(k);[i,j]];
                    k = k + 1;
                end
            end
        end
    end
    [~,maxid] = max(vote(1:k-1));
    Rvote = SO3_tb(:,:,maxid);
    R1 = Rvote;
    
    pbar = mean(p,2);
    qbar = mean(q,2);
    t1 = qbar - R1*pbar;
    t2 = qbar - R2*pbar;
    
    varargout{1} = R1;
    varargout{2} = t1;
    
    if nargout > 2
        varargout{3} = R2;
        varargout{4} = t2;
    end
end
    
function SO3s = formeSO3self(p)
    np = size(p,2);
    % mean
    center = mean(p,2);
    % diff
    dist = p - repmat(center,1,np);
    dist = diag(dist'*dist);
    % sort
    [~,sortid] = sort(dist,'ascend');
    % select from 25%-75%
    lowerbd = round(0.4*np);
    upperbd = round(0.8*np);
    
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






