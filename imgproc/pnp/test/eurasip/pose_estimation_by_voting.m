function varargout = pose_estimation_by_voting(p, q, varargin) 
    %% clustering
    nq = size(q,2);
    np = size(p,2);
    
    % form SO3
    SO3_q = formeSO3(q);
    SO3_p = formeSO3(p);

    nq1 = size(SO3_q,3);
    np1 = size(SO3_p,3);
    
%     % form so3
%     so3_q = zeros(3,nq1);
%     so3_p = zeros(3,np1);
%     for i = 1:nq1
%         so3_q(:,i) = rot2vec(SO3_q(:,:,i));
%     end
%     for i = 1:np1
%         so3_p(:,i) = rot2vec(SO3_p(:,:,i));
%     end
    
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
                if minval < 0.1
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
    

%     R1 = Mq3*Mp3';
%     R2 = Mq2*inv(Mp2)
%     R2 = Mq1*Mp1';
    
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
    







