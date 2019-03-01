function varargout = pose_estimation_correspondence_free(p, q, varargin) 
    %% clustering
    nq = size(q,2);
    np = size(p,2);
    
%     nqc = round(nq/10);
%     npc = round(np/10);
%     nn = max(nqc,npc);
%     
%     [~,cp] = kmeans(p', npc);
%     [~,cq] = kmeans(q', nqc);
%     p = cp';
%     q = cq';
%     
    %% step 1: forme SO(3) for qn and P
%     set_q = nchoosek(1:nq,2);
%     set_p = nchoosek(1:np,2);

    SO3_q = formeSO3(q);
    SO3_p = formeSO3(p);
    
%     [SO3_p, SO3_q] = formeSO32(p,q);

    %     SO3_q = rankingSO3(SO3_q);
    %     SO3_p = rankingSO3(SO3_p);
    
    wp = ones(1,size(SO3_p,3))./size(SO3_p,3);
    wq = ones(1,size(SO3_q,3))./size(SO3_q,3);
    
    maxiter = 1;
    iter = 1;
    while iter <= maxiter
        Mq1 = mean_1st_order(SO3_q, wq);
    %     Mq3 = mean_iterative_kron(SO3_q,Mq1);
        Mp1 = mean_1st_order(SO3_p, wp);
    %     Mp3 = mean_iterative_kron(SO3_p,Mp1);
    
        % update weight
        R = Mq1*Mp1';
        
        oldwp = wp;
        oldwq = wq;

%         [wp,wq] = update_weight(SO3_p,SO3_q,R,wp,wq);
        
        if iter > 1
            if norm(oldwp-wp) < 1e-6 && norm(oldwq-wq) < 1e-6
                break;
            end
        end
        
        
        iter = iter + 1;
    end
    
    % outlier
%     outlierp = abs(wp-mean(wp))>1*sqrt(var(wp));
%     outlierq = abs(wp-mean(wq))>1*sqrt(var(wq));
%     SO3_q = SO3_q(:,:,~outlierq);
%     wq = wq(~outlierq);
%     SO3_p = SO3_p(:,:,~outlierp);
%     wp = wp(~outlierp);
    % refine
    Mq1 = mean_1st_order(SO3_q,wq);
    Mp1 = mean_1st_order(SO3_p,wp);
    
    Mq3 = Mq1;%FNS_iterative2(SO3_q,Mq1,wq);
    Mp3 = Mp1;%FNS_iterative2(SO3_p,Mp1,wp);

    R1 = Mq3*Mp3';
%     R2 = Mq2*inv(Mp2)
    R2 = Mq1*Mp1';
    
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
    
function [wp,wq] = update_weight(Rp,Rq,R,wp,wq)
    m = size(Rp,3);
    n = size(Rq,3);
    
    so3q = zeros(3,n);
    for j = 1:n
        so3q(:,j) = rot2vec(Rq(:,:,j));
    end

    dists = zeros(m,n);
    for i = 1:m
        so3p = Rp(:,:,i);
        so3qp = rot2vec(R*so3p);
        err = so3q - repmat(so3qp,1,n);
        err = sqrt(diag(err'*err));
        dists(i,:) = err';
    end
    min1 = min(dists,[],2)';
    min2 = min(dists,[],1);
    gamma = 0.0001;
    min1(abs(min1) < 1e-6) = gamma;
    min2(abs(min2) < 1e-6) = gamma;
    
    wp = 1 ./ min1;
    wq = 1 ./ min2;
%     l1 = 1/m;
%     l2 = 1/n;
%     
%     beta = 0.001;
%     alpha = 3*pi/180.0;
    
%     wp = l1.*exp(-beta.*(min1-alpha))';
    wp = wp ./ norm(wp);
    
%     wq = l2.*exp(-beta.*(min2-alpha));
    wq = wq ./ norm(wq);
end

function goodSO3s = rankingSO3(SO3s)
% maybe it will be more robust to choose N salient orientation matrices
% rather than work with the whole set.
    angles = zeros(1,size(SO3s,3));
    for i = 1:size(SO3s,3)
        angles(i) = abs(norm(rot2vec(SO3s(:,:,i))));
    end
    [angless,id] = sort(angles);
    num = 100;
    if size(SO3s,3) < 100
        num = size(SO3s,3);
    end
    goodSO3s = SO3s(:,:,id(1:num));
end







