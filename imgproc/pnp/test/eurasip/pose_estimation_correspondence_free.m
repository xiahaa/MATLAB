function varargout = pose_estimation_correspondence_free(p, q, varargin) 
    nq = size(q,2);
    np = size(p,2);
    %% step 1: forme SO(3) for qn and P
    set_q = nchoosek(1:nq,2);
    set_p = nchoosek(1:np,2);

    SO3_q = formeSO3(q, set_q);
    SO3_p = formeSO3(p, set_p);

    %     SO3_q = rankingSO3(SO3_q);
    %     SO3_p = rankingSO3(SO3_p);
    Mq1 = mean_1st_order(SO3_q);
%         Mq2 = mean_iterative_kron(SO3_q,Mq1);
    Mq3 = FNS_iterative(SO3_q,Mq1);

    Mp1 = mean_1st_order(SO3_p);
%         Mp2 = mean_iterative_kron(SO3_p1,Mp1);
    Mp3 = FNS_iterative(SO3_p,Mp1);

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







