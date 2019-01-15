function varargout = homography_est(varargin)
%% homography estimation
    if nargin < 2
        error('at least two set of points!!');
    end
    
    if nargin >= 3
        refineWithRANSAC = varargin{3};
        if nargin == 4
            errorThresh = varargin{4};
            autoThresh = 0;
        else
            autoThresh = 1;
        end
    else
        refineWithRANSAC = 0;
    end
    
    %% assume homogeneous
    p1 = varargin{1};
    p2 = varargin{2};
    
    if refineWithRANSAC == 0
        [p1_normalized, T1] = point_normalization(p1);
        [p2_normalized, T2] = point_normalization(p2);
        Hn = Hestimation(p1_normalized, p2_normalized);
        H = inv(T2)*Hn*T1;
    else
        %% RANSAV estimation of homography
        sampleNum = 4;
        maxRANSACCnt = 1e5;
        bestInlier = [];
        bestCost = 1e6;
        cnt = 1;
        N = size(p1,2);
        while(cnt < maxRANSACCnt)
           %% RANDOM select 4
           id = randperm(N,sampleNum);
           %% get sample feature points
           pts1 = p1(:,id);
           pts2 = p2(:,id);

           %% get fundmental matrix
           H1 = Hestimation(pts1,pts2);

           %% compute error
           reprojerror = reprojectionError(H1,pts1,pts2);
            
           if autoThresh == 1
                mean_err = mean(reprojerror);
                std_err = std(reprojerror);
                inlier = (reprojerror-mean_err) < (1.959964 * std_err); %% 95% 1.959964, 90% 1.644854
           else
                %% get inlier
                inlier = reprojerror < errorThresh;
           end
           
           %% 
           inlier_error = reprojerror(inlier);
           avg_err = sum(inlier_error)/ length(inlier_error);

           %% update
           if (avg_err < bestCost)
               bestCost = avg_err;
               bestInlier = inlier;
               maxRANSACCnt = log(0.01) / log(1-(numel(inlier_error)/numel(reprojerror))^sampleNum);
           end
           cnt = cnt + 1; 
        end
        p1 = p1(:, bestInlier);
        p2 = p2(:, bestInlier);
        [p1_normalized, T1] = point_normalization(p1);
        [p2_normalized, T2] = point_normalization(p2);
        Hn = Hest(p1_normalized, p2_normalized);
        H = inv(T2)*Hn*T1;
    end
    varargout{1} = H;
end

function H = Hestimation(p1, p2)
    A = zeros(3*size(p1,2),9);

    for i = 1:1:size(p1,2)
        k = (i-1)*3+1;
        x1 = p1(1,i);
        x2 = p2(1,i);
        y1 = p1(2,i);
        y2 = p2(2,i);
        A(k,:)   = [0 -x1 x1*y2 0 -y1 y2*y1 0 -1 y2];
        A(k+1,:) = [x1 0 -x2*x1 y1 0 -y1*x2 1  0 -x2];
        A(k+2,:) = [-x1*y2 x2*x1 0 -y2*y1 y1*x2 0 -y2 x2 0];    
    end
    AA = A'*A;
    [~,~,V] = svd(AA);
    vmin = V(:,end);
    H = [vmin(1) vmin(4) vmin(7);...
         vmin(2) vmin(5) vmin(8);...
         vmin(3) vmin(6) vmin(9)];
end

function error = reprojectionError(H,p1,p2)
    error = zeros(1,size(p1,2));
    %% assume p1 and p2 are already in homogeous coordianates
    p1proj = H*p1;
    p2proj = inv(H)*p2;
    %% to inhomogeneous coordinate
    p1proj = p1proj ./ p1proj(3,:);
    p2proj = p2proj ./ p2proj(3,:);

    err1 = p1proj(1:2,:) - p2(1:2,:);   
    err2 = p2proj(1:2,:) - p1(1:2,:);

    err1_norm = sqrt(err1(1,:).^2+err1(2,:).^2);%
    err2_norm = sqrt(err2(1,:).^2+err2(2,:).^2);%

    error = err1_norm + err2_norm;
end