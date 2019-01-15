function varargout = Fundamental_est_8point(varargin)
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
        Fn = estimateFundmentalMatrix8Point(p1_normalized, p2_normalized);
        F = inv(T2)*Fn*T1;
    else
        [p1_normalized, T1] = point_normalization(p1);
        [p2_normalized, T2] = point_normalization(p2);
        %% RANSAV estimation of homography
        sampleNum = 8;
        maxRANSACCnt = 1e5;
        bestInlier = [];
        bestCost = 1e6;
        cnt = 1;
        N = size(p1_normalized,2);
        while(cnt < maxRANSACCnt)
           %% RANDOM select 4
           id = randperm(N,sampleNum);
           %% get sample feature points
           pts1 = p1_normalized(:,id);
           pts2 = p2_normalized(:,id);

           %% get fundmental matrix
           F1 = estimateFundmentalMatrix8Point(pts1,pts2);

           %% compute error
           reprojerror = FSampDist(F1,p1_normalized,p2_normalized);
            
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
        Fn = estimateFundmentalMatrix8Point(p1_normalized,p2_normalized);
        F = inv(T2)*Fn*T1;
    end
    
    %% constraint enforcement
    [U,S,V] = svd(F);
    F = U*diag([S(1,1),S(2,2),0])*V';
    
    varargout{1} = F;
end

function [ F ] = estimateFundmentalMatrix8Point( p1n, p2n )
    usePeusedoInverse = 0;
    if usePeusedoInverse == 1
        A = zeros(size(p1n,1),8);
        b = -ones(size(p1n,1),1);
        for i = 1:1:size(p1n,1)
            A(i,1) = p1n(i,1)*p2n(i,1);
            A(i,2) = p1n(i,1)*p2n(i,2);
            A(i,3) = p1n(i,1);
            A(i,4) = p1n(i,2)*p2n(i,1);
            A(i,5) = p1n(i,2)*p2n(i,2);
            A(i,6) = p1n(i,2);
            A(i,7) = p2n(i,1);
            A(i,8) = p2n(i,2);   
        %     A(i,9) = 1;
        end
        x = (A'*A)\(A'*b);
        F = [x(1) x(4) x(7);
         x(2) x(5) x(8);
         x(3) x(6) 1];
    else
        A = zeros(size(p1n,1),9);
        b = -ones(size(p1n,1),1);
        for i = 1:1:size(p1n,1)
            A(i,1) = p1n(i,1)*p2n(i,1);
            A(i,2) = p1n(i,1)*p2n(i,2);
            A(i,3) = p1n(i,1);
            A(i,4) = p1n(i,2)*p2n(i,1);
            A(i,5) = p1n(i,2)*p2n(i,2);
            A(i,6) = p1n(i,2);
            A(i,7) = p2n(i,1);
            A(i,8) = p2n(i,2);   
            A(i,9) = 1;
        end
        [~,~,V] = svd(A);
        x = V(:,end);
        F = [x(1) x(4) x(7);
         x(2) x(5) x(8);
         x(3) x(6) x(9)];
        F = F./F(3,3);
    end
end

function dist=FSampDist(F, p1, p2)
    den1 = F*p1;
    den2 = (p2'*F)';
    error = diag(p2'*den1);
    den1norm = (den1(1,:).^2+den1(2,:).^2);
    den2norm = (den2(1,:).^2+den2(2,:).^2);
    den = den1norm + den2norm;
    dist = error.^2 ./ den;
%     error = zeros(size(p1,1),1);
%     for i = 1:1:size(p1,2)
%         p2homo = [p2(:,i)];
%         p1homo = [p1(:,i)];    
%         err = p2homo'*F*p1homo;
%         den1 = F*p1homo;
%         den2 = p2homo'*F;
%         den = den1(1)^2+den1(2)^2+den2(1)^2+den2(2)^2;
%         error(i) = err / den;
%     end
%     dist = error;
end

