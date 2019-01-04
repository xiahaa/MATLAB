function varargout = point_normalization(varargin)
%% apply normalization for input point cloud (assume homogeneous, 3xn)
%
    p = varargin{1};
    Mean1 = mean(p, 2);
    p(1,:) = p(1,:) - Mean1(1);
    p(2,:) = p(2,:) - Mean1(2);
    % p(3,:) = p(3,:) - Mean1(1);
    p12 = p(1:2,:);
    S1 = sqrt(2)/mean(sqrt(diag(p12'*p12)));
    p(1:2,:) = p(1:2,:)*S1;
    T1=[S1 0 -Mean1(1)*S1; ... 
        0 S1 -Mean1(2)*S1; ...
        0 0 1];
    p_normalized = p;
    varargout{1} = p_normalized;
    varargout{2} = T1;
end