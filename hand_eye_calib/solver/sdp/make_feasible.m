function xf = make_feasible(x, varargin)
    n = varargin{1};
    scales = zeros(n,1);
    for i = 1:n
        A = varargin{i+1};
        scales(i) = 1./sqrt(x'*A*x);
    end
    [~,id] = max(scales);
    xf = x .* scales(id);
end