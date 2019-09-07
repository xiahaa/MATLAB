function Rreg = traj_smoothing_via_jc(Rreg, varargin)

    coeff = 1000;
    max_iter_num = 10;
    dataterm_indices = zeros(1,size(Rreg,2)/3);
    
    if nargin >= 2 
        indices = varargin{1}; 
        dataterm_indices(indices) = 1;
    else
        dataterm_indices = ones(1,size(Rreg,2)/3);
    end 
    if nargin >= 3 
        coeff = varargin{2}; 
    end
    if nargin >= 4 
        max_iter_num = varargin{3}; 
    end
    
    oldrotation = reshape(Rreg,3,3,[]);
    initrotation = oldrotation;
    
    % 2012-10-1
    % Chao Jia
    % smoothing function of a sequence of rotation matrices using manifold gradient descent, the smoothing term
    % is L2 norm manifold difference ||logm(A'*B)||
    % coeff is the weight of the L2 smoothing term

    frame_num = size(oldrotation,3);
    max_dis = 0.11; % the max Riemannian distance constraint; now it's used to crop 720*480 video into 540*360
    newrotation = initrotation;
    cost_value = [];

    % Armijo parameters
    beta = 0.5;
    sigma = 0.1;

    fprintf('starting smoothing using manifold optimization...\n');
    for i = 1:max_iter_num
        tic;
        copyrotation = newrotation;

        %%%%%%%%% compute the gradient %%%%%%%%%%%%%%%
        fprintf ('computing the gradient ... ');
        grad_f = zeros(3,3,frame_num);
        % first matrix
        x = copyrotation(:,:,1);
        p = oldrotation(:,:,1);
        xp1 = copyrotation(:,:,2);
        grad_f(:,:,1) = logm(p'*x).*dataterm_indices(1) + coeff.*logm(xp1'*x);
        % last matrix
        x = copyrotation(:,:,end);
        p = oldrotation(:,:,end);
        xm1 = copyrotation(:,:,end-1);
        grad_f(:,:,end) = logm(p'*x).*dataterm_indices(end) + coeff.*logm(xm1'*x);
        % others
        for j = 2:(frame_num-1)
            x = copyrotation(:,:,j);
            xm1 = copyrotation(:,:,j-1);
            xp1 = copyrotation(:,:,j+1);
            p = oldrotation(:,:,j);
            grad_f(:,:,j) = logm(p'*x).*dataterm_indices(j) + coeff.*logm(xm1'*x) + coeff.*logm(xp1'*x);
        end
        grad_f = real(grad_f);
        % refine it in case of accumulated error
        for j = 1:frame_num
            grad_f(:,:,j) = skew_refine(grad_f(:,:,j));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%% Compute the Hessian %%%%%%%%%%%%%%%%
        %%{
        fprintf('compute Hessian... ');
        gu = zeros(3*frame_num,1);
        for j = 1:frame_num
            gu(j*3-2:j*3) = tan_so3_dec(-grad_f(:,:,j));
        end
        Hess_f = zeros(3*frame_num, 3*frame_num);
        for j = 1:frame_num
            Hess_f(3*j-2:3*j, 3*j-2:3*j) = Hess_f(3*j-2:3*j, 3*j-2:3*j) + Hess_dist(oldrotation(:,:,j), copyrotation(:,:,j)).*dataterm_indices(j);
            if j < frame_num
                temp = coeff.*Hess_dist(copyrotation(:,:,j), copyrotation(:,:,j+1));
                Hess_f(3*j-2:3*j+3, 3*j-2:3*j+3) = Hess_f(3*j-2:3*j+3, 3*j-2:3*j+3) + [temp, -temp; -temp, temp];
            end
        end
        fprintf('compute 2-metric projection direction ... ');
        Hess_f = so3_hess_refine(Hess_f, grad_f, oldrotation, copyrotation, max_dis, dataterm_indices);
        du = Hess_f\gu;
        for j = 1:frame_num
            Hess_grad(:,:,j) = tan_so3_com(du(3*j-2:3*j));
        end
        Hess_grad = ro3_direction_refine(Hess_grad, grad_f, oldrotation, copyrotation, max_dis, dataterm_indices);
        %%}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%% Armijo line search %%%%%%%%%%%%%%%%
        fprintf('start Armijo ... ');
        m = -1;
        temp_fun_cost = inf;
        tangent_inner = 0;
        old_fun_cost = so3_smooth_fun(copyrotation, oldrotation, coeff, dataterm_indices);
        cost_value = [cost_value, old_fun_cost];
        while old_fun_cost - temp_fun_cost< -sigma*tangent_inner
            m = m+1;
            alpha = beta^m;
            d = alpha.*(Hess_grad);

            % Gradient projection for constrained optimization
            %%{
            temprotation = so3_group_update(copyrotation, d);
            d = so3_gradient_projection(oldrotation, temprotation, copyrotation, d, max_dis, dataterm_indices);
            %}

            temprotation = so3_group_update(copyrotation, d);	
            tangent_inner = so3_group_inner(grad_f,d);
            temp_fun_cost = so3_smooth_fun(temprotation, oldrotation, coeff, dataterm_indices);
            if m > 100
                break;
            end
        end
        fprintf('\n');
        if m > 100
            break;
        end
        newrotation = temprotation;
        toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    Rreg = reshape(newrotation,3,[]);
end






















