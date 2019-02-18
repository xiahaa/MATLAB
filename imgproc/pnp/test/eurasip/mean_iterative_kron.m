function M = mean_iterative_kron(SO3s,M_1st)
    M = M_1st;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % use sum (v'*v) as the cost function, v = logm(R'Rt);
    % if you do like this, then the mean is not unique.
    %
    % iter = 1;
    % maxIter = 50;
    % err = 0;
    % old_err = 1e6;
    %     
    %     while iter < maxIter
    %         e = zeros(3,1);
    %         err = 0;
    %         for i = 1:size(SO3s,3)
    %             Rbar = M'*SO3s(:,:,i);
    %             phibar = rot2vec(Rbar);
    %             err = err + norm(phibar);
    %             e = e + M*phibar;
    %         end
    %         disp([old_err, err])
    %         if abs(err - old_err) < 1e-6 
    %             break;
    %         end
    %         xhat = e./size(SO3s,3);
    %         M = vec2rot(xhat)*M;
    %         old_err = err;
    %         iter = iter + 1;
    %         
    %     end
    %     disp(iter)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxIter = 100;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % using the second order Taylor expansion and linearization via
    % Kronecker product. This corresponds to the method proposed in the
    % paper of Ma. But it works poorly since this iterative optimization
    % won't improve the result too much.
    Asum = M .* -2;
    for i=1:maxIter      % Gauss-Newton iterations
        J = zeros(9,9);
        b = zeros(3,3);
        for k=1:size(SO3s,3)
            a1 = kron(SO3s(:,:,k)', SO3s(:,:,k)*M');
            J = J + a1;
            b = b + SO3s(:,:,k)*M'*SO3s(:,:,k);
        end
        b = b ./ size(SO3s,3) ./ 2;
        bb = vec(Asum+b+1.5.*M);
        J = J./size(SO3s,3) ./ 2 - 1.5.*kron(M',eye(3));
        o = J\bb;
        Rx = [o(1:3) o(4:6) o(7:9)];Rx = (Rx-Rx').*0.5;
        
        M = expm(Rx)*M;
        
        if norm(bb) < 1e-20
            break;
        end
    end
%     Asum = M .* 2;
%     for i=1:maxIter      % Gauss-Newton iterations
%         A1 = zeros(3,3);
%         A2 = zeros(9,9);
%         for k=1:size(SO3s,3)
%             AMA = SO3s(:,:,k)*M'*SO3s(:,:,k);
%             A1 = A1 + AMA;
%             A2 = A2 + kron(AMA',SO3s(:,:,k));
%         end
%         A3 = A1./(size(SO3s,3)*2);
%         b = vec(Asum*2 - A3 - 1.5.*M);
%         J = kron(M',1.5.*eye(3)) - A2./(size(SO3s,3)*2);
%         o = J\b;
%         Rx = [o(1:3) o(4:6) o(7:9)];Rx = (Rx-Rx').*0.5;
%         M = M*expm(Rx);
%         if norm(b) < 1e-20
%             break;
%         end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % orthogonalization
    [U, ~, V] = svd(M); 
    M = U*V';
    if det(M)<0, M = U*diag([1 1 -1])*V'; end 
end