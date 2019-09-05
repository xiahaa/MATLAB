function Hess_f = so3_hess_refine(old_Hess, grad_f, oldrotation, newrotation, max_dis,indices)

% 2013-1-30
% diagonalize the Hessian matrix according to Bertsekas' 2-metric projection algorithm

% find index set to diagonalize
eps = 0.0001;
diag_set = [];
N = size(oldrotation,3);
for i = 1:N
    if indices(i) ~= 0
    	eta = ro3log(newrotation(:,:,i), oldrotation(:,:,i));
        dis = norm(eta,'fro');
        if abs(dis-max_dis) < eps
            inn_p = trace(grad_f(:,:,i)'*eta);
            if inn_p > 0
                diag_set = [diag_set, i];
            end
        end
    end
end
ind1 = (diag_set-1).*3+1;
ind2 = ind1 + 1;
ind3 = ind1 + 2;
ind_clear = reshape([ind1;ind2;ind3],[1,3*length(diag_set)]);
ind_keep = zeros(size(old_Hess));

Hess_f = old_Hess;
Hess_f(ind_clear,:) = 0;
Hess_f(:,ind_clear) = 0;
for j = 1:length(diag_set)
	temp_id = diag_set(j);
	temp_id = [3*(temp_id-1)+1, 3*(temp_id-1)+2, 3*(temp_id-1)+3];
	ind_keep(temp_id,temp_id) = 1;
end
ind_keep = (ind_keep == 1);
Hess_f(ind_keep) = old_Hess(ind_keep);
