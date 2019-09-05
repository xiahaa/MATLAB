function d = ro3_direction_refine(old_d, grad_f, oldrotation, newrotation, max_dis, indices)

% 2013-2-5
% The last step in 2-metric method by Bertsekas for general constrained optimization
% see also the paper by Dunn "a subspace decomposition principle for scaled..."
d = old_d;
eps = 0.0001;
N = size(oldrotation,3);
cnt = 0;
for i = 1:N
    if indices(i) ~= 0
        eta = ro3log(newrotation(:,:,i), oldrotation(:,:,i));
        dis = norm(eta,'fro');
        if abs(dis-max_dis) < eps
            inn_p = trace(grad_f(:,:,i)'*eta);
            if inn_p <= 0
                % do project on the cone T(see the Dunn's paper for more clear representation)
                innp2 = trace(old_d(:,:,i)'*(-eta));
                if innp2 > 0
                    d(:,:,i) = old_d(:,:,i) - (innp2/norm(eta,'fro')^2).*(-eta);
                    cnt = cnt +1;
                end
            end
        end
    end
end