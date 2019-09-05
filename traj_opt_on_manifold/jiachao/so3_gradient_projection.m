function d = so3_gradient_projection(oldrotation, newrotation, copyrotation, old_d, max_dis, indices)

% 2013-1-25
% Gradient projection for constrained manifold optimization
% constraint now is just round disc for every element in a manifold sequence

d = old_d;
N = size(oldrotation,3);
for i = 1:N
    if indices(i) ~= 0
        eta = ro3log(oldrotation(:,:,i), newrotation(:,:,i));
        dis = norm(eta,'fro');
        if dis > max_dis
            % gradient projection only when the solution is outside the convex set
            eta = (eta./dis).*max_dis;
            temp = ro3exp(oldrotation(:,:,i), eta);
            d(:,:,i) = real(skew_refine(ro3log(copyrotation(:,:,i), temp)));
        end
    end
end
		
