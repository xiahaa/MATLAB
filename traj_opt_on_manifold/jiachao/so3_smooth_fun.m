function v = so3_smooth_fun(new_seq, old_seq, coeff, indices)

% compute the value of the cost function
% the cost function is the distance between new_seq and old_seq, and also the smoothness penalty on new_seq

v = 0;
N = size(new_seq, 3);
for i = 1:N
	v = v + 0.5 * norm(logm(old_seq(:,:,i)' * new_seq(:,:,i)), 'fro')^2 .* indices(i);
	if i < N
		v = v + 0.5*coeff*norm(logm(new_seq(:,:,i)'*new_seq(:,:,i+1)), 'fro')^2;
	end
end
	