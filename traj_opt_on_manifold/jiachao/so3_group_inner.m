function v = so3_group_inner (a, b)

% compute the inner product of two groups of tangent vectors in SO(3)

N = size(a,3);
v = 0;
for i = 1:N
	v = v + trace(a(:,:,i)'*b(:,:,i));
end