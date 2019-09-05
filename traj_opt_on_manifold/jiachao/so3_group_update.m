function new_group = so3_group_update(old_group, d)

% update the sequence of SO(3) matrices based on the direction d

N = size(old_group, 3);
new_group = zeros(3,3,N);
for i = 1:N
	new_group(:,:,i) = old_group(:,:,i) * expm(d(:,:,i));
end
	
