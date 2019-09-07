function [vp, vo] = proj_vec(x, y)

% project x on the direction of y, vp is parallel to y, vo is orthognal to y, vp+ vo = x
% all the vectors are column vectors

if norm(y) == 0
	vo = x;
	vp = x-vo;
else
	vp = y.*(x'*y)./(norm(y)^2);
	vo = x - vp;
end