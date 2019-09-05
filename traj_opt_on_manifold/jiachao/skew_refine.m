function A = skew_refine(B)

% 2012-10-2
% refine the skew symmetric matrix to reduce accumulated error

x = -B(1,2);
y = B(1,3);
z = -B(2,3);
A = [0 -x y; x 0 -z; -y z 0];