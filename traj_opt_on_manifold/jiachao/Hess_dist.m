function H = Hess_dist(p,x)

% Compute the Hessian matrix of function f(x) = 1/2*||dist(p,x)||^2
% the distance if geodesic distance
% p and x are on SO(3) manifold
% the H is represented on the canonical coordinate on the vector space of x
% the orthonormal basis is 
% E_1 = 1/.sqrt(2) [ 0 1 0; -1 0 0; 0 0 0];
% E_2 = 1/.sqrt(2) [ 0 0 -1; 0 0 0; 1 0 0];
% E_3 = 1/.sqrt(2) [ 0 0 0; 0 0 -1; 0 1 0];


% constant sectional curvature
lambda = 1/4;

tan_vec = logm(p'*x);
tan_vec = real(tan_vec);
r = norm(tan_vec,'fro');

% parallel translation of log_p^x
para_log = para_so3(tan_vec, tan_vec, p);

% decompose the para_log into coordinate representation
para_log_u = tan_so3_dec(para_log);

H = zeros(3,3);
const_num = r*sqrt(lambda)/(tan(sqrt(lambda)*r)+1e-8);

for i = 1:3
	for j = i:3
		e1 = zeros(3,1);
		e2 = zeros(3,1);
		e1(i) = 1;
		e2(j) = 1;
		[e1p,e1o] = proj_vec(e1, para_log_u);
		[e2p,e2o] = proj_vec(e2, para_log_u);
		
		H(i,j) = e1p'*e2p + const_num*(e1o'*e2o);
	end
end
H(2,1) = H(1,2);
H(3,1) = H(1,3);
H(3,2) = H(2,3);
