function u = tan_so3_dec(x)

% decompose the skew symmetric matrix x in to coordinate representation u
% basis is 
% E_1 = 1/.sqrt(2) [ 0 1 0; -1 0 0; 0 0 0];
% E_2 = 1/.sqrt(2) [ 0 0 -1; 0 0 0; 1 0 0];
% E_3 = 1/.sqrt(2) [ 0 0 0; 0 0 -1; 0 1 0];

u1 = sqrt(2)*x(1,2);
u2 = sqrt(2)*x(3,1);
u3 = sqrt(3)*x(3,2);

u = [u1; u2; u3];