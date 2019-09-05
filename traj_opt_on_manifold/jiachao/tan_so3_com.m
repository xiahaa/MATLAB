function x = tan_so3_com(u)

% re-compose the skew symmetric matrix x using the coordinate representation u
% basis is 
% E_1 = 1/.sqrt(2) [ 0 1 0; -1 0 0; 0 0 0];
% E_2 = 1/.sqrt(2) [ 0 0 -1; 0 0 0; 1 0 0];
% E_3 = 1/.sqrt(2) [ 0 0 0; 0 0 -1; 0 1 0];


E_1 = 1/sqrt(2)*[ 0 1 0; -1 0 0; 0 0 0];
E_2 = 1/sqrt(2)*[ 0 0 -1; 0 0 0; 1 0 0];
E_3 = 1/sqrt(2)*[ 0 0 0; 0 0 -1; 0 1 0];

x = u(1).*E_1 + u(2).*E_2 + u(3).*E_3;