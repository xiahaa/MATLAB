function res = sylvester_resultant_x1x2(d, theta)
%% this function calcs the sylvecter resultant of two bivariant polynomial
% it is assumed that higher order coefficient stores at first.
    if size(d,1) ~= 3 || size(theta,1)~=3
        error('Wrong input');
    end
    d12 = d(1);
    d13 = d(2);
    d23 = d(3);
    a12 = theta(1);
    a13 = theta(1);
    a23 = theta(1);
    
    c1 = cos(a12)^2 + cos(a13)^2 + cos(a23)^2 - 2*cos(a12)*cos(a13)*cos(a23) - 1;
    
16*(c1)^2
(d12^4 - 2*cos(2*a23)*d12^2*d13^2 - 2*d12^2*d23^2 + d13^4 - 2*d13^2*d23^2 + d23^4)^2

end