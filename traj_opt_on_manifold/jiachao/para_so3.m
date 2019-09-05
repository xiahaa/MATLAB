function z = para_so3(x,y,p)

% The parallel translation of y along x in the tangent space of p

temp = expm(x./2);
z = temp*y*temp;
z = expm(-x)*z;