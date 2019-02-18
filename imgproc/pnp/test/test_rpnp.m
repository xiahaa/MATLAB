addpath ../../../MatrixLieGroup/barfoot_tro14/
XX = rand(3,2);
i1 = 1;
i2 = 2;

p1= XX(:,i1);
p2= XX(:,i2);
p0= (p1+p2)/2;
x= p2-p0; x= x/norm(x);
if abs([0 1 0]*x) < abs([0 0 1]*x)
    z= cross(x,[0; 1; 0]); z= z/norm(z);
    y= cross(z, x); y= y/norm(y);
else
    y= cross([0; 0; 1], x); y= y/norm(y);
    z= cross(x,y); z= z/norm(z);
end
Ro= [x y z]

b = acos(dot(x,[1,0,0]));
axis = cross([1,0,0]',x);axis = axis./norm(axis);
vec2rot(b.*axis)
