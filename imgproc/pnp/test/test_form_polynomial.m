clc;clear all;close all;
syms x y z real
syms a1 a2 a3 a4 real
syms b1 b2 b3 b4 real
syms c1 c2 c3 c4 real

p1 = a1*x*x+a2*z*z+a3*x*z+a4; 
p2 = b1*y*y+b2*z*z+b3*y*z+b4;
p3 = resultant(p1,p2,z);
p4 = c1*x^2 + c3*x*y + c2*y^2 + c4;
p5 = resultant(p3,p4,y);


cc1 = rand(1,4);
cc2 = rand(1,4);
cc3 = rand(1,5);
a1 = cc1(1);a2 = cc1(2);a3 = cc1(3);a4 = cc1(4);
b1 = cc2(1);b2 = cc2(2);b3 = cc2(3);b4 = cc2(4);
c1 = cc3(1);c2 = cc3(2);c3 = cc3(3);c4 = cc3(4);

[c,t] = coeffs(p5,x);
double(subs(c))

addpath ../utils/

poly12_4order = sylvester_resultant_22([a1 a2 a3 a4], [b1 b2 b3 b4]);
poly1_8order = sylvester_resultant_42(poly12_4order, [c1 c2 c3 c4]);
double(subs(c))-poly1_8order

%% done



 