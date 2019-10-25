clc;clear all;close all;
% derivation of endpoint representation based trajectory generation.
syms t real
a = [1 t t^2 t^3 t^4 t^5];
da = diff(a,t);
dda = diff(da,t);
A = [subs(a,0);subs(da,0);subs(dda,0);subs(a,1);subs(da,1);subs(dda,1);];
invA = inv(A);
e = sym('e',[6,1],'real');

traj = a * invA * e;trajd = diff(traj,t);trajdd = diff(trajd,t);
[c1,d1]=coeffs(traj,e);
[c2,d2]=coeffs(trajd,e);
[c3,d3]=coeffs(trajdd,e);

EV = int(c2'*c2,t,0,1);
EA = int(c3'*c3,t,0,1);


