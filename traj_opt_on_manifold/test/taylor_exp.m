clc; close all; clear all;

a = sym('a',[3,1],'real');
b = sym('b',[3,1],'real');
c = sym('c',[3,1],'real');

syms t real

rt = a.*t^3 + b.*t^2 + c.*t;
drt = diff(rt,t);
ddrt = diff(drt,t);


normr = sqrt(rt'*rt);

alphat = ddrt - rt'*drt/normr^4*(2*cos(normr)+normr*sin(normr)-2)*(crossProd(rt,drt)) - (1-cos(normr))/normr^2*(crossProd(rt,ddrt)) ...
         + rt'*drt/normr^5*(3*sin(normr)-normr*cos(normr)-2*normr)*(crossProd(rt, crossProd(rt,drt))) ...
         + (normr-sin(normr))/normr^3*(crossProd(drt, crossProd(rt,drt))+crossProd(rt, crossProd(rt,ddrt)));


energy = alphat'*alphat;



function cp = crossProd(a1,a2)
    cp = [0 -a1(3) a1(2);a1(3) 0 -a1(1);-a1(2) a1(1) 0] * a2;
end