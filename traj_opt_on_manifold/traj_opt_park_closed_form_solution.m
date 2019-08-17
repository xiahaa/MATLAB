% this file implements the trajectory generation method on SO3 by using the
% method proposed in Smooth Invariant interpolation of rotations.

clc; clear all; close all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% add Lie library
addpath ../MatrixLieGroup/barfoot_tro14/
addpath ../beautiful_plot/
%
% r0 = [0.2 0.1 0.1]';
% r1 = [0.6 0.4 0.4]';
% w0 = [0.5 0.1 0.1]';
% w1 = [0 0 0]';
% 
% R0 = vec2rot(r0);
% R1 = vec2rot(r1);
% 
% dr1 = rot2vec(R0'*R1);
% 
% % polynomial
% A = [1 0 0 1 0 0 1 0 0; ...
%      0 1 0 0 1 0 0 1 0; ...
%      0 0 1 0 0 1 0 0 1; ...
%      0 0 0 0 0 0 1 0 0; ...
%      0 0 0 0 0 0 0 1 0; ...
%      0 0 0 0 0 0 0 0 1; ...
%      3 0 0 2 0 0 1 0 0; ...
%      0 3 0 0 2 0 0 1 0; ...
%      0 0 3 0 0 2 0 0 1];
%  p = A\[dr1;w0;calcAr(dr1)\w1];
% 
%  %
%  a = p(1:3);b = p(4:6);c = p(7:9);
%  
%  t = 0:0.2:1;
%  
%  figure;
%  
%  T = [R0 zeros(3,1);[0 0 0 1]];
%  DrawAxis(T, 1, 'r', 'r', 'r');
%  hold on;
%  view(3);
%  grid on;
%  T = [R1 zeros(3,1);[0 0 0 1]];
%  DrawAxis(T, 1, 'k', 'k', 'k');
%  
%  for i = 1:length(t)
%      rs = a.*t(i)^3 + b.*t(i)^2 + c.*t(i);
%      Rs = R0*vec2rot(rs);
%      Ts = [Rs zeros(3,1);[0 0 0 1]];
%      DrawAxis(Ts, 1, 'r', 'g', 'b');
%      pause(0.1)
%  end

% here, brute-force evaluation of cost variation
N = 1000;
cost1 = zeros(1,N);
cost2 = zeros(1,N);

b = []


function Ar = calcAr(r)
    theta = norm(r);
    rskew = skew(r);
    Ar = eye(3) - (1-cos(theta))/theta^2*rskew + (theta-sin(theta))/theta^3*rskew*rskew;
end

function cost = costFunc2(a,b,c)
    a1 = a(1);
    a2 = a(2);
    a3 = a(3);
    b1 = b(1);
    b2 = b(2);
    b3 = b(3);
    c=c;% dumy
    cost = (12*a1^2 + 12*a2^2 + 12*a3^2) + (12*a1*b1 + 12*a2*b2 + 12*a3*b3) + (4*b1^2 + 4*b2^2 + 4*b3^2);
end

function cost = costFunc1(a,b,c)
    cost = integral(@intf,0,1);
    function y = intf(t)
        rt = a.*t^3 + b.*t^2 + c.*t;
        drt = 3.*a.*t^2 + 2.*b.*t + c;
        ddrt = 6.*a.*t + 2.*b;

        normr = sqrt(rt'*rt);

        alphat = ddrt - rt'*drt/normr^4*(2*cos(normr)+normr*sin(normr)-2)*(crossProd(rt,drt)) - (1-cos(normr))/normr^2*(crossProd(rt,ddrt)) ...
                 + rt'*drt/normr^5*(3*sin(normr)-normr*cos(normr)-2*normr)*(crossProd(rt, crossProd(rt,drt))) ...
                 + (normr-sin(normr))/normr^3*(crossProd(drt, crossProd(rt,drt))+crossProd(rt, crossProd(rt,ddrt)));
             
        y = alphat'*alphat;
    end
end
