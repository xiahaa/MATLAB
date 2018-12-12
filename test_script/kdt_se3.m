clc;
close all;
clear all;

% x = sym('x',[6,1],'real');
% y = sym('y',[6,1],'real');
N = 1000;
dt = zeros(N,1);
dtt = zeros(N,1);
dt1 = zeros(N,1);
dt2 = zeros(N,1);
dt3 = zeros(N,1);
dt4 = zeros(N,1);
dt5 = zeros(N,1);
for i = 1:N
    x = rand(6,1);
    y = rand(6,1);

    X = lie_hat(-x);
    Y = lie_hat(y);

    Z1 = X+Y;
    Z2 = Z1 + 0.5.*lie_bracket(X,Y);
    Z3 = Z2 + 1/12.*(lie_bracket(X,lie_bracket(X,Y))+lie_bracket(Y,lie_bracket(Y,X)));
    Z4 = Z3 - 1/24.*(lie_bracket(Y,lie_bracket(X,lie_bracket(X,Y))));
    Z5 = Z4 - 1/720.*(lie_bracket(Y,lie_bracket(Y,lie_bracket(Y,lie_bracket(Y,X))))+lie_bracket(X,lie_bracket(X,lie_bracket(X,lie_bracket(X,Y))))) ...
            + 1/360.*(lie_bracket(X,lie_bracket(Y,lie_bracket(Y,lie_bracket(Y,X))))+lie_bracket(Y,lie_bracket(X,lie_bracket(X,lie_bracket(X,Y))))) ...
            + 1/120.*(lie_bracket(Y,lie_bracket(X,lie_bracket(Y,lie_bracket(X,Y))))+lie_bracket(X,lie_bracket(Y,lie_bracket(X,lie_bracket(Y,X)))));

%     disp(lie_bracket(X,Y)-lie_bracket_hat(-x,y));

    d1 = lie_vee(Z1);
    d2 = lie_vee(Z2);
    d3 = lie_vee(Z3);
    d4 = lie_vee(Z4);
    d5 = lie_vee(Z5);


    dist1 = d1' * d1;
    dist2 = d2' * d2;
    dist3 = d3' * d3;
    dist4 = d4' * d4;
    dist5 = d5' * d5;
    
    dt1(i) = dist1;
    dt2(i) = dist2;
    dt3(i) = dist3;
    dt4(i) = dist4;
    dt5(i) = dist5;
    dt(i) = norm(lie_vee(logm(inv(expm(lie_hat(x)))*expm(Y))))^2;
    dtt(i) = norm(lie_vee(logm(expm(X)*expm(Y))))^2;
end
% figure(1)
% plot(dt,'r-o');hold on;grid on
% plot(dtt,'b-');

figure(2)
subplot(1,2,1)
plot(abs(dt - dt1),'b-');hold on;grid on
plot(abs(dt - dt2),'g-');hold on;grid on
subplot(1,2,2)
plot(abs(dt - dt3),'c-');hold on;grid on
plot(abs(dt - dt4),'m-');hold on;grid on
plot(abs(dt - dt5),'k-');hold on;grid on


function X = lie_hat(x)
    if numel(x) == 3
        X = [0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0];
    elseif numel(x) == 6
        X = [0      -x(3)	x(2)	x(4)
             x(3)    0       -x(1)   x(5)
            -x(2)   x(1)    0       x(6)
             0       0       0       0]; 
    else
        error('no hat operation!');
    end
end

function x = lie_vee(X)
    if size(X,1) == 3
        x = [X(3,2);X(1,3);-X(1,2)];
    elseif size(X,1) == 4
        x = [X(3,2);X(1,3);-X(1,2);X(1,4);X(2,4);X(3,4)];
    else
        error('no vee operation!');
    end
end

function XY = lie_bracket(X,Y)
    XY = X*Y-Y*X;
end

function XY = lie_bracket_hat(x,y)
    if numel(x) == 3
        xy = cross(x,y);
    elseif numel(x) == 6
        xy = [cross(x(1:3),y(1:3));cross(x(1:3),y(4:6))-cross(y(1:3),x(4:6))];
    end
    XY = lie_hat(xy);
end