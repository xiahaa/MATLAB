clc;close all;

% this file numerically evaluate the influence of right jacobian relative
% to so3
N = 1e4;
xin = 0.02:0.02:0.4;
errors = zeros(length(xin),N);
for j = 1:length(xin)
    xi = rand(3,N).*2-1;
    xi = xi./vecnorm(xi);
    xi = xi .* xin(j);
    for i = 1:N
        Ar = calcAr(xi(:,i));
        errors(j,i) = norm(Ar-eye(3),'fro');
    end
end

figure
cmap=lines(4);
plot(xin,mean(errors,2),'Color',cmap(1,:),'LineWidth',1.5);hold on;grid on;
plot(xin,median(errors,2),'Color',cmap(2,:),'LineWidth',1.5);
plot(xin,max(errors,[],2),'Color',cmap(3,:),'LineWidth',1.5);
plot(xin,min(errors,[],2),'Color',cmap(4,:),'LineWidth',1.5);
legend({'mean','median','max','min'},'FontName','Arial','FontSize',15);


function Ar = calcAr(r)
    theta = norm(r);
    rskew = skew(r);
    if theta < 1e-6
        Ar = eye(3);
    else
        Ar = eye(3) - (1-cos(theta))/theta^2*rskew + (theta-sin(theta))/theta^3*rskew*rskew;
    end
end
