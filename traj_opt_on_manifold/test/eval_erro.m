clc;close all;clear all;
M = 1e3;
ts = [0.001, 0.01, 0.1, 0.2, 0.4, 0.8, 1];
for i = 1:M
    tic
    phi1 = rand(3,1);
    R1 = expSO3(phi1);
    N = 1e3;
    ranphi = rand(3,N);
    ranphi = ranphi ./ vecnorm(ranphi);
    
    for j = 1:length(ts)
        ranphi_t = ranphi .* ts(j);
        ranphi_t = ranphi_t + phi1;
        for k = 1:N
            R2 = expSO3(ranphi_t(:,k));
            err1(i,j,k) = norm(logSO3(R1'*R2));
            err2(i,j,k) = sqrt(2)*norm(ranphi_t(:,k)-phi1);
        end
    end
    toc
end
save('tmp.mat','err1','err2');

% load('tmp.mat');

err1 = permute(err1,[2,1,3]);
err2 = permute(err2,[2,1,3]);

err1 = reshape(err1,length(ts),[])';
err2 = reshape(err2,length(ts),[])';

err = err2 - err1;

figure
boxplot(err,'Notch','off','Labels',convertStringsToChars(string(ts)),'Whisker',3.5);
title('Compare Random Data from Different Distributions')



