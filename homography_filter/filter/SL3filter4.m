function [Hest, tau] = SL3filter4(Hest, tau, dt, lref, lcur,omega)
% Hua, Minh-Duc, et al. 
% "Point and line feature-based observer design on SL (3) for Homography estimation and its application to image stabilization." (2017).
    % compute point correction
    w = correctionByLines(Hest, lref, lcur);
    w = -w;
    % update tau
    KI = 1; KP = 4;

    Omega = toso3(omega);
    taudot = tau*Omega+KI*adjointSL3(Hest',w);
    Hdot = (Omega+tau)-1/3*trace(tau)*eye(3)+inv(Hest)*KP*w*Hest;
    tau = tau + taudot*dt;
    Hest = Hest*expm((Hdot)*dt);
end