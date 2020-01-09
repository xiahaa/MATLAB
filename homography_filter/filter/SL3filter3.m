function [Hest, tau] = SL3filter3(Hest, tau, dt, pref, pcur,omega)
% Hamel T, Mahony R, Trumpf J, et al.
% Homography estimation on the special linear group based on direct point correspondence[C]
% //2011 50th IEEE Conference on Decision and Control and European Control Conference. IEEE, 2011: 7902-7908.
% NOTE: do not know how to validate the observor.
    % compute point correction
    w = correctionByPoints(Hest, pref, pcur);
    % update tau
    KI = 1; KP = 4;

    Omega = toso3(omega);
    taudot = tau*Omega+KI*adjointSL3(Hest',w);
    Hdot = (Omega+tau)-1/3*trace(tau)*eye(3)+inv(Hest)*KP*w*Hest;
    tau = tau + taudot*dt;
    Hest = Hest*expm((Hdot)*dt);
end