function [Hest, Aest] = SL3filter1(Hest, Hmeas, Aest, dt)
% Malis E, Hamel T, Mahony R, et al. Dynamic estimation of homography
% transformations on the special linear group for visual servo control[C]
% 2009 IEEE international conference on robotics and automation.
% IEEE, 2009: 1498-1503.
    Hpred = Hest;
    Htilde = inv(Hpred)*Hmeas;
    Hproj = tosl3(Htilde'*(eye(3) - Htilde));

    kh = 2; ka = 1;

    alpha = -kh*adjointSL3(Htilde,Hproj);
    beta  = -ka*Hproj;

    Adot = beta;

    Aest = Aest + Adot*dt;
    Hdot = (adjointSL3((Htilde),Aest) + alpha);
    Hest = Hest*expm(Hdot*dt);% alternatively, Hest = Hest*expm(Hdot*dt)
    Hest = Hest ./ (det(Hest)^(1/3));

end