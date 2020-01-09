
% the benefit of this observor is that no direct measuremnt of H is
% necessary
function Hest = SL3filter2(Hest, q, p, A, dt)
% Hamel T, Mahony R, Trumpf J, et al.
% Homography estimation on the special linear group based on direct point correspondence[C]
% //2011 50th IEEE Conference on Decision and Control and European Control Conference. IEEE, 2011: 7902-7908.

    % compute point correction
    w = correctionByPoints(Hest, p, q);
    % update tau
    KP = 0.5;

    Hdot = A+inv(Hest)*KP*w*Hest;
    Hest = Hest*(eye(3)+(Hdot)*dt);
end