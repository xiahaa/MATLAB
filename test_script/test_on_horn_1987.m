function test_on_horn_1987
    path = mfilename('fullpath');
    i = findstr(path,'\');
    folder = path(1:i(end));
    cd(folder);
    ptsrc = rand(3,1000).*10 - 5;
    yaw = 57 * pi / 180.0;
    pitch = 60 * pi / 180.0;
    roll = -30 * pi / 180.0;
    R = rot(roll,pitch,yaw);
    t = [-3;-2;5];
    ptdst = R * ptsrc;
    ptdst(1,:) = ptdst(1,:)+t(1);
    ptdst(2,:) = ptdst(2,:)+t(2);
    ptdst(3,:) = ptdst(3,:)+t(3);

    addpath('../motion_estimation/');
    R
    t
    [Ropt1,topt1] = estimate_rigid_body_transformation_Horn(ptsrc, ptdst);
    Ropt1
    topt1
    [Ropt2,topt2] = estimate_rigid_body_transformation_SVD(ptsrc, ptdst);
    Ropt2
    topt2

return

function R = q2rot(q)

    R = [q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4) 2*(q(2)*q(3)-q(1)*q(4)) 2*(q(2)*q(4)+q(1)*q(3));...
         2*(q(2)*q(3)+q(1)*q(4)) q(1)*q(1)-q(2)*q(2)+q(3)*q(3)-q(4)*q(4) 2*(q(3)*q(4)-q(1)*q(2));...
         2*(q(2)*q(4)-q(1)*q(3)) 2*(q(3)*q(4)+q(1)*q(2)) q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)];

return

function R = rot(roll, pitch, yaw)

    R1 = [cos(yaw) sin(yaw) 0;...
         -sin(yaw) cos(yaw) 0;...
         0 0 1];
    R2 = [cos(pitch) 0 -sin(pitch);0 1 0;sin(pitch) 0 cos(pitch)];

    R3 = [1 0 0;0 cos(roll) sin(roll);0 -sin(roll) cos(roll)];

    R = R3*R2*R1;

return